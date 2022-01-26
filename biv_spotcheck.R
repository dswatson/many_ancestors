### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(ranger)
library(ppcor)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(1234, kind = "L'Ecuyer-CMRG")

################################################################################

### SIMULATION ###

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param sp Sparsity of the connections from background to foreground.
#' @param r2 Proportion of variance explained by endogenous features.
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param g Expected output of an independence oracle, one of \code{"xy"},
#'   \code{"ci"}, or \code{"na"}.
#' 

# Data simulation function
sim_dat <- function(n, d_z, rho, sp, r2, lin_pr, g) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  var_z <- 1 / d_z # Does this make a difference?
  Sigma <- toeplitz(rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Random Rademacher weights
  beta <- gamma <- double(length = d_z)
  k <- round((1 - sp) * d_z)
  beta[sample(d_z, k)] <- sample(c(1, -1), k, replace = TRUE)
  gamma[sample(d_z, k)] <- sample(c(1, -1), k, replace = TRUE)
  # Nonlinear transformations?
  if (lin_pr < 1) {
    z_out <- zz <- z
    # Create matrix zz of random nonlinear transformations
    idx <- split(sample.int(d_z), sort(seq_len(d_z) %% 4))
    names(idx) <- c('sq', 'sqrt', 'sftpls', 'relu')
    zz[, idx$sq] <- z[, idx$sq]^2
    zz[, idx$sqrt] <- sqrt(abs(z[, idx$sqrt]))
    zz[, idx$sftpls] <- log(1 + exp(z[, idx$sftpls]))
    zz[, idx$relu] <- ifelse(z[, idx$relu] > 0, z[, idx$relu], 0)
    # Sample columns from zz with probability 1 - lin_pr
    nonlin_idx <- sample.int(d_z, d_z * (1 - lin_pr))
    z_out[, nonlin_idx] <- zz[, nonlin_idx]
    signal_x <- as.numeric(z_out %*% beta)
    signal_y <- as.numeric(z_out %*% gamma)
  } else {
    signal_x <- as.numeric(z %*% beta)
    signal_y <- as.numeric(z %*% gamma)
  }
  # Generate noise
  sim_noise <- function(signal, r2) {
    var_mu <- var(signal)
    var_noise <- (var_mu - r2 * var_mu) / r2
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  # X data
  x <- signal_x + sim_noise(signal_x, r2)
  # Identifiable?
  if (g == 'na') {
    shared_parents <- which(beta != 0 & gamma != 0)
    u_idx <- sample(shared_parents, size = length(shared_parents)/2)
    z <- z[, -u_idx]
    d_z <- ncol(z)
    d_u <- length(u_idx)
  } else {
    d_u <- 0
  }
  # Y data
  if (g %in% c('ci', 'na')) {
    y <- signal_y + sim_noise(signal_y, r2)
  } else if (g == 'xy') {
    signal_z_to_y <- signal_y
    xzr <- 1 / (k + 1)
    sigma_xy <- sqrt(xzr * var(signal_z_to_y))
    gamma_x <- sigma_xy / sd(x)
    signal_y <- signal_z_to_y + x * gamma_x
    y <- signal_y + sim_noise(signal_y, r2)
  }
  # Export
  params <- list(
    'n' = n, 'd_z' = d_z, 'd_u' = d_u, 'rho' = rho, 'sp' = sp, 'r2' = r2, 
    'lin_pr' = lin_pr, 'g' = g
  )
  out <- list(
    'dat' = data.table(z, 'x' = x, 'y' = y),
    'wts' = list('beta' = beta, 'gamma' = gamma), 'params' = params
  )
  return(out)
}

################################################################################

### CONFOUNDER BLANKET REGRESSION ###

#' @param m Number of nested models to fit.
#' @param max_d Number of predictors in largest model.
#' @param min_d Number of predictors in smallest model.
#' @param decay Exponential decay parameter.

# Precompute subset sizes for RFE
subsets <- function(m, max_d, min_d, decay) {
  out <- round(min_d + ((max_d - min_d) / (m + 1)^decay) * seq_len(m + 1)^decay)
  out <- na.omit(unique(out)[seq_len(m)])
  return(out)
}


#' @param x Design matrix.
#' @param y Outcome vector.
#' @param trn Training indices.
#' @param tst Test indices.
#' @param f Regression method, either \code{"lasso"} or \code{"rf"}.
#' 

# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, trn, tst, f) {
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
    betas <- coef(fit, s = fit$lambda)[-1, ]
  } else if (f == 'rf') {
    fit <- ranger(x = x[trn, ], y = y[trn], importance = 'impurity',
                  num.trees = 200, num.threads = 1)
    yhat_f0 <- predict(fit, data = x[tst, ], 
                       num.trees = 50, num.threads = 1)$predictions
    vimp <- data.frame('feature' = colnames(x), 
                       'imp' = fit$variable.importance) %>%
      arrange(desc(imp))
    s <- subsets(m = 10, max_d = ncol(x), min_d = 5, decay = 2)
    y_hat <- sapply(seq_along(s), function(k) {
      tmp_x <- x[trn, vimp$feature[seq_len(s[k])]]
      tmp_f <- ranger(x = tmp_x, y = y[trn], num.trees = 50, num.threads = 1)
      predict(tmp_f, data = x[tst, ], num.threads = 1)$predictions
    })
    y_hat <- cbind(y_hat, yhat_f0)
    beta <- double(length = ncol(x))
    names(beta) <- colnames(x)
    betas <- sapply(seq_along(s), function(k) {
      out <- beta
      keep <- vimp$feature[seq_len(s[k])]
      out[keep] <- 1
      return(out)
    })
    betas <- cbind(betas, rep(1, ncol(x)))
  }
  epsilon <- y_hat - y[tst]
  mse <- colMeans(epsilon^2)
  betas <- betas[, which.min(mse)]
  out <- ifelse(betas == 0, 0, 1)
  return(out)
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param B Number of complementary pairs to draw for stability selection.
#' 

# Compute (de)activation rates for X -> Y and Y -> X
rate_fn <- function(sim_obj, B) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  zx <- cbind(z, x)
  y <- dat$y
  zy <- cbind(z, y)
  # Linear or nonlinear?
  f <- ifelse(sim_obj$params$lin_pr == 1, 'lasso', 'rf')
  # Compute disconnections and (de)activations per subsample
  fit_fn <- function(b) {
    # Take complementary subsets
    i_set <- sample(n, round(0.5 * n))
    i_trn <- sample(i_set, round(0.8 * length(i_set)))
    i_tst <- setdiff(i_set, i_trn)
    j_set <- seq_len(n)[-i_set]
    j_trn <- sample(j_set, round(0.8 * length(j_set)))
    j_tst <- setdiff(j_set, j_trn)
    # Compute active sets
    s <- data.frame(
      y0 = c(l0(z, y, i_trn, i_tst, f), NA_real_, 
             l0(z, y, j_trn, j_tst, f), NA_real_), 
      y1 = c(l0(zx, y, i_trn, i_tst, f), l0(zx, y, j_trn, j_tst, f)),
      x0 = c(l0(z, x, i_trn, i_tst, f), NA_real_, 
             l0(z, x, j_trn, j_tst, f), NA_real_), 
      x1 = c(l0(zy, x, i_trn, i_tst, f), l0(zy, x, j_trn, j_tst, f))
    )
    # Record disconnections and (de)activations
    dis_i <- any(c(s$y1[d_z + 1], s$x1[d_z + 1]) == 0)
    dis_j <- any(c(s$y1[2 * (d_z + 1)], s$x1[2 * (d_z + 1)]) == 0)
    dis <- rep(c(dis_i, dis_j), each = d_z + 1)
    d_xy <- s$y0 == 1 & s$y1 == 0
    a_xy <- s$x0 == 0 & s$x1 == 1
    d_yx <- s$x0 == 1 & s$x1 == 0
    a_yx <- s$y0 == 0 & s$y1 == 1
    extras <- c(d_z + 1, 2 * (d_z + 1))
    d_xy[extras] <- a_xy[extras] <- d_yx[extras] <- a_yx[extras] <- NA_real_
    # Export
    out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_z + 1), 
                      dis, d_xy, a_xy, d_yx, a_yx,
                      z = rep(seq_len(d_z + 1), times = 2))
    return(out)
  }
  out <- foreach(bb = seq_len(B), .combine = rbind) %do% 
    fit_fn(bb)
  # Compute rates
  out[, disr := sum(dis) / .N]
  out[, drxy := sum(d_xy) / .N, by = z]
  out[, arxy := sum(a_xy) / .N, by = z]
  out[, dryx := sum(d_yx) / .N, by = z]
  out[, aryx := sum(a_yx) / .N, by = z]
  # Tidy up, export
  out <- unique(out[, .(z, disr, drxy, arxy, dryx, aryx)])
  return(out)
}


#' @param res Results object output by \code{rate_fn}.
#' @param B Number of complementary pairs to draw for stability selection.
#' 

# Compute consistency lower bound
lb_fn <- function(res, B) {
  # Subset the data
  df <- na.omit(res[, .(z, drxy, arxy, dryx, aryx)])
  # Loop through thresholds
  lies <- function(tau) {
    # Internal consistency
    df[, dxy := ifelse(drxy >= tau, 1, 0)]
    df[, axy := ifelse(arxy >= tau, 1, 0)]
    df[, dyx := ifelse(dryx >= tau, 1, 0)]
    df[, ayx := ifelse(aryx >= tau, 1, 0)]
    df[, int_err := ifelse((dxy + axy > 1) | (dyx + ayx > 1), 1, 0)]
    int_err <- sum(df$int_err)
    # External consistency
    sum_xy <- df[, sum(dxy + axy)]
    sum_yx <- df[, sum(dyx + ayx)]
    ext_err <- ifelse(min(c(sum_xy, sum_yx)) > 0, 1, 0)
    # Export
    out <- data.table('tau' = tau, 'int_err' = int_err, 'ext_err' = ext_err)
  }
  lie_df <- foreach(tt = seq_len(2 * B) / (2 * B), .combine = rbind) %do% 
    lies(tt)
  # Compute minimal thresholds
  min_int <- lie_df[int_err == 0, min(tau)]
  min_ext <- lie_df[ext_err == 0, min(tau)]
  min_two <- lie_df[int_err == 0 & ext_err == 0, min(tau)] # It's always ext
  return(min_two)
}


#' @param res Results object output by \code{rate_fn}.
#' @param lb Lower bound output by \code{lb_fn}.
#' @param order Assume X \preceq Y or Y \preceq X?
#' @param rule Detect via deactivation (\code{"R1"}) or activation (\code{"R2"})?
#' @param B Number of complementary pairs to draw for stability selection.
#' 

# Infer causal direction using stability selection
ss_fn <- function(res, lb, order, rule, B) {
  # Subset the data
  if (order == 'xy' & rule == 'R1') {
    r <- res$drxy 
  } else if (order == 'xy' & rule == 'R2') {
    r <- res$arxy 
  } else if (order == 'yx' & rule == 'R1') {
    r <- res$dryx 
  } else if (order == 'yx' & rule == 'R2') {
    r <- res$aryx 
  }
  r <- na.omit(r)
  if (max(r) == 0) {
    dat <- data.frame(surplus = 0)
  } else {
    # Stability selection parameters
    theta <- mean(r)
    ub <- minD(theta, B) * sum(r <= theta)
    tau <- seq_len(2 * B) / (2 * B)
    # Do any features exceed the upper bound?
    dat <- data.frame(tau, err_bound = ub) %>%
      filter(tau > lb) %>%
      rowwise() %>%
      mutate(detected = sum(r >= tau)) %>% 
      ungroup() %>%
      mutate(surplus = ifelse(detected > err_bound, 1, 0))
  }
  # Export
  out <- data.table(
    'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param gamma Omission threshold.

# Wrap it up
cbr_fn <- function(sim_obj, gamma) {
  # Compute rates for each z
  res <- rate_fn(sim_obj, B = 50)
  # Disconnected?
  if (res$disr[1] > gamma) {
    decision <- 'ci'
  } else {
    # (De)activation rates
    lb <- lb_fn(res, B = 50)
    sum_tbl <- foreach(oo = c('xy', 'yx'), .combine = rbind) %:%
      foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
      ss_fn(res, lb, oo, rr, B = 50)
    if (sum_tbl[order == 'xy', sum(decision)] > 0) {
      decision <- 'xy'
    } else if (sum_tbl[order == 'yx', sum(decision)] > 0) {
      decision <- 'yx'
    } else {
      decision <- 'na'
    }
  }
  # Export
  out <- data.table(method = 'cbr', h = decision) 
  return(out)
}

################################################################################

### ENTNER METHOD ###

#' @param x First vector.
#' @param y Second vector.
#' @param z Conditioning set.
#' @param trn Training index.
#' @param tst Test index.

# GCM subroutine (Shah & Peters 2020). Tests conditional independence of 
# x and y given z using random forest regression.
gcm_test <- function(x, y, z, trn, tst) {
  rf1 <- ranger(x = z[trn, ], y = x[trn], num.trees = 50, num.threads = 1)
  rf2 <- ranger(x = z[trn, ], y = y[trn], num.trees = 50, num.threads = 1)
  eps1 <- x[tst] - predict(rf1, z[tst, ], num.threads = 1)$predictions
  eps2 <- y[tst] - predict(rf2, z[tst, ], num.threads = 1)$predictions
  nn <- length(tst)
  R <- eps1 * eps2
  R.sq <- R^2
  meanR <- mean(R)
  z_score <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p.value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p.value)
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param alpha Significance threshold for inferring dependence.
#' @param tau Other threshold for "inferring" independence.
#' @param B Number of random subset-variable pairs for testing.

# Constraint-based: for this comparison, we presume X \preceq Y
# and assume access to the true data sparsity
constr_fn <- function(sim_obj, alpha, tau, B) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  y <- dat$y
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  # Take subsets of size equal to expected size of admissible set
  fctr <- d_z + sim_obj$params$d_u
  k <- round(sim_obj$params$sp * fctr)/2
  # Evaluate R1 10x more frequently than R2
  r2_idx <- seq_len(B/10)
  # Entner function
  entner <- function(b) {
    # Sample a random subset Z_b and variable W
    z_idx <- sample(d_z, k)
    z_b <- z[, z_idx]
    j <- sample(seq_len(d_z)[-z_idx], 1)
    w <- z[, j]
    if (linear == TRUE) {
      ### RULE 1 ###
      pmat0 <- suppressWarnings(pcor(cbind(y, w, z_b))$p.value)
      p1.i <- pmat0[1, 2]
      pmat1 <- suppressWarnings(pcor(cbind(y, w, z_b, x))$p.value)
      p1.ii <- pmat1[1, 2]
      ### RULE 2 ###
      if (b %in% r2_idx) {
        pmat2 <- suppressWarnings(pcor(cbind(y, x, z_b))$p.value)
        p2.i <- pmat2[1, 2]
        pmat3 <- suppressWarnings(pcor(cbind(x, w, z_b))$p.value)
        p2.ii <- pmat3[1, 2]
      }
    } else {
      ### RULE 1 ###
      trn <- sample(n, round(0.8 * n))
      tst <- seq_len(n)[-trn]
      p1.i <- gcm_test(y, w, z_b, trn, tst)
      p1.ii <- gcm_test(y, w, cbind(z_b, x), trn, tst)
      ### RULE 2 ###
      if (b %in% r2_idx) {
        p2.i <- gcm_test(y, x, z_b, trn, tst)
        p2.ii <- gcm_test(x, w, z_b, trn, tst)
      } 
    }
    # Apply rules
    if (!b %in% r2_idx) {
      p2.i <- 0
      p2.ii <- 1
    }
    r1 <- ifelse(p1.i <= alpha & p1.ii >= tau, 1, 0)
    r2 <- ifelse(p2.i >= tau | (p2.ii <= alpha & p1.i >= tau), 1, 0)
    # Export
    out <- data.table(r1, r2)
    return(out)
  }
  # Apply Entner's rules with B random subset-variable pairs
  df <- foreach(bb = seq_len(B), .combine = rbind) %do%
    entner(bb)
  # Selecting different decision thresholds based on experimentation
  # Note -- this is very generous!
  if (df[, sum(r1) / .N] >= 1/250) {
    h <- 'xy' 
  } else if (df[seq_len(B/10), sum(r2) / .N] >= 1/5) {
    h <- 'ci'
  } else {
    h <- 'na'
  }
  out <- data.table(method = 'constr', h)
  return(out)
}

################################################################################

### SCORE-BASED METHOD ###

#' @param x Design matrix.
#' @param y Response vector.
#' @param trn Training index.
#' @param tst Test index.
#' @param f Function class.

# Perform feature selection on a training set, return a list with 
# (i) R^2 and (ii) test residuals
scr_fn <- function(x, y, trn, tst, f) {
  # Fit models
  if (f == 'lasso') {
    fit <- cv.glmnet(x[trn, ], y[trn])
    y_hat <- as.numeric(predict(fit, newx = x[tst, ], s = 'lambda.min'))
  } else if (f == 'rf') {
    fit <- ranger(x = x[trn, ], y = y[trn], importance = 'impurity',
                  num.trees = 200, num.threads = 1)
    vimp <- data.frame('feature' = colnames(x), 
                       'imp' = fit$variable.importance) %>%
      arrange(desc(imp))
    s <- subsets(m = 10, max_d = ncol(x), min_d = 5, decay = 2)
    rf_list <- lapply(seq_along(s), function(k) {
      tmp_x <- x[trn, vimp$feature[seq_len(s[k])]]
      tmp_f <- ranger(x = tmp_x, y = y[trn], num.trees = 50, num.threads = 1)
      list('mse' = tmp_f$prediction.error, 'f' = tmp_f)
    })
    oob <- sapply(seq_along(rf_list), function(k) rf_list[[k]]$mse)
    y_hat <- predict(rf_list[[which.min(oob)]]$f, data = x[tst, ], 
                     num.threads = 1)$predictions
  }
  # Compute residuals, r2
  eps <- y[tst] - y_hat
  r2 <- 1 - (sum(eps^2) / sum((y[tst] - mean(y[tst]))^2))
  return(list('eps' = eps, 'r2' = r2))
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param alpha Significance threshold for testing.

# Evaluate graph structures using score-based method
score_fn <- function(sim_obj, alpha) {
  # Get data, parameters
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  y <- dat$y
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  f <- ifelse(linear, 'lasso', 'rf')
  # Train/test split
  trn <- sample(n, size = round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Score X -> Y
  scr_xy <- scr_fn(cbind(z, x), y, trn, tst, f)
  # Score Y -> X
  scr_yx <- scr_fn(cbind(z, y), x, trn, tst, f)
  # Score X ~ Y
  scr_x <- scr_fn(z, x, trn, tst, f)
  scr_y <- scr_fn(z, y, trn, tst, f)
  # Decision procedure
  xy_scr <- scr_x$r2 + scr_xy$r2
  yx_scr <- scr_y$r2 + scr_yx$r2
  ci_scr <- scr_x$r2 + scr_y$r2
  if (ci_scr == max(c(xy_scr, yx_scr, ci_scr))) {
    h <- 'ci'
  } else {
    if (xy_scr > yx_scr) {
      p_value <- cor.test(x[tst], scr_xy$eps)$p.value
      h <- ifelse(p_value <= alpha, 'na', 'xy')
    } else {
      p_value <- cor.test(y[tst], scr_yx$eps)$p.value
      h <- ifelse(p_value <= alpha, 'na', 'yx')
    }
  }
  # Export
  out <- data.table(method = 'score', h)
  return(out)
}

################################################################################

### PUT IT ALL TOGETHER ###


# Spot check: fix constraint-based and CBR plz
sims <- expand.grid(n = c(1000, 2000, 4000), g = c('xy', 'ci', 'na')) %>%
  mutate(sp = 0.5, d_z = 100, rho = 0.25, r2 = 2/3, lin_pr = 1/5, 
         s_id = row_number()) %>%
  as.data.table(.)
sims <- sims[n > 1000]
# Initialize
out <- data.table(
  method = NA, h = NA, s_id = NA, idx = NA
)
res_file <- './results/nl_biv_spotcheck.rds'
saveRDS(out, res_file)
# Big ol' wrapper
big_loop <- function(sims_df, sim_id, i) {
  # Housekeeping
  sdf <- sims_df[s_id == sim_id]
  linear <- FALSE
  n_reps <- 500
  # Simulate data
  sim_obj <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, sp = sdf$sp, 
                     r2 = sdf$r2, lin_pr = sdf$lin_pr, g = sdf$g)
  # Confounder blanket regression
  df_b <- cbr_fn(sim_obj, gamma = 0.75) 
  # Constraint function
  df_c <- constr_fn(sim_obj, alpha = 0.1, tau = 0.5, B = n_reps)
  # Import, export
  old <- readRDS(res_file)
  new <- rbind(df_b, df_c) %>%
    mutate(s_id = sim_id, idx = i) %>%
    as.data.table(.)
  out <- na.omit(rbind(old, new))
  saveRDS(out, res_file)
}
foreach(ii = seq_len(100)) %:%
  foreach(ss = sims$s_id) %dopar%
  big_loop(sims, ss, ii)

