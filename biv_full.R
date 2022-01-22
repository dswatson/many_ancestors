### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(ranger)
library(ppcor)
library(kpcalg)
library(tidyverse)
library(doMC)
registerDoMC(16)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param sp Sparsity of the connections from background to foreground.
#' @param r2 Proportion of variance explained by endogenous features.
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' 
#' @import data.table

# Data simulation function
sim_dat <- function(n, d_z, rho, sp, r2, lin_pr) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  var_z <- 1 / d_z # Does this make a difference?
  Sigma <- toeplitz(rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Set weights
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
  ### Null scenario: X \indep Y | Z
  x <- signal_x + sim_noise(signal_x, r2)
  y0 <- signal_y + sim_noise(signal_y, r2)
  ### Alternative scenario: X -> Y
  signal_z_to_y <- signal_y
  xzr <- 1 / (k + 1)
  sigma_xy <- sqrt(xzr * var(signal_z_to_y))
  gamma_x <- sigma_xy / sd(x)
  signal_y <- signal_z_to_y + x * gamma_x
  y1 <- signal_y + sim_noise(signal_y, r2)
  # Export
  params <- list(
    'n' = n, 'd_z' = d_z, 'rho' = rho, 'sp' = sp, 'r2' = r2, 'lin_pr' = lin_pr
  )
  out <- list(
    'dat' = data.table(z, 'x' = x, 'y0' = y0, 'y1' = y1),
    'wts' = list('beta' = beta, 'gamma' = gamma), 'params' = params
  )
  return(out)
}


#' @param m Number of nested models to fit.
#' @param max_d Number of predictors in largest model.
#' @param min_d Number of predictors in smallest model.
#' @param decay Exponential decay parameter

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
#' @import glmnet
#' @import ranger
#' @import tidyverse

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
#' @import data.table
#' @import foreach

# Compute (de)activation rates for X -> Y and Y -> X
rate_fn <- function(sim_obj, hyp, B) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  zx <- cbind(z, x)
  if (hyp == 'h0') y <- dat$y0 else y <- dat$y1
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
    out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_z + 1), h = hyp, 
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
  out <- unique(out[, .(h, z, disr, drxy, arxy, dryx, aryx)])
  return(out)
}

#' @param res Results object output by \code{rate_fn}.
#' @param hyp If X -> Y, \code{"h0"}; else if X \indep Y | Z, \code{"h1"}.
#' @param B Number of complementary pairs to draw for stability selection.
#' 
#' @import data.table
#' @import foreach

# Compute consistency lower bound
lb_fn <- function(res, hyp, B) {
  # Subset the data
  df <- na.omit(res[h == hyp, .(z, drxy, arxy, dryx, aryx)])
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
  # Export
  out <- data.table('h' = hyp, min_two)
  return(out)
}

#' @param res Results object output by \code{rate_fn}.
#' @param cons_tbl Consistency table output by \code{lb_fn}.
#' @param hyp If X -> Y, \code{"h0"}; else if X \indep Y | Z, \code{"h1"}.
#' @param order Assume X \preceq Y or Y \preceq X?
#' @param rule Detect via deactivation (\code{"R1"}) or activation (\code{"R2"})?
#' @param B Number of complementary pairs to draw for stability selection.
#' 
#' @import tidyverse
#' @import data.table

# Infer causal direction using stability selection
ss_fn <- function(res, cons_tbl, hyp, order, rule, B) {
  # Subset the data
  if (order == 'xy' & rule == 'R1') {
    r <- res[h == hyp, drxy]
  } else if (order == 'xy' & rule == 'R2') {
    r <- res[h == hyp, arxy]
  } else if (order == 'yx' & rule == 'R1') {
    r <- res[h == hyp, dryx]
  } else if (order == 'yx' & rule == 'R2') {
    r <- res[h == hyp, aryx]
  }
  r <- na.omit(r)
  # Find consistency lower bound
  lb <- cons_tbl[h == hyp, min_two]
  # Stability selection parameters
  theta <- mean(r)
  ub <- minD(theta, B) * sum(r <= theta)
  tau <- seq_len(2 * B) / (2 * B)
  # Do any features exceed the upper bound?
  dat <- data.frame(tau, err_bound = ub) %>%
    filter(tau > lb) %>%
    rowwise() %>%
    mutate(detected = sum(r >= tau)) %>% 
    ungroup(.) %>%
    mutate(surplus = ifelse(detected > err_bound, 1, 0))
  # Export
  out <- data.table(
    'h' = hyp, 'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}

# Wrap it up
cbr_fn <- function(sim, hyp, eps) {
  # Compute rates for each z, h
  res <- rate_fn(sim, hyp, B = 50)
  # Disconnected?
  if (res$disr[1] > eps) {
    decision <- 'ci'
  } else {
    # (De)activation rates
    cons_tbl <- lb_fn(res, res$h[1], B = 50)
    sum_tbl <- foreach(oo = c('xy', 'yx'), .combine = rbind) %:%
      foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
      ss_fn(res, cons_tbl, res$h[1], oo, rr, B = 50)
    if (sum_tbl[order == 'xy', sum(decision)] > 0) {
      decision <- 'xy'
    } else if (sum_tbl[order == 'yx', sum(decision)] > 0) {
      decision <- 'yx'
    } else {
      decision <- NA_character_
    }
  }
  # Export
  out <- data.table(
    method = 'cbr', h = hyp, g = decision
  ) 
  return(out)
}

################################################################################

### OTHER METHODS ###

#' @param h If X -> Y, \code{"h0"}; else if X \indep Y | Z, \code{"h1"}.
#' @param alpha Significance threshold for inferring dependence.
#' @param tau Other threshold for "inferring" independence.

# Constraint-based: for this comparison, we presume X \preceq Y
# and assume access to the true data sparsity
constr_fn <- function(sim_obj, h, alpha, tau) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  if (h == 'h0') y <- dat$y0 else y <- dat$y1
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  k <- sim_obj$params$sp * d_z
  # Entner function
  entner <- function(b) {
    z_idx <- sample(d_z, k)
    z_b <- z[, z_idx]
    in_loop <- function(j) {
      w <- z[, j]
      if (linear == TRUE) {
        ### RULE 1 ###
        pmat0 <- suppressWarnings(pcor(cbind(y, w, z_b))$p.value)
        p1.i <- pmat0[1, 2]
        pmat1 <- suppressWarnings(pcor(cbind(y, w, z_b, x))$p.value)
        p1.ii <- pmat1[1, 2]
        ### RULE 2 ###
        pmat2 <- suppressWarnings(pcor(cbind(y, x, z_b))$p.value)
        p2.i <- pmat2[1, 2]
        pmat3 <- suppressWarnings(pcor(cbind(x, w, z_b))$p.value)
        p2.ii <- pmat3[1, 2]
      } else {
        ### RULE 1 ###
        d_zb <- ncol(z_b)
        dat <- cbind(y, w, z_b)
        p1.i <- kernelCItest(
          x = 1, y = 2, S = 3:(d_zb + 2),
          suffStat = list(data = dat, ic.method = 'hsic.gamma')
        )
        dat <- cbind(y, w, z_b, x)
        p1.ii <- kernelCItest(
          x = 1, y = 2, S = 3:(d_zb + 3),
          suffStat = list(data = dat, ic.method = 'hsic.gamma')
        )
        ### RULE 2 ###
        dat <- cbind(y, x, z_b)
        p2.i <- kernelCItest(
          x = 1, y = 2, S = 3:(d_zb + 2),
          suffStat = list(data = dat, ic.method = 'hsic.gamma')
        )
        dat <- cbind(x, w, z_b)
        p2.ii <- kernelCItest(
          x = 1, y = 2, S = 3:(d_zb + 2),
          suffStat = list(data = dat, ic.method = 'hsic.gamma')
        )
      }
      # Apply rules
      r1 <- ifelse(p1.i <= alpha & p1.ii >= tau, 1, 0)
      r2 <- ifelse(p2.i >= tau | (p2.ii <= alpha & p1.i >= tau), 1, 0)
      # Export
      out <- data.table(r1, r2)
      return(out)
    }
    omega <- seq_len(d_z)[-z_idx]
    out <- foreach(jj = sample(omega, 20), .combine = rbind) %do% 
      in_loop(jj)
    return(out)
  }
  # Apply Entner's rules with 50 random subsets and 20 candidates per subset
  df <- foreach(bb = seq_len(50), .combine = rbind) %do%
    entner(bb)
  # Selecting different decision thresholds based on experimentation
  if (df[, sum(r1) / .N] > 1/200) {
    g <- 'xy'
  } else if (df[, sum(r2) / .N] > 1/4) {
    g <- 'ci'
  } else {
    g <- NA_character_
  }
  out <- data.table(method = 'constr', h, g)
  return(out)
}

#' @param x Design matrix.
#' @param y Response vector.
#' @param f Function class.

# MSE-scoring subroutine
mse_fn <- function(x, y, f) {
  # Split data
  n <- nrow(x)
  trn <- sample(n, size = round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Fit models
  if (f == 'lasso') {
    fit <- cv.glmnet(x[trn, ], y[trn])
    y_hat <- predict(fit, newx = x[tst, ], s = 'lambda.min')
  } else if (f == 'rf') {
    fit <- ranger(x = x[trn, ], y = y[trn], num.threads = 1)
    y_hat <- predict(fit, data = x[tst, ], num.threads = 1)$predictions
  }
  # Evaluate errors, export
  eps <- y_hat - y[tst]
  mse <- mean(eps^2)
  return(mse)
}

#' @param h If X -> Y, \code{"h0"}; else if X \indep Y | Z, \code{"h1"}.

# Score function evaluates three different DGPs
score_fn <- function(sim_obj, h) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  # Pick right y
  if (h == 'h0') y <- dat$y0 else y <- dat$y1
  # Pick right f
  f <- ifelse(linear, 'lasso', 'rf')
  # Score X -> Y
  mse_xy <- mse_fn(cbind(z, x), y, f)
  # Score Y -> X
  mse_yx <- mse_fn(cbind(z, y), x, f)
  # Score X ~ Y
  mse_x <- mse_fn(z, x, f)
  mse_y <- mse_fn(z, y, f)
  # Summarize
  df <- data.table(
    h, mse = c(mse_x + mse_xy, mse_y + mse_yx, mse_x + mse_y),
    g = c('xy', 'yx', 'ci')
  )
  # Export
  out <- data.table(
    method = 'score', h, g = df[which.min(mse), g]
  )
  return(out)
}

################################################################################

### PUT IT ALL TOGETHER ###

# Big ol' wrapper
big_loop <- function(sims_df, sim_id, i) {
  # Simulate data, extract ground truth
  sdf <- sims_df[s_id == sim_id]
  sim_obj <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, 
                     sp = sdf$sp, r2 = sdf$r2, lin_pr = sdf$lin_pr) 
  # Confounder blanket regression
  df_b <- foreach(hh = c('h0', 'h1'), .combine = rbind) %do%
    cbr_fn(sim_obj, hh, eps = 0.25)
  # Constraint function
  df_c <- foreach(hh = c('h0', 'h1'), .combine = rbind) %do% 
    constr_fn(sim_obj, hh, alpha = 0.1, tau = 0.5)
  # Score function
  df_s <- foreach(sim_obj, hh = c('h0', 'h1'), .combine = rbind) %do%
    score_fn(sim_obj, hh)
  # Export
  out <- rbind(df_b, df_c, df_s) 
  out[, s_id := sim_id]
  out[, idx := i]
  # Check in
  if (i == 100) {
    cat(paste('s_id =', sim_id, 'complete.\n', 
        max(sims_df$s_id) - sim_id, 'more to go...\n'))
  }
  return(out)
}

# Execute in parallel
sims <- expand.grid(n = c(1000, 2000, 4000), sp = c(0.25, 0.5, 0.75)) %>%
  mutate(d_z = 100, rho = 0.25, r2 = 2/3, lin_pr = 1/5, s_id = row_number()) %>%
  as.data.table(.)
res <- foreach(ss = sims$s_id, .combine = rbind) %:%
  foreach(ii = seq_len(100), .combine = rbind) %dopar%
  big_loop(sims, ss, ii)
fwrite(res, './results/biv_benchmark.csv')





