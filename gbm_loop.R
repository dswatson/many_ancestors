### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(ranger)
library(gbm)
library(ppcor)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

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
    beta <- beta[-u_idx]
    gamma <- gamma[-u_idx]
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


#' @param x Design matrix.
#' @param y Outcome vector.
#' @param f Regression method, either \code{"lasso"} or \code{"gbm"}.
#' 

# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, f) {
  n <- nrow(x)
  if (f == 'lasso') {
    trn <- sample(n, round(0.8 * n))
    tst <- seq_len(n)[-trn]
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
    eps <- y_hat - y[tst]
    mse <- colMeans(eps^2)
    betas <- coef(fit, s = fit$lambda)[-1, which.min(mse)]
  } else if (f == 'gbm') {
    d <- data.frame(y = y, x)
    fit <- gbm(y ~ ., data = d, distribution = 'gaussian', n.trees = 1500, 
               shrinkage = 0.1, interaction.depth = 1, train.fraction = 0.8, 
               n.cores = 1)
    n_tree <- gbm.perf(fit, plot.it = FALSE, method = 'test')
    betas <- summary(fit, plotit = FALSE, n.trees = n_tree, order = FALSE)$rel.inf
  } 
  out <- ifelse(betas == 0, 0, 1)
  return(out)
}


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param B Number of complementary pairs to draw for stability selection.
#' 

# Compute (de)activation rates for X -> Y and Y -> X
rate_fn <- function(z, x, y, linear, B) {
  # Preliminaries
  n <- nrow(z)
  d_z <- ncol(z)
  zx <- cbind(z, x)
  zy <- cbind(z, y)
  f <- ifelse(linear, 'lasso', 'gbm')
  # Compute disconnections and (de)activations per subsample
  fit_fn <- function(b) {
    # Take complementary subsets
    i_set <- sample(n, round(0.5 * n))
    j_set <- seq_len(n)[-i_set]
    # Compute active sets
    s <- data.frame(
      y0 = c(l0(z[i_set, ], y[i_set], f), NA_real_, 
             l0(z[j_set, ], y[j_set], f), NA_real_), 
      y1 = c(l0(zx[i_set, ], y[i_set], f), 
             l0(zx[j_set, ], y[j_set], f)),
      x0 = c(l0(z[i_set, ], x[i_set], f), NA_real_, 
             l0(z[j_set, ], x[j_set], f), NA_real_), 
      x1 = c(l0(zy[i_set, ], x[i_set], f), 
             l0(zy[j_set, ], x[j_set], f))
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


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param gamma Omission threshold.
#' 

# Wrap it up
cbr_fn <- function(z, x, y, linear, gamma) {
  # Compute rates for each z
  res <- rate_fn(z, x, y, linear, B = 50)
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
  out <- data.table(method = 'cbr', g_hat = decision) 
  return(out)
}

################################################################################


# Initialize
out <- data.table(
  method = NA, g_hat = NA, s_id = NA, idx = NA
)
res_file <- './results/gbm_benchmark.rds'
saveRDS(out, res_file)

sims <- expand.grid(n = c(2500, 5000, 1e4), g = c('xy', 'ci', 'na')) %>%
  mutate(sp = 0.5, d_z = 100, rho = 0.25, r2 = 2/3, lin_pr = 1/5, 
         s_id = row_number()) %>%
  as.data.table(.)

gbm_loop <- function(sims_df, sim_id, i) {
  # Housekeeping
  sdf <- sims_df[s_id == sim_id]
  # Simulate data
  sim_obj <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, sp = sdf$sp, 
                     r2 = sdf$r2, lin_pr = sdf$lin_pr, g = sdf$g)
  # Extract data
  dat <- sim_obj$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  y <- dat$y
  # Confounder blanket regression
  df_b <- cbr_fn(z, x, y, linear = FALSE, gamma = 0.25) 
  # Import, export
  old <- readRDS(res_file)
  new <- df_b %>%
    mutate(s_id = sim_id, idx = i) %>%
    as.data.table(.)
  out <- na.omit(rbind(old, new))
  saveRDS(out, res_file)
}
foreach(ii = seq_len(100)) %:%
  foreach(ss = sims$s_id) %dopar%
  gbm_loop(sims, ss, ii)













