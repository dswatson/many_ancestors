### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('./Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param sp_x Proportion of zero weights from Z to X. 
#' @param sp_y Proportion of zero weights from Z to Y.
#' @param r2_x Proportion of variance explained in the DGP for X.
#' @param r2_y Proportion of variance explained in the DGP for Y.
#' @param wt_type Type of weight, either \code{"equal"} or \code{"unequal"}.
#' @param xzr Ratio of X-signal to Z-signal when X -> Y.
#' @param form Should edges from Z encode \code{"linear"} or \code{"nonlinear"}
#'   structural equations?

# Data simulation function
sim_dat <- function(n, d_z, rho, sp_x, sp_y, r2_x, r2_y, wt_type, xzr, form) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  var_z <- 1 / d_z
  Sigma <- toeplitz(rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Nonlinear transformations?
  if (form == 'nonlinear') {
    idx <- split(sample.int(d_z), sort(seq_len(d_z) %% 5))
    names(idx) <- c('sq', 'sqrt', 'sftpls', 'relu', 'orig')
    zz <- z
    zz[, idx$sq] <- z[, idx$sq]^2
    zz[, idx$sqrt] <- sqrt(abs(z[, idx$sqrt]))
    zz[, idx$sftpls] <- log(1 + exp(z[, idx$lexp]))
    zz[, idx$relu] <- ifelse(z[, idx$relu] > 0, z[, idx$relu], 0)
  }
  # Set weights
  beta <- gamma <- double(length = d_z)
  k <- round((1 - sp_x) * d_z)
  if (wt_type == 'equal') {
    beta[sample(d_z, k)] <- 1
    gamma[sample(d_z, k)] <- 1
  } else if (wt_type == 'unequal') {
    beta[sample(d_z, k)] <- seq(from = 1, to = 10, length.out = k)
    gamma[sample(d_z, k)] <- seq(from = 1, to = 10, length.out = k)
  }
  # Random signs
  beta <- beta * sample(c(1, -1), d_z, replace = TRUE)
  gamma <- gamma * sample(c(1, -1), d_z, replace = TRUE)
  if (form == 'linear') {
    signal_x <- as.numeric(z %*% beta)
    signal_y <- as.numeric(z %*% gamma)
  } else {
    signal_x <- as.numeric(zz %*% beta)
    signal_y <- as.numeric(zz %*% gamma)
  }
  # Generate noise
  sim_noise <- function(signal, r2) {
    var_mu <- var(signal)
    var_noise <- (var_mu - r2 * var_mu) / r2
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  x <- signal_x + sim_noise(signal_x, r2_x)
  y0 <- signal_y + sim_noise(signal_y, r2_y)
  ### Alternative scenario: X -> Y
  signal_z_to_y <- signal_y
  sigma_xy <- sqrt(xzr * var(signal_z_to_y))
  alpha <- sigma_xy / sd(x)
  signal_y <- signal_z_to_y + x * alpha
  y1 <- signal_y + sim_noise(signal_y, r2_y)
  # Export
  out <- list(
    'dat' = data.table(z, 'x' = x, 'y0' = y0, 'y1' = y1),
    'wts' = list('beta' = beta, 'gamma' = gamma, 'alpha' = alpha)
  )
  return(out)
}

#' @param x Design matrix.
#' @param y Outcome vector.
#' @param trn Training indices.
#' @param tst Test indices.
#' @param d_z Dimensionality of ancestor set Z.
 
# Fit regressions, return bit vector for feature selection
f_fn <- function(x, y, trn, tst, d_z) {
  fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
  y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
  betas <- coef(fit, s = fit$lambda)[2:(d_z + 1), ]
  epsilon <- y_hat - y[tst]
  mse <- colMeans(epsilon^2)
  betas <- betas[, which.min(mse)]
  out <- ifelse(betas == 0, 0, 1)
  return(out)
}


#' @param sims Grid of simulation settings.
#' @param s_idx Index for simulation setting.
#' @param i Index for iteration of said setting.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute (de)activation rates for X -> Y and Y -> X
rate_fn <- function(sims, s_idx, i, B) {
  # Simulate data
  sdf <- sims[idx == s_idx]
  sim <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, 
                 sp_x = sdf$sp, sp_y = sdf$sp, r2_x = sdf$r2, r2_y = sdf$r2, 
                 wt_type = 'equal', xzr = 1, form = 'linear')
  d_xy_true <- sim$wts$beta == 1 & sim$wts$gamma == 0
  a_xy_true <- sim$wts$beta == 0 & sim$wts$gamma == 1
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Compute feature weights
  fit_fn <- function(b, h) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    zy <- cbind(z, y)
    # Take complementary subsets
    i_set <- sample(n, round(0.5 * n))
    i_trn <- sample(i_set, round(0.8 * length(i_set)))
    i_tst <- setdiff(i_set, i_trn)
    j_set <- seq_len(n)[-i_set]
    j_trn <- sample(j_set, round(0.8 * length(j_set)))
    j_tst <- setdiff(j_set, j_trn)
    # Compute active sets
    s <- data.frame(
      f_fn(z, y, i_trn, i_tst, d_z), f_fn(z, y, j_trn, j_tst, d_z), 
      f_fn(zx, y, i_trn, i_tst, d_z), f_fn(zx, y, j_trn, j_tst, d_z),
      f_fn(z, x, i_trn, i_tst, d_z), f_fn(z, x, j_trn, j_tst, d_z), 
      f_fn(zy, x, i_trn, i_tst, d_z), f_fn(zy, x, j_trn, j_tst, d_z)
    )
    colnames(s) <- c('y0i', 'y0j', 'y1i', 'y1j', 'x0i', 'x0j', 'x1i', 'x1j')
    # Record (de)activations
    d_xy_i <- s$y0i == 1 & s$y1i == 0
    d_xy_j <- s$y0j == 1 & s$y1j == 0
    a_xy_i <- s$x0i == 0 & s$x1i == 1
    a_xy_j <- s$x0j == 0 & s$x1j == 1
    a_yx_i <- s$y0i == 0 & s$y1i == 1
    a_yx_j <- s$y0j == 0 & s$y1j == 1
    d_yx_i <- s$x0i == 1 & s$x1i == 0
    d_yx_j <- s$x0j == 1 & s$x1j == 0
    # Export
    out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_z), h, 
                      d_xy = c(d_xy_i, d_xy_j), a_xy = c(a_xy_i, a_xy_j), 
                      d_yx = c(d_yx_i, d_yx_j), a_yx = c(a_yx_i, a_yx_j), 
                      d_xy_true, a_xy_true, z = rep(seq_len(d_z), times = 2))
    return(out)
  }
  out <- foreach(aa = c('h0', 'h1'), .combine = rbind) %:%
    foreach(bb = seq_len(B), .combine = rbind) %do% 
    fit_fn(h = aa, b = bb)
  # Compute rates
  out[, drxy := sum(d_xy) / .N, by = .(h, z)]
  out[, arxy := sum(a_xy) / .N, by = .(h, z)]
  out[, dryx := sum(d_yx) / .N, by = .(h, z)]
  out[, aryx := sum(a_yx) / .N, by = .(h, z)]
  # Tidy up, export
  out[, sim_idx := s_idx]
  out[, idx := i]
  out <- unique(out[, .(sim_idx, idx, h, z, 
                        drxy, arxy, dryx, aryx, 
                        d_xy_true, a_xy_true)])
  return(out)
}


# Compute consistency lower bound
lb_fn <- function(res, i, hyp) {
  # Subset the data
  df <- res[idx == i & h == hyp, .(z, drxy, arxy, dryx, aryx)]
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
  lie_df <- foreach(tt = seq(0.01, 1, 0.01), .combine = rbind) %do% 
    lies(tt)
  # Compute minimal thresholds
  min_int <- lie_df[int_err == 0, min(tau)]
  min_ext <- lie_df[ext_err == 0, min(tau)]
  min_two <- lie_df[int_err == 0 & ext_err == 0, min(tau)] # It's always ext
  # Export
  out <- data.table('idx' = i, 'h' = hyp, min_int, min_ext, min_two)
  return(out)
}


# Infer causal direction using stability selection
ss_fn <- function(res, i, hyp, order, rule, B) {
  # Subset the data
  if (order == 'xy' & rule == 'R1') {
    r <- res[idx == i & h == hyp, drxy]
  } else if (order == 'xy' & rule == 'R2') {
    r <- res[idx == i & h == hyp, arxy]
  } else if (order == 'yx' & rule == 'R1') {
    r <- res[idx == i & h == hyp, dryx]
  } else if (order == 'yx' & rule == 'R2') {
    r <- res[idx == i & h == hyp, aryx]
  }
  # Find consistency lower bound
  lb <- cons_tbl[idx == i & h == hyp, min_two]
  # Stability selection parameters
  q <- sum(r)
  theta <- q / length(r)
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
    'idx' = i, 'h' = hyp, 'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}

# Simulation grid
sims <- expand.grid(
  n = c(500, 2000), d_z = c(50, 200), rho = c(0, 0.75),
  sp = c(0.25, 0.75), r2 = c(1/3,  2/3)
)
sims <- expand.grid(
  n = c(500, 1000, 2000), d_z = c(50, 100, 200), rho = c(0, 0.25, 0.75),
  sp = c(0.25, 0.5, 0.75), r2 = c(1/3, 1/2, 2/3)
)
sims$idx <- seq_len(nrow(sims))


# Big wrapper
big_loop <- function(iters, B) {
  res <- foreach(ss = sims$idx, .combine = rbind) %:%
    foreach(ii = seq_len(iters), .combine = rbind) %dopar% 
    rate_fn(sims, ss, ii, B)
  cons_tbl <- foreach(ii = seq_len(iters), .combine = rbind) %:%
    foreach(hh = c('h0', 'h1'), .combine = rbind) %dopar%
    lb_fn(res, ii, hh)
  sum_tbl <- foreach(ii = seq_len(100), .combine = rbind) %:%
    foreach(hh = c('h0', 'h1'), .combine = rbind) %:%
    foreach(oo = c('xy', 'yx'), .combine = rbind) %:%
    foreach(rr = c('R1', 'R2'), .combine = rbind) %dopar%
    ss_fn(ii, hh, oo, rr, B = 50)
  out <- sum_tbl[, sum(decision) / .N, by = .(h, order, rule)]
}



