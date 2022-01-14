### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('./Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(randomForest)
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
#' @param max_x Number of predictors in largest model.
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
#' @param d_z Dimensionality of ancestor set Z.
#' @param f Regression method, either \code{"lasso"} or \code{"rf"}.
 
# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, trn, tst, d_z, f) {
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
    betas <- coef(fit, s = fit$lambda)[2:(d_z + 1), ]
  } else if (f == 'step') {
    fit <- fs(x[trn, ], y[trn], intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = x[tst, ])
    betas <- coef(fit)[1:d_z, ]
  } else if (f == 'rf') {
    fit <- randomForest(x[trn, ], y[trn], ntree = 200)
    vimp <- data.frame('feature' = colnames(x), 
                       'imp' = as.numeric(importance(fit))) %>%
      arrange(desc(imp))
    beta <- double(length = d_z)
    names(beta) <- paste0('z', seq_len(d_z))
    s <- subsets(m = 10, max_d = d_z, min_d = 5, decay = 2)
    rfe_loop <- function(k) {
      tmp_x <- x[trn, vimp$feature[seq_len(s[k])]]
      tmp_f <- randomForest(tmp_x, y[trn], ntree = 50)
      tmp_v <- data.frame('feature' = colnames(tmp_x), 
                          'imp' = as.numeric(importance(tmp_f))) %>%
        filter(grepl('z', feature))
      beta[tmp_v$feature] <- tmp_v$imp
      out <- list('y_hat' = predict(tmp_f, newdata = x[tst, ]), 'beta' = beta)
      return(out)
    }
    rf_out <- foreach(kk = seq_along(s)) %do% rfe_loop(kk)
    y_hat <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$y_hat)
    y_hat <- cbind(y_hat, predict(fit, newdata = x[tst, ]))
    betas <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$beta)
    betas <- cbind(betas, as.numeric(importance(fit))[seq_len(d_z)])
  }
  epsilon <- y_hat - y[tst]
  mse <- colMeans(epsilon^2)
  betas <- betas[, which.min(mse)]
  out <- ifelse(betas == 0, 0, 1)
  return(out)
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute (de)activation rates for X -> Y and Y -> X
rate_fn <- function(sim_obj, B) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  zx <- cbind(z, x)
  # Linear or nonlinear?
  f <- ifelse(sim_obj$params$lin_pr == 1, 'lasso', 'rf')
  # Compute feature weights
  fit_fn <- function(h, b) {
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
      y0 = c(l0(z, y, i_trn, i_tst, d_z, f), l0(z, y, j_trn, j_tst, d_z, f)), 
      y1 = c(l0(zx, y, i_trn, i_tst, d_z, f), l0(zx, y, j_trn, j_tst, d_z, f)),
      x0 = c(l0(z, x, i_trn, i_tst, d_z, f), l0(z, x, j_trn, j_tst, d_z, f)), 
      x1 = c(l0(zy, x, i_trn, i_tst, d_z, f), l0(zy, x, j_trn, j_tst, d_z, f))
    )
    # Record (de)activations
    d_xy <- s$y0 == 1 & s$y1 == 0
    a_xy <- s$x0 == 0 & s$x1 == 1
    d_yx <- s$x0 == 1 & s$x1 == 0
    a_yx <- s$y0 == 0 & s$y1 == 1
    # Export
    out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_z), h, 
                      d_xy, a_xy, d_yx, a_yx,
                      z = rep(seq_len(d_z), times = 2))
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
  out <- unique(out[, .(h, z, drxy, arxy, dryx, aryx)])
  return(out)
}


# Compute consistency lower bound
lb_fn <- function(res, hyp, B) {
  # Subset the data
  df <- res[h == hyp, .(z, drxy, arxy, dryx, aryx)]
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
  # Find consistency lower bound
  lb <- cons_tbl[h == hyp, min_two]
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
    'h' = hyp, 'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}


# Big ol' wrapper
big_loop <- function(sims_df, sim_id, i, B) {
  # Simulate data, extract ground truth
  sdf <- sims_df[s_id == sim_id]
  sim <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, 
                 sp = sdf$sp, r2 = sdf$r2, lin_pr = sdf$lin_pr)
  d_xy_true <- sim$wts$beta == 1 & sim$wts$gamma == 0
  a_xy_true <- sim$wts$beta == 0 & sim$wts$gamma == 1
  # Compute (de)activation rates for each z, h
  res <- rate_fn(sim, B)
  # Consistent lower bound
  cons_tbl <- foreach(hh = c('h0', 'h1'), .combine = rbind) %do%
    lb_fn(res, hh, B)
  # Stable upper bound
  sum_tbl <- foreach(hh = c('h0', 'h1'), .combine = rbind) %:%
    foreach(oo = c('xy', 'yx'), .combine = rbind) %:%
    foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
    ss_fn(res, cons_tbl, hh, oo, rr, B)
  # Export
  sum_tbl[, s_id := sim_id]
  sum_tbl[, idx := i]
  return(sum_tbl)
}


### SIMULATION GRID ###
sims <- expand.grid(
  n = c(500, 1000, 2000), d_z = c(50, 100, 200), rho = c(0, 0.5),
  sp = c(0.25, 0.5, 0.75), r2 = c(1/3, 1/2, 2/3), lin_pr = 1
)
# Linear?
if (linear == TRUE) {
  sims$lin_pr <- 1
  lab <- 'linear_sim.csv'
} else {
  sims$lin_pr <- 1/5
  lab <- 'nonlinear_sim.csv'
}
# Index, data table-ify
sims$s_id <- seq_len(nrow(sims))
sims <- as.data.table(sims)

# Compute in parallel
res <- foreach(ss = sims$s_id, .combine = rbind) %:%
  foreach(ii = seq_len(100), .combine = rbind) %dopar%
  big_loop(sims, ss, ii, B = 50)
res[, hit_rate := sum(decision) / .N, by = .(h, order, rule, s_id)]
res <- unique(res[, .(s_id, h, order, rule, hit_rate)])
res <- merge(res, sims, by = 's_id')
fwrite(res, lab)

# Polish for plotting
res[, hit := ifelse(h == 'h1' & order == 'xy', TRUE, FALSE)]
res[, setting := paste(h, order, rule, sep = ',')]
res[, setting := factor(setting, levels = c('h0,yx,R1', 'h0,yx,R2', 'h0,xy,R1', 'h0,xy,R2',
                                            'h1,yx,R1', 'h1,yx,R2', 'h1,xy,R1', 'h1,xy,R2'))]
res[, n := paste0('n=', n)]
res[, n := factor(n, levels = c('n=500', 'n=1000', 'n=2000'))]
res[, d_z := paste0('d=', d_z)]
res[, d_z := factor(d_z, levels = c('d=50', 'd=100', 'd=200'))]
plot_loop <- function(r, s, p) {
  # Labels
  if (r == 1/3) {
    rlab <- 'lo' 
  } else if (r == 1/2) {
    rlab <- 'me'
  } else if (r == 2/3) {
    rlab <- 'hi'
  }
  if (s == 0.25) {
    slab <- 'lo'
  } else if (s == 0.5) {
    slab <- 'me'
  } else if (s == '0.75') {
    slab <- 'hi'
  }
  if (p == 0) {
    plab <- 'lo'
  } else if (p == 0.5) {
    plab <- 'hi'
  }
  # Plot
  p <- ggplot(res[r2 == r & sp == s & rho == p], 
              aes(setting, hit_rate, color = hit)) + 
    geom_bar(stat = 'identity') + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45)) + 
    facet_grid(d_z ~ n)
  # Export
  ggsave(paste0('./plots/linear_bivariate_r2=', rlab, ',sp=', slab, 
                ',rho=', plab, '.png'))
}
foreach(aa = c(1/3, 1/2, 2/3)) %:%
  foreach(bb = c(0.25, 0.5, 0.75)) %:%
  foreach(cc = c(0, 0.5)) %do%
  plot_loop(aa, bb, cc)






