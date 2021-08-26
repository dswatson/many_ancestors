### Simulations for subgraph discovery with many ancestors ###

# Load libraries, register cores
library(data.table)
library(glmnet)
library(bestsubset)
library(randomForest)
library(tidyverse)
library(ggsci)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param k Number of nonzero weights from Z to children. Either a scalar, in  
#'   which case it applies to both X and Y, or a length-two vector.
#' @param snr Signal-to-noise ratio for true data generating functions. Either a 
#'   scalar, in which case it applies to both X and Y, or a length-two vector.
#' @param xzr X-to-Z ratio, i.e. the ratio of X-signal to Z-signal for 
#'   generating Y when X -> Y.
#' @param form Should structural equations be \code{"linear"} or 
#'   \code{"nonlinear"}?

# Data simulation function
sim_dat <- function(n, d_z, rho, k, snr, xzr, form) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  if (rho != 0) {
    Sigma <- toeplitz(rho^(0:(d_z - 1)))
    z <- z %*% chol(Sigma)
  } 
  dimnames(z) <- list(NULL, paste0('z', seq_len(d_z)))
  # Nonlinear transformations?
  if (form != 'linear') {
    idx <- split(sample.int(d_z), sort(seq_len(d_z) %% 5))
    names(idx) <- c('lexp', 'sqrt', 'sq', 'relu', 'orig')
    zz <- z
    zz[, idx$lexp] <- log(1 + exp(z[, idx$lexp]))
    zz[, idx$sqrt] <- sqrt(abs(z[, idx$sqrt]))
    zz[, idx$sq] <- z[, idx$sq]^2
    zz[, idx$relu] <- ifelse(z[, idx$relu] > 0, z[, idx$relu], 0)
  }
  # Draw weights for Z
  sim_wts <- function(n, d_z, k) {
    nonzero <- sample(d_z, k)
    amplitude <- seq(1, 20, length.out = k) # Why 20?
    signs <- sample(c(1, -1), size = d_z, replace = TRUE)
    beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
    return(beta)
  }
  # Generate noise
  sim_noise <- function(signal, snr) {
    noise_factor <- sqrt(var(signal) / snr)
    noise <- rnorm(n, sd = noise_factor) 
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  beta_z_to_x <- sim_wts(n, d_z, k)
  beta_z_to_y <- sim_wts(n, d_z, k)
  if (form == 'linear') {
    signal_x <- as.numeric(z %*% beta_z_to_x)
    signal_y <- as.numeric(z %*% beta_z_to_y)
  } else {
    signal_x <- as.numeric(zz %*% beta_z_to_x)
    signal_y <- as.numeric(zz %*% beta_z_to_y)
  }
  x <- signal_x + sim_noise(signal_x, snr)
  y0 <- signal_y + sim_noise(signal_y, snr)
  ### Alternative scenario: X -> Y
  signal_z_to_y <- signal_y
  sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
  beta_x_to_y <- sigma_x_to_y / sd(x)
  signal_y <- signal_z_to_y + x * beta_x_to_y
  y1 <- signal_y + sim_noise(signal_y, snr)
  # Export
  out <- list('dat' = data.table(z, 'x' = x, 'y0' = y0, 'y1' = y1),
              'wts' = data.table('x' = beta_z_to_x, 'y' = beta_z_to_y))
  return(out)
}

#' @param m Number of nested models to fit.
#' @param max_x Number of predictors in largest model.
#' @param min_d Number of predictors in smallest model.
#' @param decay Exponential decay parameter

subsets <- function(m, max_d, min_d, decay) {
  unique(round(min_d + ((max_d - min_d) / m^decay) * seq_len(m)^decay))
}


#' @param trn_x Training set of predictors.
#' @param trn_y Training outcomes.
#' @param tst_x Test set of predictors.
#' @param tst_y Test outcomes.
#' @param f Regression function to use. 

f_fn <- function(trn_x, trn_y, tst_x, tst_y, f) {
  if (f == 'lasso') {
    fit <- glmnet(trn_x, trn_y, intercept = FALSE)
    y_hat <- predict(fit, newx = tst_x, s = fit$lambda)
    betas <- coef(fit, s = fit$lambda)[2:(d_z + 1), ]
  } else if (f == 'step') {
    fit <- fs(trn_x, trn_y, intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = tst_x)
    betas <- coef(fit)[1:d_z, ]
  } else if (f == 'rf') {
    fit <- randomForest(trn_x, trn_y, ntree = 500)
    vimp <- data.frame('feature' = colnames(trn_x), 
                           'imp' = as.numeric(importance(fit))) %>%
      arrange(desc(imp))
    beta <- double(length = d_z)
    names(beta) <- paste0('z', seq_len(d_z))
    s <- subsets(m = 10, max_d = d_z, min_d = 5, decay = 2)
    rfe_loop <- function(k) {
      tmp_x <- trn_x[, vimp$feature[seq_len(s[k])]]
      tmp_f <- randomForest(tmp_x, trn_y, ntree = 200)
      beta[colnames(tmp_x)] <- as.numeric(importance(tmp_f))
      out <- list('y_hat' = predict(tmp_f, newdata = tst_x), 'beta' = beta)
      return(out)
    }
    rf_out <- foreach(kk = seq_along(s)) %do% rfe_loop(kk)
    y_hat <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$y_hat)
    betas <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$beta)
  }
  epsilon <- y_hat - tst_y
  mse <- colMeans(epsilon^2)
  out <- list('beta' = betas[, which.min(mse)], 
               'eps' = epsilon[, which.min(mse)])
  return(out)
}


#' @param b Simulation index.
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param k Number of nonzero weights from Z to children.
#' @param snr Signal-to-noise ratio. 
#' @param xzr X-to-Z ratio.
#' @param form Functional form for structural equations.
#' @param l Norm to use, either 0, 1, or 2.

our_fn <- function(b, n, d_z, rho, k, snr, xzr, form, l) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho, k, snr, xzr, form)
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Compute coefficients
  fit_fn <- function(h, f) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    zy <- cbind(z, y)
    betas <- list(
      f_fn(z[trn, ], y[trn], z[tst, ], y[tst], f)$beta,
      f_fn(zx[trn, ], y[trn], zx[tst, ], y[tst], f)$beta,
      f_fn(z[trn, ], x[trn], z[tst, ], x[tst], f)$beta,
      f_fn(zy[trn, ], x[trn], zy[tst, ], x[tst], f)$beta
    )
    v <- sapply(seq_len(4), function(f_idx) {
      switch(l, 'l0' = sum(betas[[f_idx]] != 0), 
                'l1' = sum(abs(betas[[f_idx]])),
                'l2' = sum(betas[[f_idx]]^2))
    })
    w0 <- v[1] - v[2]
    w1 <- v[3] - v[4]
    delta <- w0 - w1
    # Export
    out <- data.table('h' = h, 'f' = f, 'delta' = delta)
    return(out)
  }
  out <- foreach(a = c('h0', 'h1'), .combine = rbind) %:%
    foreach(b = c('lasso', 'step'), .combine = rbind) %do% fit_fn(h = a, f = b)
  return(out)
}
res <- foreach(i = seq_len(500), .combine = rbind) %dopar% 
  our_fn(i, n = 500, d_z = 50, rho = 0.3, k = 25, snr = 3, xzr = 1, l = 'l0')


sample_sizes <- c(200, 500, 1000)
dim_z <- c(20, 50, 100)
auto_cor <- c(0, 0.3, 0.5)
sparsity <- c(0.1, 0.5, 0.9)
signals <- c(1, 3, 5)
x_strength <- c(0.5, 1, 2)

# Entner
entner_fn <- function(b, n, d_z, rho, k, snr, xzr, form, alpha) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho, k, snr, xzr, form)
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Fit models
  fit_fn <- function(h, f) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    # Three different designs
    f0 <- f_fn(z[trn, ], y[trn], z[tst, ], y[tst], f)
    f1 <- f_fn(zx[trn, ], y[trn], zx[tst, ], y[tst], f)
    f2 <- f_fn(z[trn, ], x[trn], z[tst, ], x[tst], f)
    # Evaluate rules
    delta <- abs(f0$eps) - abs(f1$eps)
    loco_p <- wilcox.test(delta, alt = 'greater')$p.value
    r1 <- f0$beta != 0 & f1$beta == 0
    r2 <- loco_p <= alpha | (f2$beta != 0 & f0$beta == 0)
    # Export results
    out <- data.table('b' = b, 'h' = h, 'f' = f, 'r1' = r1, 'r2' = r2)
    return(out)
  }
  out <- foreach(aa = c('h0', 'h1'), .combine = rbind) %:%
    foreach(bb = c('lasso', 'step'), .combine = rbind) %do% 
    fit_fn(aa, bb)
  return(out)
}
res <- foreach(i = seq_len(200), .combine = rbind) %dopar% 
  entner_fn(i, n = 500, d_z = 50, rho = 0.3, k = 25, snr = 5, xzr = 1, 
            form = 'linear', alpha = 0.05)

n <- 500; d_z <- 50; rho <- 0.3; k <- 25; snr <- 5; xzr <- 1; 
form <- 'linear'; alpha <- 0.05


# Need to record which Z's are "instrumental", i.e. satisfy
# Z -> X -> Y and Z \ind Y | X. Define k := number of such Z's.
# Then we should expect k many Z's to satisfy R1 under h1 and 
# 0 to satisfy R1 under h0. We should also expect d_z many Z's
# to satisfy R2 under h0 and 0 to satisfy under h1.


# Score-based benchmark?













