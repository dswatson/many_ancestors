### Simulations for subgraph discovery with many ancestors ###

# Load libraries, register cores
library(data.table)
library(glmnet)
library(bestsubset)
library(tidyverse)
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

# Data simulation function
sim_dat <- function(n, d_z, rho, k, snr, xzr) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  if (rho > 0) {
    Sigma <- toeplitz(rho^(0:(d_z - 1)))
    z <- z %*% chol(Sigma)
  } 
  dimnames(z) <- list(NULL, paste0('z', seq_len(d_z)))
  # Draw Z weights
  sim_z_wts <- function(n, d_z, k) {
    nonzero <- sample(d_z, k)
    amplitude <- seq(1, 20, length.out = k) # Why 20?
    signs <- sample(c(1, -1), size = d_z, replace = TRUE)
    beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
    return(beta)
  }
  # Generate noise
  sim_u <- function(signal, snr) {
    noise_factor <- sqrt(var(signal) / snr)
    noise <- rnorm(n, sd = noise_factor) 
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  beta_z_to_x <- sim_z_wts(n, d_z, k)
  signal_x <- as.numeric(z %*% beta_z_to_x)
  noise_x <- sim_u(signal_x, snr)
  x <- signal_x + noise_x
  beta_z_to_y <- sim_z_wts(n, d_z, k)
  signal_y <- as.numeric(z %*% beta_z_to_y)
  noise_y <- sim_u(signal_y, snr)
  y0 <- signal_y + noise_y
  ### Alternative scenario: X -> Y
  beta_z_to_y <- sim_z_wts(n, d_z, k)
  signal_z_to_y <- as.numeric(z %*% beta_z_to_y)
  sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
  beta_x_to_y <- sigma_x_to_y / sd(x)
  signal_y <- signal_z_to_y + x * c(beta_x_to_y)
  noise_y <- sim_u(signal_y, snr)
  y1 <- signal_y + noise_y
  # Export
  out <- data.table(z, 'x' = x, 'y0' = y0, 'y1' = y1)
  return(out)
}


#' @param trn_x Training set of predictors.
#' @param trn_y Training outcomes.
#' @param tst_x Test set of predictors.
#' @param tst_y Test outcomes.
#' @param f Regression function to use, either lasso or forward selection. 

beta_fn <- function(trn_x, trn_y, tst_x, tst_y, f) {
  if (f == 'lasso') {
    fit <- glmnet(trn_x, trn_y, intercept = FALSE)
    y_hat <- predict(fit, newx = tst_x, s = fit$lambda)
    mse <- colMeans((y_hat - tst_y)^2)
    beta <- as.numeric(coef(fit, s = fit$lambda[which.min(mse)]))
  } else if (f == 'step') {
    fit <- fs(trn_x, trn_y, intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = tst_x)
    mse <- colMeans((y_hat - tst_y)^2)
    beta <- coef(fit)[, which.min(mse)]
  }
  return(beta)
}


#' @param b Simulation index.
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param k Number of nonzero weights from Z to children.
#' @param snr Signal-to-noise ratio. 
#' @param xzr X-to-Z ratio.
#' @param l Norm to use, either 0, 1, or 2.

# Big ol' wrapper
test_fn <- function(b, n, d_z, rho, k, snr, xzr, l) {
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Simulate data
  dat <- sim_dat(n, d_z, rho, k, snr, xzr)
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Compute coefficients
  fit_fn <- function(h, f) {
    if (h == 0) y <- dat$y0 else y <- dat$y1
    zy <- cbind(z, y)
    betas <- list(
      beta_fn(z[trn, ], y[trn], z[tst, ], y[tst], f),
      beta_fn(zx[trn, ], y[trn], zx[tst, ], y[tst], f),
      beta_fn(z[trn, ], x[trn], z[tst, ], x[tst], f),
      beta_fn(zy[trn, ], x[trn], zy[tst, ], x[tst], f)
    )
    v <- sapply(seq_len(4), function(i) {
      if (l == 0) {
        sum(betas[[i]] != 0)
      } else if (l == 1) {
        sum(abs(betas[[i]]))
      } else if (l == 2) {
        sum(betas[[i]]^2)
    }})
    # Entner's rules
    r1_x <- v[1] - v[2]
    r2_x <- v[3] - v[1]
    x_score <- r1_x - r2_x
    r1_y <- v[3] - v[4]
    r2_y <- v[1] - v[3]
    y_score <- r1_y - r2_y
    delta1 <- x_score - y_score
    # Our method
    w0 <- v[1] - v[2]
    w1 <- v[3] - v[4]
    delta2 <- w0 - w1
    out <- data.table(
      'h' = h, 'f' = f, 'entner' = delta1, 'ours' = delta2
    )
    # Export
    return(out)
  }
  out <- foreach(a = c(0, 1), .combine = rbind) %:%
    foreach(b = c('lasso', 'step'), .combine = rbind) %do% fit_fn(h = a, f = b)
  return(out)
}
res <- foreach(i = seq_len(100), .combine = rbind) %dopar% 
  test_fn(i, n = 200, d_z = 50, rho = 0.3, k = 25, snr = 5, xzr = 1, l = 0)

#Plot
res %>%
  pivot_longer(cols = c(entner, ours), 
               names_to = 'method', values_to = 'value') %>%
  ggplot(aes(value, fill = method)) +
  geom_histogram(bins = 40, alpha = 0.75) + 
  scale_fill_npg() + 
  theme_bw() +
  facet_grid(h ~ f, scales = 'free') 

# Benchmark against naive solutions: complete discovery and ignoring Z







