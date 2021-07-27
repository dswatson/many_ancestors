### Simulation for Causal Discovery with High-Dimensional Ancestors ###

# Load libraries, register cores
library(data.table)
library(glmnet)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param k Number of nonzero weights from Z to child node (\code{<= d_z}).

# Fix Z weights
sim_z_wts <- function(n, d_z, k) {
  nonzero <- sample(d_z, k)
  amplitude <- seq(1, 20, length.out = k) # Why 20?
  signs <- sample(c(1, -1), size = d_z, replace = TRUE)
  beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
  return(beta)
}
beta_z_to_x <- sim_z_wts(n = 1000, d_z = 100, k = 50)
beta_z_to_y <- sim_z_wts(n = 1000, d_z = 100, k = 50)

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param snr Signal-to-noise ratio. Either a scalar value, in which case it 
#'   applies to both x and y, or a length-two vector.
#' @param xzr X-to-Z ratio, i.e. the ratio of X-signal to Z-signal when X's
#'   are not conditionally independent.

# Data simulation function
sim_dat <- function(n, d_z, rho, snr, xzr) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  if (rho > 0) {
    Sigma <- toeplitz(rho^(0:(d_z - 1)))
    z <- z %*% chol(Sigma)
  }
  dimnames(z) <- list(NULL, paste0('z', seq_len(d_z)))
  # Generate noise
  sim_u <- function(signal, snr) {
    eps <- rnorm(n)
    eps <- eps / sd(eps)
    noise_factor <- sqrt(var(signal) / snr)
    noise <- eps * c(noise_factor)
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  signal_x <- z %*% beta_z_to_x
  u_x <- sim_u(signal_x, snr)
  x0 <- signal_x + u_x
  signal_y <- z %*% beta_z_to_y
  u_y <- sim_u(signal_y, snr)
  y0 <- signal_y + u_y
  ### Alternative scenario: X -> Y
  signal_x <- z %*% beta_z_to_x
  u_x <- sim_u(signal_x, snr)
  x1 <- signal_x + u_x
  signal_z_to_y <- z %*% beta_z_to_y
  sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
  beta_x_to_y <- sigma_x_to_y / sd(x1)
  signal_y <- signal_z_to_y + x1 * c(beta_x_to_y)
  u_y <- sim_u(signal_y, snr)
  y1 <- signal_y + u_y
  # Export
  out <- list(
    'z' = z, 
    'h0' = data.frame('x' = x0, 'y' = y0), 
    'h1' = data.frame('x' = x1, 'y' = y1)
  )
  return(out)
}

# Draw lots of W's
w_fn <- function(b, n, d_z, rho, snr, xzr) {
  df <- sim_dat(n, d_z, rho, snr, xzr)
  f0 <- cv.glmnet(x = df$z, y = df$h0$y)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$x)),
                  y = df$h0$y)
  v0 <- abs(coef(f0, s = 'lambda.min')[2:(d_z + 1), ])
  v1 <- abs(coef(f1, s = 'lambda.min')[2:(d_z + 1), ])
  w0 <- sum(sign(v0 - v1))
  f0 <- cv.glmnet(x = df$z, y = df$h1$y)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$x)),
                  y = df$h1$y)
  v0 <- abs(coef(f0, s = 'lambda.min')[2:(d_z + 1), ])
  v1 <- abs(coef(f1, s = 'lambda.min')[2:(d_z + 1), ])
  w1 <- sum(sign(v0 - v1))
  out <- data.table('w0' = w0, 'w1' = w1)
  return(out)
}
res <- foreach(i = seq_len(2000), .combine = rbind) %dopar% 
  w_fn(i, n = 1000, d_z = 100, rho = 0, snr = 5, xzr = 1)

# Wrong direction
w_fn <- function(b, n, d_z, rho, snr, xzr) {
  df <- sim_dat(n, d_z, rho, snr, xzr)
  f0 <- cv.glmnet(x = df$z, y = df$h1$x)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$y)),
                  y = df$h1$x)
  v0 <- abs(coef(f0, s = 'lambda.min')[2:(d_z + 1), ])
  v1 <- abs(coef(f1, s = 'lambda.min')[2:(d_z + 1), ])
  w1 <- sum(sign(v0 - v1))
  return(w1)
}
res2 <- foreach(i = seq_len(2000), .combine = c) %dopar% 
  w_fn(i, n = 1000, d_z = 100, rho = 0, snr = 5, xzr = 1)






# Variance about zero is a function of rho! (At least holding all else constant)
# Classifying who's a parent of Y and who is not. 
# The classification rule will be a threshold of abs(w)
# Just plot ROC curve -- is AUC higher for right direction than wrong direction

# Dividing by MSE or RMSE could help with scaling





























