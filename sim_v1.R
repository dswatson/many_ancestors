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
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param max_wt Maximum amplitude of 
#' @param snr Signal-to-noise ratio. Either a scalar value, in which case it 
#'   applies to both x and y, or a length-two vector.
#' @param xzr X-to-Z ratio, i.e. the ratio of X-signal to Z-signal when X's
#'   are not conditionally independent.

# Data simulation function
sim_dat <- function(n, d_z, k, rho, snr, xzr) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  if (rho > 0) {
    Sigma <- toeplitz(rho^(0:(d_z - 1)))
    z <- z %*% chol(Sigma)
  }
  dimnames(z) <- list(NULL, paste0('z', seq_len(d_z)))
  # Generate Z-weights
  sim_z_wts <- function(d_z, k) {
    nonzero <- sample(d_z, k)
    amplitude <- seq(1, 20, length.out = k) # Why 20?
    signs <- sample(c(1, -1), size = d_z, replace = TRUE)
    beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
    return(beta)
  }
  # Generate noise
  sim_u <- function(signal, snr) {
    eps <- rnorm(n)
    eps <- eps / sd(eps)
    noise_factor <- sqrt(var(signal) / snr)
    noise <- eps * c(noise_factor)
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  beta_z_to_x <- sim_z_wts(d_z, k)
  signal_x <- z %*% beta_z_to_x
  u_x <- sim_u(signal_x, snr)
  x0 <- signal_x + u_x
  beta_z_to_y <- sim_z_wts(d_z, k)
  signal_y <- z %*% beta_z_to_y
  u_y <- sim_u(signal_y, snr)
  y0 <- signal_y + u_y
  ### Alternative scenario: X -> Y
  beta_z_to_x <- sim_z_wts(d_z, k)
  signal_x <- z %*% beta_z_to_x
  u_x <- sim_u(signal_x, snr)
  x1 <- signal_x + u_x
  beta_z_to_y <- sim_z_wts(d_z, k)
  signal_z_to_y <- z %*% beta_z_to_y
  sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
  beta_x_to_y <- sigma_x_to_y / sd(x1)
  signal_y <- signal_z_to_y + x1 * c(beta_x_to_y)
  u_y <- sim_u(signal_y, snr)
  y1 <- signal_y + u_y
  # Export
  out <- list(
    'data' = list(
      'z' = z, 
      'h0' = data.frame('x' = x0, 'y' = y0), 
      'h1' = data.frame('x' = x1, 'y' = y1)
    ),
    'beta' = list(
      'h0' = c(beta_z_to_x, beta_z_to_y), 
      'h1' = c(beta_z_to_x, beta_z_to_y, beta_x_to_y)
    )
  )
  return(out)
}

# The sparser the Z to X connections, the greater differences we should expect
# between pa(x) and pa(y). Do we want both x and y to have as much 
# Z signal as each other? Should the X_j -> X_k signal swamp the Z signal?


# Try to learn!
df <- sim_dat(n = 1000, d_z = 100, k = 50, rho = 0.1, snr = 5, xzr = 2)

# Quartet of models
f0 <- cv.glmnet(x = as.matrix(cbind(df$data$z, df$data$h0$x[sample.int(n)])), 
                y = df$data$h0$y)
f1 <- cv.glmnet(x = as.matrix(cbind(df$data$z, df$data$h0$x)),
                y = df$data$h0$y)
b0 <- coef(f0, s = 'lambda.min')[2:(d_z + 1), ]
b1 <- coef(f1, s = 'lambda.min')[2:(d_z + 1), ] 
w <- abs(b0) - abs(b1)

































