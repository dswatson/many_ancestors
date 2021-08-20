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
#' @param s Number of nonzero weights from Z to child node (\code{<= d_z}).

# Fix Z weights
sim_z_wts <- function(n, d_z, s) {
  nonzero <- sample(d_z, s)
  amplitude <- seq(1, 20, length.out = s) # Why 20?
  signs <- sample(c(1, -1), size = d_z, replace = TRUE)
  beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
  return(beta)
}
beta_z_to_x <- sim_z_wts(n = 500, d_z = 20, s = 10)
beta_z_to_y <- sim_z_wts(n = 500, d_z = 20, s = 10)



#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param snr Signal-to-noise ratio. Either a scalar value, in which case it 
#'   applies to both X and Y, or a length-two vector.
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



#' @param b Simulation index.
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param snr Signal-to-noise ratio. Either a scalar value, in which case it 
#'   applies to both X and Y, or a length-two vector.
#' @param xzr X-to-Z ratio, i.e. the ratio of X-signal to Z-signal when X's
#'   are not conditionally independent.
#'   

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

# Big ol' wrapper
test_fn <- function(b, n, d_z, rho, snr, xzr, l = 0) {
  # Simulate data (with fixed weights)
  trn <- sim_dat(0.8 * n, d_z, rho, snr, xzr)
  tst <- sim_dat(0.2 * n, d_z, rho, snr, xzr)
  # Linear weights
  fit_fn <- function(h, f) {
    trn_zx <- as.matrix(cbind(trn$z, trn[[h]]$x))
    tst_zx <- as.matrix(cbind(tst$z, tst[[h]]$x))
    trn_zy <- as.matrix(cbind(trn$z, trn[[h]]$y))
    tst_zy <- as.matrix(cbind(tst$z, tst[[h]]$y))
    betas <- list(
      beta_fn(trn$z, trn[[h]]$y, tst$z, tst[[h]]$y, f),
      beta_fn(trn_zx, trn[[h]]$y, tst_zx, tst[[h]]$y, f),
      beta_fn(trn$z, trn[[h]]$x, tst$z, tst[[h]]$x, f),
      beta_fn(trn_zy, trn[[h]]$x, tst_zy, tst[[h]]$x, f)
    )
    if (l == 0) {
      l <- sapply(seq_len(4), function(i) {
        sum(betas[[i]] != 0)
      })
    } else if (l == 1) {
      l <- sapply(seq_len(4), function(i) {
        sum(abs(betas[[i]]))
      })
    } else if (l == 2) {
      l <- sapply(seq_len(4), function(i) {
        sum(betas[[i]]^2)
      })
    }
    # Entner's rules
    r1_x <- l[1] - l[2]
    r2_x <- l[3] - l[1]
    x_score <- r1_x - r2_x
    r1_y <- l[3] - l[4]
    r2_y <- l[1] - l[3]
    y_score <- r1_y - r2_y
    delta1 <- x_score - y_score
    # Our method
    w0 <- l[1] - l[2]
    w1 <- l[3] - l[4]
    delta2 <- w0 - w1
    out <- data.table(
      'h' = h, 'f' = f, 'entner' = delta1, 'ours' = delta2
    )
    # Export
    return(out)
  }
  out <- foreach(a = c('h0', 'h1'), .combine = rbind) %:%
    foreach(b = c('lasso', 'step'), .combine = rbind) %do% fit_fn(h = a, l = b)
  return(out)
}
res <- foreach(i = seq_len(1000), .combine = rbind) %dopar% 
  test_fn(i, n = 500, d_z = 20, rho = 0.3, snr = 5, xzr = 1)


# Consider taking different norms?
# Benchmark against naive solutions: complete discovery and ignoring Z

#Plot
res %>%
  pivot_longer(cols = c(entner, ours), 
               names_to = 'method', values_to = 'value') %>%
  ggplot(aes(value, fill = method)) +
  geom_histogram(bins = 40, alpha = 0.75) + 
  scale_fill_npg() + 
  theme_bw() +
  facet_grid(h ~ f, scales = 'free') 









