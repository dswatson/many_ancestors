### Simulations for subgraph discovery with many ancestors ###

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

# Big ol' wrapper
test_fn <- function(b, n, d_z, rho, snr, xzr) {
  # Simulate data (with fixed weights)
  df <- sim_dat(n, d_z, rho, snr, xzr)
  # Linear weights
  fit_fn <- function(h, l) {
    if (l == 1) {
      x1 <- as.matrix(cbind(df$z, df[[h]]$x))
      x2 <- as.matrix(cbind(df$z, df[[h]]$y))
      f1 <- cv.glmnet(x = df$z, y = df[[h]]$y)
      f2 <- cv.glmnet(x = x1, y = df[[h]]$y)
      f3 <- cv.glmnet(x = df$z, y = df[[h]]$x)
      f4 <- cv.glmnet(x = x2, y = df[[h]]$x)
      models <- list(f1, f2, f3, f4)
      betas <- lapply(seq_len(4), function(i) {
        as.numeric(coef(models[[i]], s = 'lambda.min')[2:(d_z + 1)])
      })
    } else if (l == 0) {
      tmp <- data.frame(x = df[[h]]$x, y = df[[h]]$y, df$z)
      max_f1 <- lm(y ~ ., data = select(tmp, -x))
      max_f2 <- lm(y ~ ., data = tmp)
      max_f3 <- lm(x ~ ., data = select(tmp, -y))
      max_f4 <- lm(x ~ ., data = tmp)
      f1 <- step(max_f1, direction = 'backward', trace = 0, k = log(n))
      f2 <- step(max_f2, direction = 'backward', trace = 0, k = log(n))
      f3 <- step(max_f3, direction = 'backward', trace = 0, k = log(n))
      f4 <- step(max_f4, direction = 'backward', trace = 0, k = log(n))
      models <- list(f1, f2, f3, f4)
      betas <- lapply(seq_len(4), function(i) {
        out <- double(length = d_z)
        names(out) <- paste0('z', seq_len(d_z))
        b <- coef(models[[i]])[-1]
        out[names(b)] <- b
        return(out)
      })
    }
    l0 <- sapply(seq_len(4), function(i) {
      sum(betas[[i]] != 0)
    })
    # Entner's rules
    r1_x <- l0[1] - l0[2]
    r2_x <- l0[3] - l0[1]
    x_score <- r1_x - r2_x
    r1_y <- l0[3] - l0[4]
    r2_y <- l0[1] - l0[3]
    y_score <- r1_y - r2_y
    delta1 <- x_score - y_score
    # Our method
    w0 <- l0[1] - l0[2]
    w1 <- l0[3] - l0[4]
    delta2 <- w0 - w1
    out <- data.table(
      'h' = h, 'l' = l, 'entner' = delta1, 'ours' = delta2
    )
    # Export
    return(out)
  }
  out <- foreach(a = c('h0', 'h1'), .combine = rbind) %:%
    foreach(b = c(0, 1), .combine = rbind) %do% fit_fn(h = a, l = b)
  return(out)
}
res <- foreach(i = seq_len(1000), .combine = rbind) %dopar% 
  test_fn(i, n = 1000, d_z = 100, rho = 0.1, snr = 5, xzr = 1)








