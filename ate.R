### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(lightgbm)
library(tidyverse)
library(doMC)
registerDoMC(16)

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
    z_out <- apply(z_out, 2, function(z) z - mean(z))
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
  x <- rbinom(n, 1, plogis(signal_x))
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
    signal_y <- signal_y + x
    y <- signal_y + sim_noise(signal_y, r2)
  }
  y <- y - mean(y)
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

### ESTIMATE CAUSAL EFFECTS ### 

# Big ol' wrapper
ate_loop <- function(pve, i) { # or g = 'na'?
  # Simulate data
  sim_obj <- sim_dat(n = 1e4, d_z = 100, rho = 1/4, sp = 1/2, 
                     r2 = pve, lin_pr = 1/5, g = 'xy')
  # Extract data
  dat <- sim_obj$dat
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  y <- dat$y
  zx <- cbind(z, x)
  trn <- sample(n, round(n * 0.8))
  tst <- seq_len(n)[-trn]
  # Fit models (these will be recycled)
  f_hat <- fit_fn(z[trn, ], x[trn], f, prms)
  g_hat <- fit_fn(z[trn, ], y[trn], f, prms)
  h_hat <- fit_fn(zx[trn, ], y[trn], f, prms)
  # Prediction vectors
  pi_hat <- predict(f_hat, z[tst, ])
  mu_hat <- predict(g_hat, z[tst, ])
  al_hat <- predict(h_hat, zx[tst, ])
  # DML
  eps_x <- x[tst] - pi_hat
  eps_y <- y[tst] - mu_hat
  ate_dml <- as.numeric(coef(lm(eps_y ~ 0 + eps_x)))
  # IPW
  ate_ipw <- mean((x[tst] * y[tst])/pi_hat - ((1 - x[tst]) * y[tst])/(1 - pi_hat))
  # TMLE
  trgt <- x[tst]/pi_hat - (1 - x[tst])/(1 - pi_hat)
  trgt1 <- 1 / pi_hat
  trgt0 <- -1 / (1 - pi_hat)
  d <- data.frame(pseudo_y = y[tst] - al_hat, trgt)
  delta <- as.numeric(coef(lm(pseudo_y ~ 0 + trgt, data = d)))
  y1_star <- y1_hat + delta * trgt1
  y0_star <- y0_hat + delta * trgt0
  ate_tmle <- mean(y1_star - y0_star)
  # Export
  out <- data.table(
    Method = c('DML', 'IPW', 'TMLE'),
    ATE = c(ate_dml, ate_ipw, ate_tmle),
    rho = ac, r2 = pve, idx = i
  )
  return(out)
}
df <- foreach(rr = c(1/3, 1/2, 2/3), .combine = rbind) %:%
  foreach(ii = seq_len(1000), .combine = rbind) %dopar%
  ate_loop(rr, ii)
saveRDS(df, './results/ate.rds')

# Polish
df[r2 == 1/3, SNR := 'SNR = 0.5']
df[r2 == 1/2, SNR := 'SNR = 1']
df[r2 == 2/3, SNR := 'SNR = 2']

# Plot
ggplot(tmp, aes(Estimator, ATE, fill = Estimator)) + 
  geom_violin(alpha = 0.85) + 
  scale_fill_d3() +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') +
  theme_bw() +  
  theme(legend.position = 'bottom',
        legend.box.background = element_rect(colour = 'black')) +
  facet_wrap(~ SNR)





