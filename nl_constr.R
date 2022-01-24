### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(ranger)
library(ppcor)
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
#' @param oracle Expected output of an independence oracle, one of \code{"xy"},
#'   \code{"ci"}, or \code{"na"}.
#' 
#' @import data.table

# Data simulation function
sim_dat <- function(n, d_z, rho, sp, r2, lin_pr, oracle) {
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
  # X data
  x <- signal_x + sim_noise(signal_x, r2)
  # Identifiable?
  if (oracle == 'na') {
    shared_parents <- which(beta != 0 & gamma != 0)
    u_idx <- sample(shared_parents, size = length(shared_parents)/2)
    z <- z[, -u_idx]
    d_z <- ncol(z)
    d_u <- length(u_idx)
  } else {
    d_u <- 0
  }
  # Y data
  if (oracle %in% c('ci', 'na')) {
    y <- signal_y + sim_noise(signal_y, r2)
  } else if (oracle == 'xy') {
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
    'lin_pr' = lin_pr, 'oracle' = oracle
  )
  out <- list(
    'dat' = data.table(z, 'x' = x, 'y' = y),
    'wts' = list('beta' = beta, 'gamma' = gamma), 'params' = params
  )
  return(out)
}

################################################################################

### ENTNER METHOD ###

#' @param x First vector.
#' @param y Second vector.
#' @param z Conditioning set.
#' @param trn Training index.
#' @param tst Test index.

# GCM subroutine (Shah & Peters 2020)
gcm_test <- function(x, y, z, trn, tst) {
  rf1 <- ranger(x = z[trn, ], y = x[trn], num.trees = 50, num.threads = 1)
  rf2 <- ranger(x = z[trn, ], y = y[trn], num.trees = 50, num.threads = 1)
  eps1 <- x[tst] - predict(rf1, z[tst, ], num.threads = 1)$predictions
  eps2 <- y[tst] - predict(rf2, z[tst, ], num.threads = 1)$predictions
  nn <- length(tst)
  R <- eps1 * eps2
  R.sq <- R^2
  meanR <- mean(R)
  z_score <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p.value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p.value)
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param alpha Significance threshold for inferring dependence.
#' @param tau Other threshold for "inferring" independence.
#' @param B Number of random subset-variable pairs for testing.

# Constraint-based: for this comparison, we presume X \preceq Y
# and assume access to the true data sparsity
constr_fn <- function(sim_obj, alpha, tau, B) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  y <- dat$y
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  # Take subsets of size equal to expected size of admissible set
  fctr <- d_z + sim_obj$params$d_u
  k <- round(sim_obj$params$sp * fctr)/2
  # Evaluate R1 10x more frequently than R2
  r2_idx <- seq_len(B/10)
  # Entner function
  entner <- function(b) {
    # Sample a random subset Z_b and variable W
    z_idx <- sample(d_z, k)
    z_b <- z[, z_idx]
    j <- sample(seq_len(d_z)[-z_idx], 1)
    w <- z[, j]
    if (linear == TRUE) {
      ### RULE 1 ###
      pmat0 <- suppressWarnings(pcor(cbind(y, w, z_b))$p.value)
      p1.i <- pmat0[1, 2]
      pmat1 <- suppressWarnings(pcor(cbind(y, w, z_b, x))$p.value)
      p1.ii <- pmat1[1, 2]
      ### RULE 2 ###
      if (b %in% r2_idx) {
        pmat2 <- suppressWarnings(pcor(cbind(y, x, z_b))$p.value)
        p2.i <- pmat2[1, 2]
        pmat3 <- suppressWarnings(pcor(cbind(x, w, z_b))$p.value)
        p2.ii <- pmat3[1, 2]
      }
    } else {
      ### RULE 1 ###
      trn <- sample(n, round(0.8 * n))
      tst <- seq_len(n)[-trn]
      p1.i <- gcm_test(y, w, z_b, trn, tst)
      p1.ii <- gcm_test(y, w, cbind(z_b, x), trn, tst)
      ### RULE 2 ###
      if (b %in% r2_idx) {
        p2.i <- gcm_test(y, x, z_b, trn, tst)
        p2.ii <- gcm_test(x, w, z_b, trn, tst)
      } 
    }
    # Apply rules
    if (!b %in% r2_idx) {
      p2.i <- 0
      p2.ii <- 1
    }
    r1 <- ifelse(p1.i <= alpha & p1.ii >= tau, 1, 0)
    r2 <- ifelse(p2.i >= tau | (p2.ii <= alpha & p1.i >= tau), 1, 0)
    # Export
    out <- data.table(r1, r2)
    return(out)
  }
  # Apply Entner's rules with B random subset-variable pairs
  df <- foreach(bb = seq_len(B), .combine = rbind) %do%
    entner(bb)
  # Selecting different decision thresholds based on experimentation
  # Note -- this is very generous!
  if (df[, sum(r1) / .N] >= 1/200) {
    g <- 'xy' 
  } else if (df[seq_len(B/10), sum(r2) / .N] >= 1/5) {
    g <- 'ci'
  } else {
    g <- 'na'
  }
  out <- data.table(method = 'constr', g)
  return(out)
}

################################################################################

### PUT IT ALL TOGETHER ###

# Big ol' wrapper
big_loop <- function(sims_df, sim_id, i) {
  # Simulate data, extract ground truth
  sdf <- sims_df[s_id == sim_id]
  n_reps <- ifelse(sdf$lin_pr == 1, 1000, 500)
  sim_obj <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, sp = sdf$sp, 
                     r2 = sdf$r2, lin_pr = sdf$lin_pr, oracle = sdf$oracle)
  # Constraint function
  constr_fn(sim_obj, alpha = 0.1, tau = 0.5, B = n_reps) %>%
    mutate(s_id = sim_id, idx = i) %>%
    as.data.table(.) %>%
    return(.)
}

# Execute in parallel
sims <- expand.grid(n = c(1000, 2000, 4000),
                    oracle = c('xy', 'ci', 'na')) %>%
  mutate(sp = 0.5, d_z = 100, rho = 0.25, r2 = 2/3, lin_pr = 1/5, 
         s_id = row_number()) %>%
  as.data.table(.)

# How long will this take?
library(microbenchmark)
microbenchmark(
  n1 = big_loop(sims, 1, 1),
  n2 = big_loop(sims, 2, 1),
  n4 = big_loop(sims, 3, 1), times = 1
)

# Nonlinear:
res <- foreach(ss = sims$s_id, .combine = rbind) %:%
  foreach(ii = seq_len(50), .combine = rbind) %dopar%
  big_loop(sims, ss, ii)
fwrite(res, './results/nl_constr.csv')

