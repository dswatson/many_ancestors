### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(lightgbm)
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
    'lin_pr' = lin_pr, 'g' = g
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
#' @param val Validation index.
#' @param tst Test index.
#' @param prms List of model parameters.
#' 

# Generalized covariance measure (Shah & Peters, 2020). 
# Tests conditional independence of x and y given z using gradient boosting.
gcm_test <- function(x, y, z, trn, val, tst, prms) {
  # Model 1
  d1_trn <- lgb.Dataset(z[trn, ], label = y[trn])
  d1_val <- lgb.Dataset.create.valid(d1_trn, z[val, ], label = y[val])
  f1 <- lgb.train(params = prms, data = d1_trn, valids = list(val = d1_val),
                  nrounds = 3000, early_stopping_rounds = 10, verbose = 0)
  eps1 <- y[tst] - predict(f1, z[tst, ])
  # Model 2
  d2_trn <- lgb.Dataset(z[trn, ], label = x[trn])
  d2_val <- lgb.Dataset.create.valid(d2_trn, z[val, ], label = x[val])
  f2 <- lgb.train(params = prms, data = d2_trn, valids = list(val = d2_val),
                  nrounds = 3000, early_stopping_rounds = 10, verbose = 0)
  eps2 <- x[tst] - predict(f2, z[tst, ])
  # Inference
  nn <- length(tst)
  R <- eps1 * eps2
  R.sq <- R^2
  meanR <- mean(R)
  z_score <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p_value)
}


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param k Size of expected admissible set.
#' @param alpha Significance threshold for inferring dependence.
#' @param tau Other threshold for "inferring" independence.
#' 

# Constraint-based: for this comparison, we presume X \preceq Y
# and assume access to the true data sparsity
constr_fn <- function(z, x, y, linear, k, alpha, tau) {
  # Preliminaries
  n <- nrow(z)
  d_z <- ncol(z)
  if (linear) {
    B <- 1000
    r1_thresh <- 1/200
  } else {
    B <- 500
    r1_thresh <- 0
    prms <- list(
      objective = 'regression', max_depth = 1, 
      bagging.fraction = 0.5, feature_fraction = 0.8, 
      num_threads = 1, force_col_wise = TRUE
    )
  }
  # Evaluate R1 10x more frequently than R2
  r2_idx <- seq_len(B/10)
  # Entner function
  entner <- function(b) {
    # Sample a random subset Z_b and variable W
    z_idx <- sample(d_z, k)
    z_b <- z[, z_idx]
    j <- sample(seq_len(d_z)[-z_idx], 1)
    w <- z[, j]
    if (linear) {
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
      # Train/validation/test split
      trn <- sample(n, round(0.7 * n))
      val <- sample(setdiff(seq_len(n), trn), round(0.15 * n))
      tst <- seq_len(n)[-c(trn, val)]
      ### RULE 1 ###
      p1.i <- gcm_test(w, y, z_b, trn, val, tst, prms)
      p1.ii <- gcm_test(w, y, cbind(z_b, x), trn, val, tst, prms)
      ### RULE 2 ###
      if (b %in% r2_idx) {
        p2.i <- gcm_test(x, y, z_b, trn, val, tst, prms)
        p2.ii <- gcm_test(w, x, z_b, trn, val, tst, prms)
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
  if (df[, sum(r1)] > r1_thresh) {
    g_hat <- 'xy' 
  } else if (df[seq_len(B/10), sum(r2) / .N] >= 1/5) {
    g_hat <- 'ci'
  } else {
    g_hat <- 'na'
  }
  out <- data.table(method = 'constr', g_hat)
  return(out)
}

################################################################################

### PUT IT ALL TOGETHER ###

# Initialize
out <- data.table(
  method = NA, g_hat = NA, n = NA, g = NA, idx = NA
)
res_file <- './results/lin_constr.rds'
saveRDS(out, res_file)

# Big ol' wrapper
big_loop <- function(n, g, i) {
  # Simulate data
  sim_obj <- sim_dat(n = n, d_z = 100, rho = 0.25, sp = 0.5, 
                     r2 = 2/3, lin_pr = 1/5, g = g)
  # Extract data
  dat <- sim_obj$dat
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  y <- dat$y
  k <- round(d_z/4)
  # Constraint function
  new <- constr_fn(z, x, y, linear = TRUE, k, alpha = 0.1, tau = 0.5) %>%
    mutate('n' = n, 'g' = g, 'idx' = i)
  # Import, export
  old <- readRDS(res_file)
  out <- na.omit(rbind(old, new))
  saveRDS(out, res_file)
}

# Linear
foreach(ii = seq_len(100)) %:%
  foreach(nn = c(5000, 1e4, 2e4)) %:%
  foreach(gg = c('xy', 'ci', 'na')) %dopar%
  big_loop(nn, gg, ii)



ggplot(df, aes(rule, rate, fill = rule)) + 
  geom_boxplot() + 
  scale_fill_npg() + 
  theme_bw() + 
  facet_grid(g ~ n)


# What we want to see: highest incidence of R1 in g = 'xy', 
# highest incidence of R2 in g = 'ci'. Clearer preferences in larger sample
# sizes. Interestingly, in the linear case, R2 fires much more frequently
# than R1. Not so here, although that could obviously change with different
# values for alpha and tau.

# Experiments are suggesting that perhaps we need very different values of B
# for different sample sizes? More runs for smaller n. R1 dominates in 
# frequency when g = 'xy', but R2 doesn't dominate for g = 'ci'. Margin of 
# victory varies with sample size. I don't think we've got time for tons more
# runs, so if power is weak at n = 5000, then so be it. 



























