# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(geneNetBP)
library(matrixStats)
library(glmnet)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import data
data(hdl)

# Create noisy SNPs
n <- nrow(hdl)
probs <- runif(50)
snps <- matrix(rbinom(n * 50, 1, prob = probs), ncol = 50, byrow = TRUE,
               dimnames = list(NULL, paste0('snp', seq_len(50))))

# Background variables
z <- hdl %>%
  select(starts_with('c')) %>%
  cbind(snps)

# Foreground variables
x <- hdl %>%
  select(-starts_with('c'), -HDL)

# Subdag discovery algorithm
subdag <- function(z, x, gamma = 0.25, maxiter = 100, B = 50) {
  ### PRELIMINARIES ###
  # Get data, hyperparameters, train/test split
  n <- nrow(z)
  d_z <- ncol(z)
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  f <- 'lasso'
  # Initialize
  adj_list <- list(
    matrix(NA_real_, nrow = d_x, ncol = d_x, 
           dimnames = list(xlabs, xlabs))
  )
  converged <- FALSE
  iter <- 0
  ### LOOP IT ###
  while(converged == FALSE & iter <= maxiter) {
    # Extract relevant adjacency matrices
    if (iter == 0) {
      adj0 <- adj1 <- adj_list[[1]]
    } else {
      adj0 <- adj_list[[iter]]
      adj1 <- adj_list[[iter + 1]]
    }
    # Subsampling loop
    sub_loop <- function(b, i, j, a1) {
      z_t <- cbind(z, x[, a1])
      d_zt <- ncol(z_t)
      # Take complementary subsets
      a_set <- sample(n, round(0.5 * n))
      a_trn <- sample(a_set, round(0.8 * length(a_set)))
      a_tst <- setdiff(a_set, a_trn)
      b_set <- seq_len(n)[-a_set]
      b_trn <- sample(b_set, round(0.8 * length(b_set)))
      b_tst <- setdiff(b_set, b_trn)
      # Fit reduced models
      s0 <- sapply(c(i, j), function(k) {
        c(l0(z_t, x[, k], a_trn, a_tst, f), 
          l0(z_t, x[, k], b_trn, b_tst, f))
      })
      # Fit expanded models
      s1 <- sapply(c(i, j), function(k) {
        not_k <- setdiff(c(i, j), k)
        c(l0(cbind(z_t, x[, not_k]), x[, k], a_trn, a_tst, f),
          l0(cbind(z_t, x[, not_k]), x[, k], b_trn, b_tst, f))
      })
      # Record disconnections and (de)activations
      dis_a <- any(s1[d_zt + 1, ] == 0)
      dis_b <- any(s1[2 * (d_zt + 1), ] == 0)
      dis <- rep(c(dis_a, dis_b), each = d_zt)
      d_ji <- s0[, 1] == 1 & s1[seq_len(d_zt), 1] == 0
      a_ji <- s0[, 2] == 0 & s1[seq_len(d_zt), 2] == 1
      d_ij <- s0[, 2] == 1 & s1[seq_len(d_zt), 2] == 0
      a_ij <- s0[, 1] == 0 & s1[seq_len(d_zt), 1] == 1
      # Export
      out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_zt), i, j,
                        z = rep(colnames(z_t), times = 2),
                        dis, d_ji, a_ji, d_ij, a_ij)
      return(out)
    }
    # Pairwise test loop
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        # Only continue if relationship is unknown
        if (is.na(adj1[i, j]) & is.na(adj1[j, i])) { 
          preceq_i <- which(adj0[i, ] > 0)
          preceq_j <- which(adj0[j, ] > 0)
          a0 <- intersect(preceq_i, preceq_j) 
          preceq_i <- which(adj1[i, ] > 0)
          preceq_j <- which(adj1[j, ] > 0)
          a1 <- intersect(preceq_i, preceq_j) 
          # Only continue if the set of nondescendants has increased since last 
          # iteration (i.e., have we learned anything new?)
          if (iter == 0 | length(a1) > length(a0)) {
            df <- foreach(bb = seq_len(B), .combine = rbind) %do%
              sub_loop(bb, i, j, a1)
            # Compute rates
            df[, disr := sum(dis) / .N]
            if (df$disr[1] >= gamma) { # Totally arbitrary threshold?
              adj1[i, j] <- adj1[j, i] <- 0
            } else {
              df[, drji := sum(d_ji) / .N, by = z]
              df[, arji := sum(a_ji) / .N, by = z]
              df[, drij := sum(d_ij) / .N, by = z]
              df[, arij := sum(a_ij) / .N, by = z]
              df <- unique(df[, .(i, j, z, disr, drji, arji, drij, arij)])
              # Consistent lower bound
              lb <- lb_fn(df, B)
              # Stable upper bound
              out <- foreach(oo = c('ji', 'ij'), .combine = rbind) %:%
                foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
                ss_fn(df, lb, oo, rr, B)
              # Update adjacency matrix
              if (sum(out$decision) == 1) {
                if (out[decision == 1, order == 'ji' & rule == 'R1']) {
                  adj1[i, j] <- 1
                } else if (out[decision == 1, order == 'ji' & rule == 'R2']) {
                  adj1[i, j] <- 0.5
                } else if (out[decision == 1, order == 'ij' & rule == 'R1']) {
                  adj1[j, i] <- 1
                } else if (out[decision == 1, order == 'ij' & rule == 'R2']) {
                  adj1[j, i] <- 0.5
                }
              } else if (sum(out$decision == 2)) {
                if (out[order == 'ji', sum(decision) == 2]) {
                  adj1[i, j] <- 0.5
                } else if (out[order == 'ij', sum(decision) == 2]) {
                  adj1[j, i] <- 0.5
                }
              }
            }
          }
        } 
      }
    }
    # Store that iteration's adjacency matrix
    iter <- iter + 1
    adj_list <- append(adj_list, list(adj1))
    # Check for convergence
    if (identical(adj0, adj1)) {
      converged <- TRUE
    }
  }
  # Export final adjacency matrix
  adj_mat <- adj_list[[length(adj_list)]]
  return(adj_mat)
}









