### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('./Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(pcalg)
library(RBGL)
library(glmnet)
library(randomForest)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param d_x Dimensionality of X.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param r2_params Shape parameters for the beta distribution from which 
#'   proportion of variance explained for each X variable is drawn.
#' @param wt_type Type of weight, either \code{"equal"} or \code{"unequal"}.
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param sp Average sparsity of the graph.
#' @param method Method used for generating the graph structure. Options are
#'   \code{"er"} for Erdós-Rényi and \code{"barabasi"} for Barabási-Albert.
#' @param pref Strength of preferential attachment if \code{method = "barabasi"}.
#' 

# Data simulation function
sim_dat <- function(n, d_z, d_x, rho, r2_params, wt_type, lin_pr, 
                    sp, method, pref) {
  # Simulate background variables
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  var_z <- 1 / d_z
  Sigma <- toeplitz(rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Optionally apply nonlinear transformations
  prep <- function(dat, pr) {
     out <- dat
     if (pr < 1) {
       # Pick features to transform
       n_nl <- round((1 - pr) * ncol(dat))
       if (n_nl > 0) {
         tmp <- data.table(idx = sample.int(ncol(dat), size = n_nl))
         tmp[, nl := sample(c('sq', 'sqrt', 'sftpls', 'relu'), 
                            size = n_nl, replace = TRUE)]
         out[, tmp[nl == 'sq', idx]] <- dat[, tmp[nl == 'sq', idx]]^2
         out[, tmp[nl == 'sqrt', idx]] <- sqrt(abs(dat[, tmp[nl == 'sqrt', idx]]))
         out[, tmp[nl == 'sftpls', idx]] <- log(1 + exp(dat[, tmp[nl == 'sftpls', idx]]))
         out[, tmp[nl == 'relu', idx]] <- ifelse(dat[, tmp[nl == 'relu', idx]] > 0, 
                                                 dat[, tmp[nl == 'relu', idx]], 0)
       }
     }
     return(out)
  }
  # Generate noise
  sim_noise <- function(signal, r2) {
    var_mu <- var(signal)
    var_noise <- (var_mu - r2 * var_mu) / r2
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  # Simulate graph
  m <- (1 - sp) * (d_z + d_x - 1)
  g <- randDAG(d_z + d_x, m, method = method, par1 = pref, weighted = FALSE)
  t_srt <- as.numeric(tsort(g))
  z_idx <- data.table(z = seq_len(d_z), g = t_srt[seq_len(d_z)])
  x_idx <- data.table(x = seq_len(d_x), g = t_srt[(d_z + 1):(d_z + d_x)])
  # Compute X recursively
  x <- matrix(nrow = n, ncol = d_x, 
              dimnames = list(NULL, paste0('x', seq_len(d_x))))
  for (j in seq_len(d_x)) {
    pa <- c()
    for (i in seq_len(d_z + j - 1)) {
      if (any(grepl(t_srt[d_z + j], as.numeric(g@edgeL[[t_srt[i]]]$edges)))) {
        pa <- c(pa, t_srt[i])
      }
    }
    n_pa <- length(pa)
    if (wt_type == 'equal') {
      beta <- sample(c(1, -1), size = n_pa, replace = TRUE)
    } else if (wt_type == 'unequal') {
      beta <- seq(1, 10, length.out = n_pa) * 
        sample(c(1, -1), size = n_pa, replace = TRUE)
    }
    pa_dat <- cbind(z[, z_idx[g %in% pa, z]], x[, x_idx[g %in% pa, x]])
    pa_dat <- prep(pa_dat, lin_pr)
    signal_xj <- as.numeric(pa_dat %*% beta) 
    r2_xj <- rbeta(1, r2_params[1], r2_params[2])
    x[, j] <- signal_xj + sim_noise(signal_xj, r2_xj)
  }
  # Export
  params <- list(
    'n' = n, 'd_z' = d_z, 'd_x' = d_x, 'rho' = rho, 
    'r2_params' = r2_params, 'wt_type' = wt_type, 'lin_pr' = lin_pr, 
    'sp' = sp, 'method' = method, 'pref' = pref
  )
  out <- list('dat' = data.table(z, x), 'params' = params)
  return(out)
}















