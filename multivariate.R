### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('./Documents/UCL/many_ancestors')

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(pcalg)
library(RBGL)
library(matrixStats)
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
#' @param r2 Proportion of variance explained for all foreground variables. 
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param method Method used for generating the graph structure. Options are
#'   \code{"er"} for Erdós-Rényi and \code{"barabasi"} for Barabási-Albert.
#' @param sp Average sparsity of the graph. Note that this must be high for 
#'   \code{method = "barabasi"} or else you'll run into errors.
#' @param pref Strength of preferential attachment if \code{method = "barabasi"}.
#' 

# Data simulation function
sim_dat <- function(n, d_z, d_x, rho, r2, lin_pr, sp, method, pref) {
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
  
  # QUESTION: SHOULD WE HAVE UNOBSERVED CONFOUNDERS???
  
  # Compute X recursively, record adjacency matrix
  x_labs <- paste0('x', seq_len(d_x))
  x <- matrix(nrow = n, ncol = d_x, dimnames = list(NULL, x_labs))
  adj_mat <- matrix(0, nrow = d_x, ncol = d_x, dimnames = list(x_labs, x_labs))
  for (j in seq_len(d_x)) {
    # Index the parents
    pa <- c()
    for (i in seq_len(d_z + j - 1)) {
      if (any(grepl(t_srt[d_z + j], as.numeric(g@edgeL[[t_srt[i]]]$edges)))) {
        pa <- c(pa, t_srt[i])
      }
    }
    # Compute Z signal with Rademacher weights
    pa_z <- prep(z[, z_idx[g %in% pa, z]], lin_pr)
    beta_z <- sample(c(1, -1), size = ncol(pa_z), replace = TRUE)
    signal_z <- as.numeric(pa_z %*% beta_z)
    # Compute X signal, if applicable
    if (any(x_idx$g %in% pa)) {
      pa_x <- as.matrix(prep(x[, x_idx[g %in% pa, x]], lin_pr))
      adj_mat[j, x_idx[g %in% pa, x]] <- 1
      causal_wt <- 1 / length(pa)
      sigma_xij <- sqrt(causal_wt * var(signal_z))
      beta_x <- sigma_xij / colSds(pa_x)
      signal_x <- as.numeric(pa_x %*% beta_x)
    } else {
      signal_x <- 0
    }
    signal_xj <- signal_z + signal_x
    # Add appropriate noise and export
    x[, j] <- signal_xj + sim_noise(signal_xj, r2)
  }
  # Export
  params <- list(
    'n' = n, 'd_z' = d_z, 'd_x' = d_x, 'rho' = rho, 'r2' = r2, 'lin_pr' = lin_pr, 
    'sp' = sp, 'method' = method, 'pref' = pref
  )
  out <- list('dat' = data.table(z, x), 'adj_mat' = adj_mat, 'params' = params)
  return(out)
}
# Note: lower triangular adj_mat means that column is a parent of row

#' @param m Number of nested models to fit.
#' @param max_x Number of predictors in largest model.
#' @param min_d Number of predictors in smallest model.
#' @param decay Exponential decay parameter

# Precompute subset sizes for RFE
subsets <- function(m, max_d, min_d, decay) {
  out <- round(min_d + ((max_d - min_d) / (m + 1)^decay) * seq_len(m + 1)^decay)
  out <- na.omit(unique(out)[seq_len(m)])
  return(out)
}


#' @param x Design matrix.
#' @param y Outcome vector.
#' @param trn Training indices.
#' @param tst Test indices.
#' @param f Regression method, either \code{"lasso"} or \code{"rf"}.

# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, trn, tst, f) {
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
    betas <- coef(fit, s = fit$lambda)[-1, ]
  } else if (f == 'rf') {
    fit <- randomForest(x[trn, ], y[trn], ntree = 200)
    vimp <- data.frame('feature' = colnames(x), 
                       'imp' = as.numeric(importance(fit))) %>%
      arrange(desc(imp))
    beta <- double(length = ncol(x))
    names(beta) <- paste0('z', seq_len(ncol(x)))
    s <- subsets(m = 10, max_d = ncol(x), min_d = 5, decay = 2)
    rfe_loop <- function(k) {
      tmp_x <- x[trn, vimp$feature[seq_len(s[k])]]
      tmp_f <- randomForest(tmp_x, y[trn], ntree = 50)
      tmp_v <- data.frame('feature' = colnames(tmp_x), 
                          'imp' = as.numeric(importance(tmp_f))) %>%
        filter(grepl('z', feature))
      beta[tmp_v$feature] <- tmp_v$imp
      out <- list('y_hat' = predict(tmp_f, newdata = x[tst, ]), 'beta' = beta)
      return(out)
    }
    rf_out <- foreach(kk = seq_along(s)) %do% rfe_loop(kk)
    y_hat <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$y_hat)
    y_hat <- cbind(y_hat, predict(fit, newdata = x[tst, ]))
    betas <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$beta)
    betas <- cbind(betas, as.numeric(importance(fit))[seq_len(ncol(x))])
  }
  epsilon <- y_hat - y[tst]
  mse <- colMeans(epsilon^2)
  betas <- betas[, which.min(mse)]
  out <- ifelse(betas == 0, 0, 1)
  return(out)
}


#' @param sim_obj Simulation object output by \code{sim_dat}.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute (de)activation rates
rate_fn <- function(sim_obj, B) {
  ### PRELIMINARIES ###
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- as.matrix(select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  # Linear or nonlinear?
  f <- ifelse(sim_obj$params$lin_pr == 1, 'lasso', 'rf')
  # Initialize adjacency matrix
  adj_mat <- matrix(NA_character_, nrow = d_x, ncol = d_x, 
                    dimnames = list(x_labs, x_labs))
  
  ### LOOP THROUGH ROUNDS ###
  t_loop <- function(...) {
    # Update dimensionality
    d_z <- ncol(z_plus)
    d_x <- ncol(x_minus)
    # Random train/test split
    trn <- sample(n, round(0.8 * n))
    tst <- seq_len(n)[-trn]
    # Fit reduced models: d_z x d_x matrix
    s0_mat <- sapply(seq_len(d_x), function(i) {
      l0(z, x[, i], trn, tst, f)
    })
    # Fit expanded models: list of length d_x, 
    # each element a (d_z + 1) x (d_x - 1) matrix
    s1_list <- lapply(seq_len(d_x), function(i) {
      m <- sapply(seq_len(d_x)[-i], function(j) {
        l0(cbind(z, x[, j]), x[, i], trn, tst, f)
      })
      colnames(m) <- paste0('x', seq_len(d_x)[-i])
      return(m)
    })
    # Apply rules: count instances of R1-R3 for all pairs
    out <- expand.grid(c = seq_len(d_x), e = seq_len(d_x), 
                       r1 = 0, r2 = 0, r3 = 0)
    out <- as.data.table(out)
    out <- out[c != e]
    for (i in seq_len(d_x)) {
      for (j in seq_len(d_x)[-i]) {
        out[c == i & e == j, r1 := 
              sum(s0_mat[, j] == 1 & s1_list[[j]][1:d_z, paste0('x', i)] == 0)]
        out[c == i & e == j, r2 := 
              sum(s0_mat[, i] == 0 & s1_list[[i]][1:d_z, paste0('x', j)] == 1)]
        out[c == i & e == j, r3 := 
              (s1_list[[i]][d_z + 1, paste0('x', j)] == 0) + 
              (s1_list[[j]][d_z + 1, paste0('x', i)] == 0)]
      }
    }
    # Now compare results for i \preceq j and j \preceq i
    # Heuristic: if r3 > 0, i ~ j. 
    # Otherwise, max(c(r1, r2)) in either direction
    # If nothing applies, then we admit our ignorance
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        tmp <- out[c %in% c(i, j) & e %in% c(i, j)]
        if (tmp[, sum(r3) > 0]) {
          adj_mat[i, j] <- '0'
        } else if (tmp[, max(r1)] > tmp[, max(r2)]) {
          adj_mat[tmp[which.max(r1), e], tmp[which.max(r1), c]] <- 'prec'
        } else if (tmp[, max(r2)] > tmp[, max(r1)]) {
          adj_mat[tmp[which.max(r2), e], tmp[which.max(r2), c]] <- 'preceq'
        }
      }
    }
    iter <- iter + 1
    
    # Need some way to update and loop back
    
    
  }
  
  

  
  
}



subdag <- function(sim_obj, maxiter) {
  ### PRELIMINARIES ###
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- as.matrix(select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  # Linear or nonlinear?
  f <- ifelse(sim_obj$params$lin_pr == 1, 'lasso', 'rf')
  # Initialize adjacency matrix
  adj_list <- list(
    matrix(NA_character_, nrow = d_x, ncol = d_x, 
           dimnames = list(x_labs, x_labs))
  )
  # Stopping parameters
  converge <- FALSE
  iter <- 0
  ### LOOP IT ###
  while(converge == FALSE & iter <= maxiter) {
    # Extract relevant adjacency matrices
    if (iter == 0) {
      adj0 <- adj1 <- adj_list[[1]]
    } else {
      adj0 <- adj_list[[iter - 1]]
      adj1 <- adj_list[[iter]]
    }
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        # Only continue if the relationship is unknown
        if (is.na(adj1[i, j]) & is.na(adj1[j, i])) {
          # Only continue if the set of nondescendants has increased since last 
          # iteration (i.e., have we learned anything new?)
          preceq_i <- which(grepl('prec', adj0[i, ]))
          preceq_j <- which(grepl('prec', adj0[j, ]))
          a0 <- intersect(preceq_i, preceq_j) 
          preceq_i <- which(grepl('prec', adj1[i, ]))
          preceq_j <- which(grepl('prec', adj1[j, ]))
          a1 <- intersect(preceq_i, preceq_j) 
          if (length(a1) > length(a0)) {
            z_t <- cbind(z, x[, a])
            d_zt <- ncol(z_t)
            # Fit reduced models
            s0 <- sapply(c(i, j), function(k) {
              l0(z_t, x[, k], trn, tst, f)
            })
            # Fit expanded models
            s1 <- sapply(c(i, j), function(k) {
              not_k <- setdiff(c(i, j), k)
              l0(cbind(z_t, x[, not_k]), x[, k], trn, tst, f)
            })
            # Apply rules
            if (any(s1[d_zt + 1, ] == 0)) {
              # If either Xi or Xj receives zero weight in the other's regression,
              # we stop there
              adj1[i, j] <- adj1[j, i] <- '0'
            } else {
              # Sum (de)activations across all nondescendants
              d_ji <- sum(s0[, 1] == 1 & s1[, 1] == 0)
              a_ji <- sum(s0[, 2] == 0 & s1[, 2] == 1)
              d_ij <- sum(s0[, 2] == 1 & s1[, 2] == 0)
              a_ij <- sum(s0[, 1] == 0 & s1[, 1] == 1)
              events <- c(d_ji, a_ji, d_ij, a_ij)
              # Heuristic: go with the max
              if (sum(events) > 0) {
                if (d_ji > max(events[-1])) {
                  adj1[j, i] <- 'prec'
                } else if (a_ji > max(events[-2]) | 
                           min(c(d_ji, a_ji)) > max(d_ij, a_ij)) {
                  adj1[j, i] <- 'preceq'
                } else if (d_ij > max(events[-3])) {
                  adj1[i, j] <- 'prec'
                } else if (a_ij > max(events[-4]) | 
                           min(c(d_ij, a_ij)) > max(d_ji, a_ji)) {
                  adj1[i, j] <- 'preceq'
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
    adj0 <- adj_list[[iter - 1]]
    adj1 <- adj_list[[iter]]
    # Check for convergence
    if (identical(adj0, adj1)) {
      converge <- TRUE
    }
  }
}











