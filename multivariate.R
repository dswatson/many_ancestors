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
    signal_zxj <- as.numeric(pa_z %*% beta_z)
    # Compute X signal, if applicable
    pa_x <- prep(x[, x_idx[g %in% pa, x]], lin_pr)
    if (ncol(pa_x) >= 1) {
      adj_mat[x_idx[g %in% pa, x], j] <- 1
      causal_wt <- 1 / length(pa)
      sigma_xij <- sqrt(causal_wt * var(signal_zxj))
      beta_x <- sigma_xij / colSds(pa_x)
      signal_xij <- as.numeric(pa_x %*% beta_x)
    } else {
      signal_xij <- 0
    }
    signal_xj <- signal_zxj + signal_xij
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
  # Linear or nonlinear?
  f <- ifelse(sim_obj$params$lin_pr == 1, 'lasso', 'rf')
  
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
    # Apply them rules
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
    # I say: if r3 > 0, it's i ~ j
    # Otherwise, max(c(r1, r2)) unless it's goose eggs across the board
    # Need some way to update and loop back
    
    
  }
  
  

  
  
}














