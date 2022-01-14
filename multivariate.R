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


#' @param df Table of (de)activation rates.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute consistency lower bound
lb_fn <- function(df, B) {
  # Loop through thresholds
  lies <- function(tau) {
    # Internal consistency
    df[, dji := ifelse(drji >= tau, 1, 0)]
    df[, aji := ifelse(arji >= tau, 1, 0)]
    df[, dij := ifelse(drij >= tau, 1, 0)]
    df[, aij := ifelse(arij >= tau, 1, 0)]
    df[, int_err := ifelse((dji + aji > 1) | (dij + aij > 1), 1, 0)]
    int_err <- sum(df$int_err)
    # External consistency
    sum_ji <- df[, sum(dji + aji)]
    sum_ij <- df[, sum(dij + aij)]
    ext_err <- ifelse(min(c(sum_ji, sum_ij)) > 0, 1, 0)
    # Export
    out <- data.table('tau' = tau, 'int_err' = int_err, 'ext_err' = ext_err)
  }
  lie_df <- foreach(tt = seq_len(2 * B) / (2 * B), .combine = rbind) %do% 
    lies(tt)
  # Compute minimal thresholds
  min_int <- lie_df[int_err == 0, min(tau)]
  min_ext <- lie_df[ext_err == 0, min(tau)]
  min_two <- lie_df[int_err == 0 & ext_err == 0, min(tau)] # It's always ext
  # Export
  return(min_two)
}


#' @param df Table of (de)activation rates.
#' @param lb Consistency lower bound, as computed by \code{lb_fn}.
#' @param order Causal order of interest, either \code{"ij"} or \code{"ji"}.
#' @param rule Inference rule, either \code{"R1"} or \code{"R2"}.
#' @param B Number of complementary pairs to draw for stability selection.

# Infer causal direction using stability selection
ss_fn <- function(df, lb, order, rule, B) {
  # Find the right rate
  if (order == 'ji' & rule == 'R1') {
    r <- df[, drji]
  } else if (order == 'ji' & rule == 'R2') {
    r <- df[, arji]
  } else if (order == 'ij' & rule == 'R1') {
    r <- df[, drij]
  } else if (order == 'ij' & rule == 'R2') {
    r <- df[, arij]
  } 
  # Stability selection parameters
  q <- sum(r)
  theta <- q / length(r)
  ub <- minD(theta, B) * sum(r <= theta)
  tau <- seq_len(2 * B) / (2 * B)
  # Do any features exceed the upper bound?
  dat <- data.frame(tau, err_bound = ub) %>%
    filter(tau > lb) %>%
    rowwise() %>%
    mutate(detected = sum(r >= tau)) %>% 
    ungroup(.) %>%
    mutate(surplus = ifelse(detected > err_bound, 1, 0))
  # Export
  out <- data.table(
    'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}


#' @param sim_obj Simulation object as computed by \code{sim_dat}.
#' @param maxiter Maximum number of iterations to loop through if convergence
#'   is elusive.
#' @param B Number of complementary pairs to draw for stability selection.

subdag <- function(sim_obj, maxiter = 100, B = 50) {
  ### PRELIMINARIES ###
  # Get data, hyperparameters, train/test split
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- as.matrix(select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  f <- ifelse(sim_obj$params$lin_pr == 1, 'lasso', 'rf')
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
            df <- foreach(bb = seq_len(B), .combine = rbind) %dopar%
              sub_loop(bb, i, j, a1)
            # Compute rates
            df[, disr := sum(dis) / .N]
            if (df$disr[1] >= 0.25) { # Totally arbitrary threshold?
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

# Evaluate performance
sim <- sim_dat(n = 2000, d_z = 100, d_x = 10, rho = 0, r2 = 2/3, lin_pr = 1,
               sp = 0.5, method = 'er', pref = 1)
ahat <- subdag(sim)

# Sensitivity and specificity
df <- data.table(y = as.numeric(sim$adj_mat), yhat = as.numeric(ahat))
df[y == 0, sum(yhat, na.rm = TRUE) / .N]     # False positive rate
df[y == 1, sum(yhat, na.rm = TRUE) / .N]     # Power




# TODO:
# -Add unobserved confounders



