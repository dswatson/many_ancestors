### Simulations for subgraph discovery with many ancestors ###

# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries and code, register cores
library(data.table)
library(ppcor)
library(kpcalg)
library(glmnet)
library(ranger)
library(tidyverse)
library(doMC)
registerDoMC(8)
source('shah_ss.R')
source('bivariate.R')


#' @param h If X -> Y, \code{"h0"}; else if X \indep Y | Z, \code{"h1"}.
#' @param alpha Significance threshold for inferring dependence.
#' @param tau Other threshold for "inferring" independence.

# Constraint-based: for this comparison, we presume X \preceq Y
constr_fn <- function(sim_obj, h, alpha, tau) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  # Apply Entner's rules once per Z
  entnerish <- function(j) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    if (linear == TRUE) {
      pcor_pmat <- pcor(cbind(x, y, z[, j]))$p.value
      ### RULE 1 ###
      p1.i <- cor.test(y, z[, j])$p.value
      p1.ii <- pcor_pmat[2, 3]
      ### RULE 2 ###
      p2.i <- pcor_pmat[1, 2]
      p2.ii <- cor.test(x, z[, j])$p.value
    } else {
      dat <- cbind(x, y, z[, j])
      ### RULE 1 ###
      p1.i <- hsic.gamma(y, z[, j])$p.value
      p1.ii <- kernelCItest(
        x = 2, y = 3, S = 1, 
        suffStat = list(data = dat, ic.method = 'hsic.gamma')
      )
      ### RULE 2 ###
      p2.i <- kernelCItest(
        x = 1, y = 2, S = 3, 
        suffStat = list(data = dat, ic.method = 'hsic.gamma')
      )
      p2.ii <- hsic.gamma(x, z[, j])$p.value
    }
    # Apply rules
    r1 <- ifelse(p1.i <= alpha & p1.ii >= tau, 1, 0)
    r2 <- ifelse(p2.i >= tau | (p2.ii <= alpha & p1.i >= tau), 1, 0)
    # Export
    out <- data.table(j, r1, r2)
  }
  df <- foreach(jj = seq_len(ncol(z)), .combine = rbind) %do%
    entnerish(jj)
  out <- data.table(
    method = 'constr', h, 
    g = ifelse(df[, sum(r1) > sum(r2)], 'xy', 'ci')
  )
  return(out)
}

#' @param x Design matrix.
#' @param y Response vector.
#' @param trn Training index.
#' @param tst Test index.
#' @param f Function class.

# RMSE-scoring subroutine
rmse_fn <- function(x, y, trn, tst, f) {
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
  } else if (f == 'rf') {
    fit <- ranger(x = x[trn, ], y = y[trn], importance = 'impurity',
                  num.trees = 200, num.threads = 1)
    yhat_f0 <- predict(fit, data = x[tst, ], 
                       num.trees = 50, num.threads = 1)$predictions
    vimp <- data.frame('feature' = colnames(x), 
                       'imp' = fit$variable.importance) %>%
      arrange(desc(imp))
    s <- subsets(m = 10, max_d = d_z, min_d = 5, decay = 2)
    y_hat <- sapply(seq_along(s), function(k) {
      tmp_x <- x[trn, vimp$feature[seq_len(s[k])]]
      tmp_f <- ranger(x = tmp_x, y = y[trn], num.trees = 50, num.threads = 1)
      predict(tmp_f, data = x[tst, ], num.threads = 1)$predictions
    })
    y_hat <- cbind(y_hat, yhat_f0)
  }
  epsilon <- y_hat - y[tst]
  mse <- colMeans(epsilon^2)
  rmse <- sqrt(mse[which.min(mse)])
  return(rmse)
}

#' @param h If X -> Y, \code{"h0"}; else if X \indep Y | Z, \code{"h1"}.

# Score function evaluates three different DGPs
score_fn <- function(sim_obj, h) {
  # Get data
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  linear <- ifelse(sim_obj$params$lin_pr == 1, TRUE, FALSE)
  # Pick right y
  if (h == 'h0') y <- dat$y0 else y <- dat$y1
  # Pick right f
  f <- ifelse(linear, 'lasso', 'rf')
  # Split training and test
  trn <- sample(n, size = round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Score X -> Y
  rmse_xy <- rmse_fn(cbind(z, x), y, trn, tst, f)
  # Score Y -> X
  rmse_yx <- rmse_fn(cbind(z, y), x, trn, tst, f)
  # Score X ~ Y
  rmse_x <- rmse_fn(z, x, trn, tst, f)
  rmse_y <- rmse_fn(z, y, trn, tst, f)
  # Summarize
  df <- data.table(
    h, rmse = c(rmse_x + rmse_xy, rmse_y + rmse_yx, rmse_x + rmse_y),
    g = c('xy', 'yx', 'ci')
  )
  # Export
  out <- data.table(
    method = 'score', h, g = df[which.min(rmse), g]
  )
  return(out)
}

################################################################################

# Big ol' wrapper
big_loop <- function(sims_df, sim_id, i) {
  # Simulate data, extract ground truth
  sdf <- sims_df[s_id == sim_id]
  sim_obj <- sim_dat(n = sdf$n, d_z = sdf$d_z, rho = sdf$rho, 
                     sp = sdf$sp, r2 = sdf$r2, lin_pr = sdf$lin_pr)  
  # Constraint function
  df_c <- foreach(hh = c('h0', 'h1'), .combine = rbind) %do% 
    constr_fn(sim_obj, hh, alpha = 0.05, tau = 0.75)
  # Score function
  df_s <- foreach(sim_obj, hh = c('h0', 'h1'), .combine = rbind) %do%
    score_fn(sim_obj, hh)
  # Export
  out <- rbind(df_c, df_s) 
  out[, s_id := sim_id]
  out[, idx := i]
  return(out)
}

### SIMULATION GRID ###
sims <- expand.grid(
  n = c(500, 1000, 2000), d_z = c(50, 100, 200), rho = c(0, 0.5),
  sp = c(0.25, 0.5, 0.75), r2 = c(1/3, 1/2, 2/3)
)
# Linear?
if (linear == TRUE) {
  sims$lin_pr <- 1
  lab <- 'linear_sim.csv'
} else {
  sims$lin_pr <- 1/5
  lab <- 'nonlinear_sim.csv'
}
# Index, data table-ify
sims$s_id <- seq_len(nrow(sims))
sims <- as.data.table(sims)

# Compute in parallel
res <- foreach(ss = sims$s_id, .combine = rbind) %:%
  foreach(ii = seq_len(100), .combine = rbind) %dopar%
  big_loop(sims, ss, ii)




