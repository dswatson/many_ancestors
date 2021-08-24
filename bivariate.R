### Simulations for subgraph discovery with many ancestors ###

# Load libraries, register cores
library(data.table)
library(glmnet)
library(bestsubset)
library(randomForest)
library(tidyverse)
library(ggsci)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param k Number of nonzero weights from Z to children. Either a scalar, in  
#'   which case it applies to both X and Y, or a length-two vector.
#' @param snr Signal-to-noise ratio for true data generating functions. Either a 
#'   scalar, in which case it applies to both X and Y, or a length-two vector.
#' @param xzr X-to-Z ratio, i.e. the ratio of X-signal to Z-signal for 
#'   generating Y when X -> Y.
#' @param form Should structural equations be \code{"linear"} or 
#'   \code{"nonlinear"}?

# Data simulation function
sim_dat <- function(n, d_z, rho, k, snr, xzr, form) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  if (rho != 0) {
    Sigma <- toeplitz(rho^(0:(d_z - 1)))
    z <- z %*% chol(Sigma)
  } 
  dimnames(z) <- list(NULL, paste0('z', seq_len(d_z)))
  # Nonlinear transformations?
  if (form != 'linear') {
    idx <- split(sample.int(d_z), sort(seq_len(d_z) %% 5))
    names(idx) <- c('lexp', 'sqrt', 'sq', 'thresh', 'orig')
    zz <- z
    zz[, idx$lexp] <- log(1 + exp(z[, idx$lexp]))
    zz[, idx$sqrt] <- sqrt(abs(z[, idx$sqrt]))
    zz[, idx$sq] <- z[, idx$sq]^2
    zz[, idx$thresh] <- ifelse(z[, idx$thresh] > 0, z[, idx$thresh], 0)
  }
  # Draw Z weights
  sim_z_wts <- function(n, d_z, k) {
    nonzero <- sample(d_z, k)
    amplitude <- seq(1, 20, length.out = k) # Why 20?
    signs <- sample(c(1, -1), size = d_z, replace = TRUE)
    beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
    return(beta)
  }
  # Generate noise
  sim_u <- function(signal, snr) {
    noise_factor <- sqrt(var(signal) / snr)
    noise <- rnorm(n, sd = noise_factor) 
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  beta_z_to_x <- sim_z_wts(n, d_z, k)
  beta_z_to_y <- sim_z_wts(n, d_z, k)
  if (form == 'linear') {
    signal_x <- as.numeric(z %*% beta_z_to_x)
    signal_y <- as.numeric(z %*% beta_z_to_y)
  } else {
    signal_x <- as.numeric(zz %*% beta_z_to_x)
    signal_y <- as.numeric(zz %*% beta_z_to_y)
  }
  noise_x <- sim_u(signal_x, snr)
  noise_y <- sim_u(signal_y, snr)
  x <- signal_x + noise_x
  y0 <- signal_y + noise_y
  ### Alternative scenario: X -> Y
  beta_z_to_y <- sim_z_wts(n, d_z, k)
  if (form == 'linear') {
    signal_z_to_y <- as.numeric(z %*% beta_z_to_y)
    sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
    beta_x_to_y <- sigma_x_to_y / sd(x)
    signal_y <- signal_z_to_y + x * c(beta_x_to_y)
  } else {
    signal_z_to_y <- as.numeric(zz %*% beta_z_to_y)
    xx <- log(1 + exp(x))
    sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
    beta_x_to_y <- sigma_x_to_y / sd(xx)
    signal_y <- signal_z_to_y + xx * c(beta_x_to_y)
  }
  noise_y <- sim_u(signal_y, snr)
  y1 <- signal_y + noise_y
  # Export
  out <- data.table(z, 'x' = x, 'y0' = y0, 'y1' = y1)
  return(out)
}


#' @param trn_x Training set of predictors.
#' @param trn_y Training outcomes.
#' @param tst_x Test set of predictors.
#' @param tst_y Test outcomes.
#' @param f Regression function to use. 

beta_fn <- function(trn_x, trn_y, tst_x, tst_y, f) {
  if (f == 'lasso') {
    fit <- glmnet(trn_x, trn_y, intercept = FALSE)
    y_hat <- predict(fit, newx = tst_x, s = fit$lambda)
    mse <- colMeans((y_hat - tst_y)^2)
    beta <- as.numeric(coef(fit, s = fit$lambda[which.min(mse)]))[1:d_z]
  } else if (f == 'step') {
    fit <- fs(trn_x, trn_y, intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = tst_x)
    mse <- colMeans((y_hat - tst_y)^2)
    beta <- coef(fit)[1:d_z, which.min(mse)]
  } else if (f == 'rf') {
    fit <- randomForest(trn_x, trn_y, ntree = 500)
    vimp <- data.frame('feature' = colnames(trn_x), 
                           'imp' = as.numeric(importance(fit))) %>%
      arrange(desc(imp))
    # Consider 10 models, min d = 5, quadratic decay on dimensionality
    m <- unique(round(5 + ((d_z - 5) / 10^2) * seq_len(10)^2))
    err <- sapply(seq_along(m), function(i) {
      tmp_x <- trn_x[, vimp$feature[seq_len(m[i])]]
      f <- randomForest(tmp_x, trn_y, ntree = 200)
      y_hat <- predict(f, newdata = tst_x)
      mse <- mean((y_hat - tst_y)^2)
      return(mse)
    })
    k_hat <- m[which.min(err)]
    beta <- double(length = d_z)
    names(beta) <- paste0('z', seq_len(d_z))
    keep <- vimp$feature[seq_len(k_hat)]
    beta[keep] <- vimp$imp[seq_len(k_hat)]
  }
  return(beta)
}


#' @param b Simulation index.
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param k Number of nonzero weights from Z to children.
#' @param snr Signal-to-noise ratio. 
#' @param xzr X-to-Z ratio.
#' @param form Functional form for structural equations.
#' @param l Norm to use, either 0, 1, or 2.

# Big ol' wrapper
test_fn <- function(b, n, d_z, rho, k, snr, xzr, form, l) {
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Simulate data
  dat <- sim_dat(n, d_z, rho, k, snr, xzr, form)
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Compute coefficients
  fit_fn <- function(h, f) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    zy <- cbind(z, y)
    betas <- list(
      beta_fn(z[trn, ], y[trn], z[tst, ], y[tst], f),
      beta_fn(zx[trn, ], y[trn], zx[tst, ], y[tst], f),
      beta_fn(z[trn, ], x[trn], z[tst, ], x[tst], f),
      beta_fn(zy[trn, ], x[trn], zy[tst, ], x[tst], f)
    )
    v <- sapply(seq_len(4), function(i) {
      switch(l, 
             'l0' = sum(betas[[i]] != 0),
             'l1' = sum(abs(betas[[i]])),
             'l2' = sum(betas[[i]]^2)
      )})
    w0 <- v[1] - v[2]
    w1 <- v[3] - v[4]
    delta <- w0 - w1
    out <- data.table(
      'h' = h, 'f' = f, 'delta' = delta
    )
    # Export
    return(out)
  }
  out <- foreach(a = c('h0', 'h1'), .combine = rbind) %:%
    foreach(b = c('lasso', 'step'), .combine = rbind) %do% fit_fn(h = a, f = b)
  return(out)
}
res <- foreach(i = seq_len(500), .combine = rbind) %dopar% 
  test_fn(i, n = 500, d_z = 50, rho = 0.3, k = 25, snr = 3, xzr = 1, l = 'l0')

#Plot
res %>%
  pivot_longer(cols = c(entner, ours), 
               names_to = 'method', values_to = 'value') %>%
  ggplot(aes(value, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.75) + 
  scale_fill_npg() + 
  theme_bw() +
  facet_grid(h ~ f, scales = 'free') 



# Possible benchmarks: true score-based
# Revised Entner: each Z gets treated as W once, 
# see if anything satisfies R1 or R2


sample_sizes <- c(200, 500, 1000)
dim_z <- c(20, 50, 100)
auto_cor <- c(0, 0.3, 0.5)
sparsity <- c(0.1, 0.5, 0.9)
signals <- c(1, 3, 5)
x_strength <- c(0.5, 1, 2)



# To implement Entner's method with GCM or LOCO test
eps_fn <- function(trn_x, trn_y, tst_x, tst_y, f) {
  if (f == 'lasso') {
    fit <- glmnet(trn_x, trn_y, intercept = FALSE)
    y_hat <- predict(fit, newx = tst_x, s = fit$lambda)
  } else if (f == 'step') {
    fit <- fs(trn_x, trn_y, intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = tst_x)
  } else if (f == 'rf') {
    fit <- randomForest(trn_x, trn_y, ntree = 500)
    vimp <- data.frame('feature' = colnames(trn_x), 
                       'imp' = as.numeric(importance(fit))) %>%
      arrange(desc(imp))
    # Consider 10 models, min d = 5, quadratic decay on dimensionality
    m <- unique(round(5 + ((d_z - 5) / 10^2) * seq_len(10)^2))
    y_hat <- sapply(seq_along(m), function(i) {
      tmp_x <- trn_x[, vimp$feature[seq_len(m[i])]]
      f <- randomForest(tmp_x, trn_y, ntree = 200)
      y_hat <- predict(f, newdata = tst_x)
      return(y_hat)
    })
  }
  eps_mat <- y_hat - tst_y
  mse <- colMeans(eps_mat^2)
  eps <- res_mat[, which.min(mse)]
  return(eps)
}
gcm_test <- function(x, y) {
  nn <- length(x)
  R <- x * y
  R.sq <- R^2
  meanR <- mean(R)
  z <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p.value <- 2 * pnorm(abs(z), lower.tail = FALSE)
  return(p.value)
}
loco_test <- function(x, y, z) {
  eps_x <- eps_fn(z, y)
  eps_y <- eps_fn(cbind(x, z), y)
  delta <- abs(eps_x) - abs(eps_y)
  p.value <- wilcox.test(delta, alt = 'greater')$p.value
  return(p.value)
}



f_x <- lm(x ~ z)
f_y <- lm(y ~ z)
for (j in seq_len(d_z)) {
  ### Rule 1 ###
  # Rule 1(i)
  f0_w <- lm(z[, j] ~ z[, -j])
  f0_y <- lm(y ~ z[, -j])
  p1.i <- gcm_test(residuals(f0_w), residuals(f0_y))
  rule1.i <- ifelse(p1.i <= alpha, TRUE, FALSE)
  # Rule 1(ii)
  f1_w <- lm(z[, j] ~ zx[, -j])
  f1_y <- lm(y ~ zx[, -j])
  p1.ii <- cor.test(residuals(f1_w), residuals(f1_y))$p.value
  rule1.ii <- ifelse(p1.ii <= alpha, FALSE, TRUE)
  # Therefore
  rule1 <- ifelse(rule1.i & rule1.ii, TRUE, FALSE)
  
  ### Rule 2 ###
  # Rule 2(i)
  p2.i <- cor.test(residuals(f_x), residuals(f_y))$p.value
  rule2.i <- ifelse(p2.i <= alpha, FALSE, TRUE)
  # Rule 2(ii)
  f0_x <- lm(x ~ z[, -j])
  p2 <- cor.test(residuals(f0_w), residuals(f0_x))$p.value
  rule2.ii <- ifelse(p2 <= alpha, TRUE, FALSE)
  # Rule 2(iii)
  rule2.iii <- !rule1.i
  # Therefore
  rule2 <- ifelse(rule2.i | (rule2.ii & rule2.iii), TRUE, FALSE)
}





















