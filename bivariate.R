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
#' @param sp Proportion of nonzero weights from Z to children. Either a 
#'   scalar, in which case it applies to both X and Y, or a length-two vector.
#' @param snr Signal-to-noise ratio for true data generating functions. Either a 
#'   scalar, in which case it applies to both X and Y, or a length-two vector.
#' @param xzr X-to-Z ratio, i.e. the ratio of X-signal to Z-signal for 
#'   generating Y when X -> Y.
#' @param form Should structural equations be \code{"linear"} or 
#'   \code{"nonlinear"}?

# Data simulation function
sim_dat <- function(n, d_z, rho, sp, snr, xzr, form) {
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
    names(idx) <- c('lexp', 'sqrt', 'sq', 'relu', 'orig')
    zz <- z
    zz[, idx$lexp] <- log(1 + exp(z[, idx$lexp]))
    zz[, idx$sqrt] <- sqrt(abs(z[, idx$sqrt]))
    zz[, idx$sq] <- z[, idx$sq]^2
    zz[, idx$relu] <- ifelse(z[, idx$relu] > 0, z[, idx$relu], 0)
  }
  # Draw weights for Z
  k <- round(sp * d_z)
  sim_wts <- function(n, d_z, k) {
    nonzero <- sample(d_z, k)
    amplitude <- seq(1, 20, length.out = k) # Why 20?
    signs <- sample(c(1, -1), size = d_z, replace = TRUE)
    beta <- amplitude * (seq_len(d_z) %in% nonzero) / sqrt(n) * signs
    return(beta)
  }
  # Generate noise
  sim_noise <- function(signal, snr) {
    noise_factor <- sqrt(var(signal) / snr)
    noise <- rnorm(n, sd = noise_factor) 
    return(noise)
  }
  ### Null scenario: X \indep Y | Z
  beta_z_to_x <- sim_wts(n, d_z, k)
  beta_z_to_y <- sim_wts(n, d_z, k)
  if (form == 'linear') {
    signal_x <- as.numeric(z %*% beta_z_to_x)
    signal_y <- as.numeric(z %*% beta_z_to_y)
  } else {
    signal_x <- as.numeric(zz %*% beta_z_to_x)
    signal_y <- as.numeric(zz %*% beta_z_to_y)
  }
  x <- signal_x + sim_noise(signal_x, snr)
  y0 <- signal_y + sim_noise(signal_y, snr)
  ### Alternative scenario: X -> Y
  signal_z_to_y <- signal_y
  sigma_x_to_y <- sqrt(xzr * var(signal_z_to_y))
  beta_x_to_y <- sigma_x_to_y / sd(x)
  signal_y <- signal_z_to_y + x * beta_x_to_y
  y1 <- signal_y + sim_noise(signal_y, snr)
  # Export
  out <- list('dat' = data.table(z, 'x' = x, 'y0' = y0, 'y1' = y1),
              'wts' = data.table('x' = beta_z_to_x, 'y' = beta_z_to_y))
  return(out)
}

#' @param m Number of nested models to fit.
#' @param max_x Number of predictors in largest model.
#' @param min_d Number of predictors in smallest model.
#' @param decay Exponential decay parameter

# Precompute subset sizes for RFE
subsets <- function(m, max_d, min_d, decay) {
  unique(round(min_d + ((max_d - min_d) / m^decay) * seq_len(m)^decay))
}


#' @param trn_x Training set of predictors.
#' @param trn_y Training outcomes.
#' @param tst_x Test set of predictors.
#' @param tst_y Test outcomes.
#' @param f Regression function to use. 
#' @param d_z Dimensionality of Z.

# Fit regressions, return estimated weights and residuals
f_fn <- function(trn_x, trn_y, tst_x, tst_y, f, d_z) {
  if (f == 'lasso') {
    fit <- glmnet(trn_x, trn_y, intercept = FALSE)
    y_hat <- predict(fit, newx = tst_x, s = fit$lambda)
    betas <- coef(fit, s = fit$lambda)[2:(d_z + 1), ]
  } else if (f == 'step') {
    fit <- fs(trn_x, trn_y, intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = tst_x)
    betas <- coef(fit)[1:d_z, ]
  } else if (f == 'rf') {
    fit <- randomForest(trn_x, trn_y, ntree = 500)
    vimp <- data.frame('feature' = colnames(trn_x), 
                           'imp' = as.numeric(importance(fit))) %>%
      arrange(desc(imp))
    beta <- double(length = d_z)
    names(beta) <- paste0('z', seq_len(d_z))
    s <- subsets(m = 10, max_d = d_z, min_d = 5, decay = 2)
    rfe_loop <- function(k) {
      tmp_x <- trn_x[, vimp$feature[seq_len(s[k])]]
      tmp_f <- randomForest(tmp_x, trn_y, ntree = 200)
      beta[colnames(tmp_x)] <- as.numeric(importance(tmp_f))
      out <- list('y_hat' = predict(tmp_f, newdata = tst_x), 'beta' = beta)
      return(out)
    }
    rf_out <- foreach(kk = seq_along(s)) %do% rfe_loop(kk)
    y_hat <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$y_hat)
    betas <- sapply(seq_along(rf_out), function(k) rf_out[[k]]$beta)
  }
  # Find best fit, export results
  epsilon <- y_hat - tst_y
  mse <- colMeans(epsilon^2)
  out <- list('beta' = betas[, which.min(mse)], 
               'eps' = epsilon[, which.min(mse)])
  return(out)
}


#' @param b Simulation index.
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param sp Proportion of nonzero weights from Z to children.
#' @param snr Signal-to-noise ratio. 
#' @param xzr X-to-Z ratio.
#' @param form Functional form for structural equations.
#' @param l Norm to use, either 0, 1, or 2.

# Estimate causal directions by taking difference-in-differences between norms
our_fn <- function(b, n, d_z, rho, sp, snr, xzr, form, l) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho, sp, snr, xzr, form)
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Compute coefficients
  fit_fn <- function(h, f) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    zy <- cbind(z, y)
    betas <- list(
      f_fn(z[trn, ], y[trn], z[tst, ], y[tst], f, d_z)$beta,
      f_fn(zx[trn, ], y[trn], zx[tst, ], y[tst], f, d_z)$beta,
      f_fn(z[trn, ], x[trn], z[tst, ], x[tst], f, d_z)$beta,
      f_fn(zy[trn, ], x[trn], zy[tst, ], x[tst], f, d_z)$beta
    )
    v <- sapply(seq_len(4), function(f_idx) {
      switch(l, 'l0' = sum(betas[[f_idx]] != 0), 
                'l1' = sum(abs(betas[[f_idx]])),
                'l2' = sum(betas[[f_idx]]^2))
    })
    w0 <- v[1] - v[2]
    w1 <- v[3] - v[4]
    delta <- w0 - w1
    # Export
    out <- data.table('h' = h, 'f' = f, 'delta' = delta)
    return(out)
  }
  out <- foreach(aa = c('h0', 'h1'), .combine = rbind) %:%
    foreach(bb = c('lasso', 'step'), .combine = rbind) %do% 
    fit_fn(h = aa, f = bb)
  return(out)
}
res <- foreach(i = seq_len(500), .combine = rbind) %dopar% 
  our_fn(i, n = 500, d_z = 50, rho = 0.3, sp = 0.5, snr = 3, xzr = 1, 
         form = 'linear', l = 'l0')

# Simplified Entner test: evaluate R1 and R2 per Z using sparse regression
entner_fn <- function(b, n, d_z, rho, sp, snr, xzr, form, alpha) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho, sp, snr, xzr, form)
  z <- as.matrix(select(sim$dat, starts_with('z')))
  x <- sim$dat$x
  zx <- cbind(z, x)
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Fit models
  fit_fn <- function(h, f) {
    if (h == 'h0') y <- sim$dat$y0 else y <- sim$dat$y1
    # Three different designs
    f0 <- f_fn(z[trn, ], y[trn], z[tst, ], y[tst], f, d_z)
    f1 <- f_fn(zx[trn, ], y[trn], zx[tst, ], y[tst], f, d_z)
    f2 <- f_fn(z[trn, ], x[trn], z[tst, ], x[tst], f, d_z)
    # Evaluate rules
    #delta <- abs(f0$eps) - abs(f1$eps)
    #loco_p <- wilcox.test(delta, alt = 'greater')$p.value
    delta <- abs(f0$eps) - abs(f1$eps)
    loco_p <- binom.test(sum(sign(delta) > 0), length(delta), 
                         alt = 'greater')$p.value
    r1 <- f0$beta != 0 & f1$beta == 0
    r2 <- loco_p <= alpha | (f2$beta != 0 & f0$beta == 0)
    # Ground truth?
    R1 <- sim$wts$x != 0 & sim$wts$y == 0
    # Export results
    out <- data.table(
      'b' = b, 'h' = h, 'f' = f, 'R1' = R1, 'r1' = r1, 'r2' = r2
    )
    return(out)
  }
  out <- foreach(aa = c('h0', 'h1'), .combine = rbind) %:%
    foreach(bb = c('lasso', 'step'), .combine = rbind) %do% 
    fit_fn(aa, bb)
  return(out)
}
res <- foreach(i = seq_len(200), .combine = rbind) %dopar% 
  entner_fn(i, n = 500, d_z = 50, rho = 0.3, sp = 0.5, snr = 5, xzr = 1, 
            form = 'linear', alpha = 0.05)

n <- 500; d_z <- 50; rho <- 0.3; k <- 25; snr <- 5; xzr <- 1; 
form <- 'linear'; alpha <- 0.05


# Need to record which Z's are "instrumental", i.e. satisfy
# Z -> X -> Y and Z \ind Y | X. Define k := number of such Z's.
# Then we should expect k many Z's to satisfy R1 under h1 and 
# 0 to satisfy R1 under h0. We should also expect d_z many Z's
# to satisfy R2 under h0 and 0 to satisfy under h1.


# Score-based benchmark?


# Simulation settings
ns <- c(200, 500, 1000)
dim_zs <- c(20, 50, 100)
rhos <- c(0, 0.3, 0.5)
sps <- c(0.1, 0.5, 0.9)
snrs <- c(1, 3, 5)
xzrs <- c(0.5, 1, 2)
forms <- c('linear', 'nonlinear')
algs <- c('lasso', 'step', 'rf')



eps_fn <- function(des) {
  tmp <- as.data.frame(cbind(des, 'y' = y))
  f <- lm(y ~ ., data = tmp[trn, ])
  eps <- tmp$y[tst] - predict(f, tmp[tst, ])
  return(eps)
}
# Simplified Entner 2.0: loop through Z's (linear prototype)
entner2_fn <- function(b, n, d_z, rho, sp, snr, xzr, form, alpha) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho, sp, snr, xzr, form)
  dat <- as.data.frame(sim$dat)
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  zx <- cbind(z, x)
  # Split training and validation sets
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  # Fit models
  fit_fn <- function(h) {
    if (h == 'h0') y <- dat$y0 else y <- dat$y1
    # Residuals for full models
    eps_1i <- eps_fn(z)
    eps_1ii <- eps_fn(zx)
    # Residuals for null models
    p_i <- p_ii <- double(length = d_z)
    for (j in seq_len(d_z)) {
      eps_0i <- eps_fn(z[, -j])
      eps_0ii <- eps_fn(zx[, -j])
      delta_i <- abs(eps_0i) - abs(eps_1i)
      delta_ii <- abs(eps_0ii) - abs(eps_1ii)
      #p_i[j] <- wilcox.test(delta_i, alt = 'greater')$p.value
      #p_ii[j] <- wilcox.test(delta_ii, alt = 'less')$p.value
      #delta_i <- sign(abs(eps_0i) - abs(eps_1i))
      #delta_ii <- sign(abs(eps_0ii) - abs(eps_1ii))
      p_i[j] <- binom.test(sum(delta_i > 0), length(delta_i), 
                           alt = 'greater')$p.value
      p_ii[j] <- binom.test(sum(delta_ii > 0), length(delta_ii), 
                            alt = 'less')$p.value
    }
    r1 <- p_i <= alpha & p_ii <= alpha
    R1 <- sim$wts$x != 0 & sim$wts$y == 0
    # Export results
    out <- data.table('b' = b, 'h' = h, 'R1' = R1, 'r1' = r1)
    return(out)
  }
  out <- foreach(hh = c('h0', 'h1'), .combine = rbind) %do% fit_fn(hh)
  return(out)
}
res <- foreach(i = seq_len(200), .combine = rbind) %dopar% 
  entner2_fn(i, n = 500, d_z = 50, rho = 0.3, sp = 0.5, snr = 5, xzr = 1, 
             form = 'linear', alpha = 0.05)










