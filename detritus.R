





# Four possible outputs: no edge, X -> Y, Y -> X, or ?
# Every method should be allowed to say each
# Strictly speaking, faithfulness coheres only with the L0 norm for our method
# Possible to establish some weaker assumptions that work with L1/L2?

# Define w as a difference in norms on Z weights
eval_fn <- function(b, n, d_z, rho, snr, xzr, l) {
  # Simulate data (with fixed weights)
  df <- sim_dat(n, d_z, rho, snr, xzr)
  # Function for computing redundancy of Z
  s_fn <- function(models, l) {
    # Compute norm
    v <- lapply(seq_len(4), function(i) {
      betas <- coef(models[[i]], s = 'lambda.min')[2:(d_z + 1)]
      if (l == 0) {
        sum(betas[[i]] != 0)
      } else if (l == 1) {
        sum(abs(betas[[i]]))
      } else if (l == 2) {
        sum(betas[[i]]^2)
      }
    })
    # Take difference of differences
    w0 <- v[[1]] - v[[2]]
    w1 <- v[[3]] - v[[4]]
    out <- w0 - w1
    # Export
    return(out)
  }
  ### H0: X \indep Y | Z ###
  f0 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$x[sample.int(n)])),
                  y = df$h0$y)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$x)),
                  y = df$h0$y)
  f2 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$y[sample.int(n)])), 
                  y = df$h0$x)
  f3 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$y)),
                  y = df$h0$x)
  s_null <- s_fn(list(f0, f1, f2, f3), l)
  ### H1: X -> Y
  f0 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$x[sample.int(n)])),
                  y = df$h1$y)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$x)),
                  y = df$h1$y)
  f2 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$y[sample.int(n)])), 
                  y = df$h1$x)
  f3 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$y)),
                  y = df$h1$x)
  s_alt <- s_fn(list(f0, f1, f2, f3), l)
  # Export
  out <- data.table('s_null' = s_null, 's_alt' = s_alt)
  return(out)
}
res <- foreach(i = seq_len(2000), .combine = rbind) %dopar% 
  eval_fn(i, n = 1000, d_z = 100, rho = 0.1, snr = 5, xzr = 1, l = 0)


# Variance about zero is a function of rho! (At least holding all else constant)
# Classifying who's a parent of Y and who is not. 
# The classification rule will be a threshold of abs(w)
# Just plot ROC curve -- is AUC higher for right direction than wrong direction

# Dividing by MSE or RMSE could help with scaling



# Benchmarking against:
# (i) Entner's test (partial correlation test), assume 1 IV
# (ii) Score-based: literally MSE with additive errors
# (iii) Score-based: BIC or similar

# First up, partial correlation test
trn <- sample(n, size = round(0.8 * n))
tst <- seq_len(n)[-trn]
rho_fn <- function(h) {
  f0 <- cv.glmnet(x = df$z[trn], y = df[[h]]$x[trn])
  f1 <- cv.glmnet(x = df$z[trn], y = df[[h]]$y[trn])
  x_tst <- df[[h]]$x[tst]
  y_tst <- df[[h]]$y[tst]
  eps_x <- x_tst - predict(f0, newx = df$z[tst], s = 'lambda.min')
  eps_y <- y_tst - predict(f1, newx = df$z[tst], s = 'lambda.min')
  pcor <- cor.test(eps_x, eps_y)
  if (pcor$p.value <= 0.05) {
    eps_xtoy <- residuals(lm(eps_y ~ eps_x))
    r2_xy <- cor(eps_x, eps_xtoy)^2
    eps_ytox <- residuals(lm(eps_x ~ eps_y))
    yx <- cor(eps_y, eps_ytox)^2
    if (r2_xy > r2_yx) {
      ruling <- 'y_to_x'
    } else {
      ruling <- 'x_to_y'
    }
  } else {
    ruling <- 'independent'
  }
  return(ruling)
}

out <- data.table(
  'h' = rep(c('h0', 'h1'), each = 4), 
  'x' = c('z', 'z', 'zy', 'zx'), 
  'y' = rep(c('x', 'y'), times = 2)
)
# Evaluate BIC and out-of-sample MSE
score_fn <- function(i) {
  # Select correct x and y
  if (nchar(out$x[i]) == 1) {
    x <- df$z
  } else {
    x <- cbind(df$z, df[[out$h[i]]][[substring(out$x[i], 2)]])
  }
  y <- df[[out$h[i]]][[out$y[i]]]
  # BIC
  f <- cv.glmnet(x = x, y = y)
  k <- sum(coef(f, s = 'lambda.min') != 0)
  y_hat <- predict(f, newx = x, s = 'lambda.min')
  eps_y <- y - y_hat
  rmse <- sqrt(mean(eps_y^2))
  bic <- -2 * sum(dnorm(eps_y, sd = rmse, log = TRUE)) + k * log(n)
  # MSE
  trn <- sample(n, size = round(0.8 * n))
  tst <- seq_len(n)[-trn]
  f <- cv.glmnet(x = x[trn], y = y[trn])
  y_hat <- predict(f, newx = x[tst], s = 'lambda.min')
  eps_y <- y[tst] - y_hat
  mse <- mean(eps_y^2)
  # Export
  return(data.table('bic' = bic, 'mse' = mse))
}
scores <- foreach(i = seq_len(nrow(out)), .combine = rbind) %dopar%
  score_fn(i)
out <- cbind(out, scores)

# If Y <- Z -> X, then p(x,y|z) factorizes as p(x|z)p(y|z)
# If Z -> X -> Y <- Z, then p(x,y|z) factorizes as p(x|z)p(y|x,z)
# If Z -> Y -> X <- Z, then p(x,y|z) factorizes as p(y|z)p(x|y,z)
# Calculate BIC for each and pick winner


# On Entner's rules:
# R1 says that, when predicting Y, if adding X to Z *activates any Z_j* that was 
# previously dormant, then infer X -> Y.
# R2 says that, if (a) X receives zero weight in E[Y|X,Z] OR
# (b) there exists some Z_j with zero weight in E[X|Z] but nonzero weight
# in E[Y|Z], then infer X \indep Y | Z.

# Question: Why just activating in R1, why not *deactivating* as well?













tmp <- data.table(
  'h' = rep(c('h0', 'h1'), each = 4), 
  'x' = rep(c('z', 'z', 'zy', 'zx'), times = 2), 
  'y' = rep(c('x', 'y'), times = 4)
)
score_fn <- function(i) {
  # Select correct x and y
  if (nchar(tmp$x[i]) == 1) {
    x <- df$z
  } else {
    x <- cbind(df$z, df[[tmp$h[i]]][[substring(tmp$x[i], 2)]])
  }
  y <- df[[tmp$h[i]]][[tmp$y[i]]]
  # BIC
  
  # Use forward selection with BIC, possibly initialized by lasso solution
  # If both X and Y get dropped from their respective regression, we 
  # declare conditional independence. Otherwise we orient edge according to 
  # sparsity pattern. Using forward selection for L0.
  f <- cv.glmnet(x = x, y = y)
  k <- sum(coef(f, s = 'lambda.min') != 0)
  y_hat <- predict(f, newx = x, s = 'lambda.min')
  eps_y <- y - y_hat
  rmse <- sqrt(mean(eps_y^2))
  bic <- -2 * sum(dnorm(eps_y, sd = rmse, log = TRUE)) + k * log(n)
  
  
}
scores <- foreach(i = seq_len(nrow(tmp)), .combine = rbind) %do%
  score_fn(i)
tmp <- cbind(tmp, scores)
bic_fn <- function(h) {
  tmp2 <- tmp[h == h]
  indep <- tmp2[, sum(bic)[1:2]]
  xtoy <- tmp2[, ]
}

bic_out <- foreach(h = c('h0', 'h1'), .combine = c) %do% bic_fn
indep_bic <- tmp[h == h, sum(bic[1:2])]
indep_mse <- tmp[h == h, sum(mse[1:2])]



# If Y <- Z -> X, then p(x,y|z) factorizes as p(x|z)p(y|z)
# If Z -> X -> Y <- Z, then p(x,y|z) factorizes as p(x|z)p(y|x,z)
# If Z -> Y -> X <- Z, then p(x,y|z) factorizes as p(y|z)p(x|y,z)

# Once model is chosen (either by L0 or L1), go back to Entner's rules
# For linear case, lasso implements L1 and greedy search implements L0
# For nonlinear case, deep lasso for L1 and RFE for L0



# Function for computing redundancy of Z
s_fn <- function(models, l) {
  # Extract Z weights
  betas <- lapply(seq_len(4), function(i) {
    coef(models[[i]], s = 'lambda.min')[2:(d_z + 1), ]
  })
  # Compute relevant norm
  v <- lapply(seq_len(4), function(i) {
    if (l == 0) {
      sum(betas[[i]] == 0)
    } else if (l == 1) {
      sum(abs(betas[[i]]))
    } else if (l == 2) {
      sum(betas[[i]]^2)
    }
  })
  # Take difference of differences
  w0 <- v[[1]] - v[[2]]
  w1 <- v[[3]] - v[[4]]
  out <- w0 - w1
  # Export
  return(out)
}

### H0: X \indep Y | Z ###
f0 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$x[sample.int(n)])),
                y = df$h0$y)
f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$x)),
                y = df$h0$y)
f2 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$y[sample.int(n)])), 
                y = df$h0$x)
f3 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$y)),
                y = df$h0$x)
s_null <- s_fn(list(f0, f1, f2, f3), l)
### H1: X -> Y
f0 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$x[sample.int(n)])),
                y = df$h1$y)
f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$x)),
                y = df$h1$y)
f2 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$y[sample.int(n)])), 
                y = df$h1$x)
f3 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$y)),
                y = df$h1$x)
s_alt <- s_fn(list(f0, f1, f2, f3), l)
# Export
out <- data.table('s_null' = s_null, 's_alt' = s_alt)
return(out)
}
res <- foreach(i = seq_len(2000), .combine = rbind) %dopar% 
  eval_fn(i, n = 1000, d_z = 100, rho = 0.1, snr = 5, xzr = 1, l = 2)


# Initialization: lasso E[X|Z] and E[Y|Z].
# For every edge Z -> X: either add or delete (depending on which is a change)
# and score results. Keep on adding/deleting until it doesn't help anymore.
# 







# Random forest version
eval_fn <- function(b, n, d_z, rho, snr, xzr) {
  # Simulate data (with fixed weights)
  df <- sim_dat(n, d_z, rho, snr, xzr)
  ### H0: X \indep Y | Z ###
  f0 <- randomForest(x = as.matrix(cbind(df$z, df$h0$x[sample.int(n)])),
                     y = df$h0$y, ntree = 100)
  f1 <- randomForest(x = as.matrix(cbind(df$z, df$h0$x)),
                     y = df$h0$y, ntree = 100)
  f2 <- randomForest(x = as.matrix(cbind(df$z, df$h0$y[sample.int(n)])), 
                     y = df$h0$x, ntree = 100)
  f3 <- randomForest(x = as.matrix(cbind(df$z, df$h0$y)),
                     y = df$h0$x, ntree = 100)
  v0 <- sum(importance(f0)[seq_len(d_z), ])
  v1 <- sum(importance(f1)[seq_len(d_z), ])
  v2 <- sum(importance(f2)[seq_len(d_z), ])
  v3 <- sum(importance(f3)[seq_len(d_z), ])
  w0 <- v0 - v1
  w1 <- v2 - v3
  s_null <- w0 - w1
  ### H1: X \dep Y | Z ###
  f0 <- randomForest(x = as.matrix(cbind(df$z, df$h1$x[sample.int(n)])),
                     y = df$h1$y, ntree = 100)
  f1 <- randomForest(x = as.matrix(cbind(df$z, df$h1$x)),
                     y = df$h1$y, ntree = 100)
  f2 <- randomForest(x = as.matrix(cbind(df$z, df$h1$y[sample.int(n)])), 
                     y = df$h1$x, ntree = 100)
  f3 <- randomForest(x = as.matrix(cbind(df$z, df$h1$y)),
                     y = df$h1$x, ntree = 100)
  v0 <- sum(importance(f0)[seq_len(d_z), ])
  v1 <- sum(importance(f1)[seq_len(d_z), ])
  v2 <- sum(importance(f2)[seq_len(d_z), ])
  v3 <- sum(importance(f3)[seq_len(d_z), ])
  w0 <- v0 - v1
  w1 <- v2 - v3
  s_alt <- w0 - w1
  # Export
  out <- data.table('s_null' = s_null, 's_alt' = s_alt)
  return(out)
}
res <- foreach(i = seq_len(2000), .combine = rbind) %dopar% 
  eval_fn(i, n = 1000, d_z = 100, rho = 0.1, snr = 5, xzr = 1)





# One idea: heteroskedastic errors as a failure mode for other methods
# 





# Ricardo's version: w := sum(sign(abs(coef0) - abs(coef1)))
w_fn <- function(b, n, d_z, rho, snr, xzr) {
  df <- sim_dat(n, d_z, rho, snr, xzr)
  f0 <- cv.glmnet(x = df$z, y = df$h0$y)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h0$x)),
                  y = df$h0$y)
  v0 <- abs(coef(f0, s = 'lambda.min')[2:(d_z + 1), ])
  v1 <- abs(coef(f1, s = 'lambda.min')[2:(d_z + 1), ])
  w0 <- sum(sign(v0 - v1))
  f0 <- cv.glmnet(x = df$z, y = df$h1$y)
  f1 <- cv.glmnet(x = as.matrix(cbind(df$z, df$h1$x)),
                  y = df$h1$y)
  v0 <- abs(coef(f0, s = 'lambda.min')[2:(d_z + 1), ])
  v1 <- abs(coef(f1, s = 'lambda.min')[2:(d_z + 1), ])
  w1 <- sum(sign(v0 - v1))
  out <- data.table('w0' = w0, 'w1' = w1)
  return(out)
}
res <- foreach(i = seq_len(2000), .combine = rbind) %dopar% 
  w_fn(i, n = 1000, d_z = 100, rho = 0, snr = 5, xzr = 1)




