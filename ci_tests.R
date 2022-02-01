# Want to test speed and power of GCM VS. LOCO
# using RF and GBM subroutines

loco_test <- function(x, y, z, trn, val, tst, prms) {
  zx <- cbind(z, x)
  # Global params

  # Reduced model
  d0_trn <- lgb.Dataset(z[trn, ], label = y[trn])
  d0_val <- lgb.Dataset.create.valid(d0_trn, z[val, ], label = y[val])
  f0 <- lgb.train(params = prms, data = d0_trn, valids = list(val = d0_val), 
                  nrounds = 500, early_stopping_rounds = 5, verbose = 0)
  err0 <- abs(y[tst] - predict(f0, z[tst, ]))
  # Maximal model
  d1_trn <- lgb.Dataset(zx[trn, ], label = y[trn])
  d1_val <- lgb.Dataset.create.valid(d0_trn, zx[val, ], label = y[val])
  f1 <- lgb.train(params = prms, data = d1_trn, valids = list(val = d1_val), 
                  nrounds = 500, early_stopping_rounds = 5, verbose = 0)
  err1 <- abs(y[tst] - predict(f1, zx[tst, ]))
  # Wilcox test
  delta <- err0 - err1
  p_value <- wilcox.test(delta, alt = 'greater')$p.value
  return(p_value)
}


gcm_test <- function(x, y, z, trn, tst, f) {
  if (f == 'rf') {
    rf1 <- ranger(x = z[trn, ], y = x[trn], num.trees = 50, num.threads = 1)
    rf2 <- ranger(x = z[trn, ], y = y[trn], num.trees = 50, num.threads = 1)
    eps1 <- x[tst] - predict(rf1, z[tst, ], num.threads = 1)$predictions
    eps2 <- y[tst] - predict(rf2, z[tst, ], num.threads = 1)$predictions
  } else if (f == 'gbm') {
    # Need validation set for early stopping
    val <- sample(trn, round(0.1 * length(trn)))
    trn <- trn[-val]
    # Global params
    prms <- list(
      objective = 'regression', max_depth = 1, 
      bagging.fraction = 0.5, feature_fraction = 0.8, 
      num_threads = 1, force_col_wise = TRUE
    )
    # Model 1
    d1_trn <- lgb.Dataset(z[trn, ], label = y[trn])
    d1_val <- lgb.Dataset.create.valid(d1_trn, z[val, ], label = y[val])
    f1 <- lgb.train(params = prms, data = d1_trn, valids = list(val = d1_val), 
                    nrounds = 200, early_stopping_rounds = 5, verbose = 0)
    eps1 <- y[tst] - predict(f1, z[tst, ])
    # Model 2
    d2_trn <- lgb.Dataset(z[trn, ], label = x[trn])
    d2_val <- lgb.Dataset.create.valid(d2_trn, z[val, ], label = x[val])
    f2 <- lgb.train(params = prms, data = d2_trn, valids = list(val = d2_val), 
                    nrounds = 200, early_stopping_rounds = 5, verbose = 0)
    eps2 <- y[tst] - predict(f2, z[tst, ])
  }
  nn <- length(tst)
  R <- eps1 * eps2
  R.sq <- R^2
  meanR <- mean(R)
  z_score <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p_value)
}


# Test speed
library(microbenchmark)
loop_fn <- function(n, d_z = 100, method, f) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho = 0.25, sp = 0.5, r2 = 2/3, 
                 lin_pr = 1/5, g = 'xy')
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  z_b <- z[, sample(d_z, d_z/2)]
  x <- dat$x
  y <- dat$y
  # Test loop
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  if (method == 'loco') {
    out <- loco_test(x, y, z_b, trn, tst, f)
  } else if (method == 'gcm') {
    out <- gcm_test(x, y, z_b, trn, tst, f)
  }
  return(out)
}

microbenchmark(
  loco_rf = loop_fn(1000, 100, 'loco', 'rf'),
  loco_gb = loop_fn(1000, 100, 'loco', 'gbm'),
  gcm_rf = loop_fn(1000, 100, 'gcm', 'rf'),
  gcm_gb = loop_fn(1000, 100, 'gcm', 'gbm'),
  times = 75
)

# Result: Boosting much faster than random forest (>2x speedup)
# GCM very slightly faster than LOCO (~10 milliseconds on avg)
# Not worth it if one's much more powerful than the other! Let's see...

# LOCO clear winner, but try out different testing subroutines, ya?

loco_test <- function(x, y, z, trn, tst, b) {
  zx <- cbind(z, x)
  # Reduced model
  d0_trn <- lgb.Dataset(z[trn, ], label = y[trn])
  f0 <- lgb.train(params = prms, data = d0_trn, nrounds = b, verbose = 0)
  err0 <- abs(y[tst] - predict(f0, z[tst, ]))
  # Maximal model
  d1_trn <- lgb.Dataset(zx[trn, ], label = y[trn])
  f1 <- lgb.train(params = prms, data = d1_trn, nrounds = b, verbose = 0)
  err1 <- abs(y[tst] - predict(f1, zx[tst, ]))
  # Wilcox test
  delta <- err0 - err1
  out <- data.table(
         't' = t.test(delta, alt = 'greater')$p.value,
    'wilcox' = wilcox.test(delta, alt = 'greater')$p.value,
     'binom' = binom.test(x = sum(delta > 0), n = length(delta),
                         alt = 'greater')$p.value
  )
  return(out)
}

prms <- list(
  objective = 'regression', max_depth = 1, 
  bagging.fraction = 0.5, feature_fraction = 0.8, 
  num_threads = 1, force_col_wise = TRUE
)

pwr_fn <- function(n, b, d_z = 100, idx) {
  # Simulate data
  sim <- sim_dat(n, d_z, rho = 0.25, sp = 0.5, r2 = 2/3, 
                 lin_pr = 1/5, g = 'ci')
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  pa_x <- ifelse(sim$wts$beta != 0, 1, 0)
  pa_y <- ifelse(sim$wts$gamma != 0, 1, 0)
  shared_pa <- ifelse(pa_x + pa_y == 2, TRUE, FALSE)
  z_b <- z[, shared_pa]
  # Increasing X influence 
  w <- seq(0, 1, by = 0.05)
  y <- sapply(seq_along(w), function(k) dat$y + w[k] * x)
  # Test loop
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  out <- foreach(k = seq_along(w), .combine = rbind) %do%
    loco_test(x, y[, k], z_b, trn, tst, b) %>%
    mutate(n = n, w = w, b = b, idx = idx)
  return(out)
}
df <- foreach(ii = seq_len(100), .combine = rbind) %:%
  foreach(nn = c(2500, 5000, 1e4), .combine = rbind) %:%
  foreach(bb = c(25, 50, 100), .combine = rbind) %dopar%
  pwr_fn(n = nn, b = bb, idx = ii)


tmp <- df %>% 
  pivot_longer(t:binom, names_to = 'test', values_to = 'p_value') %>%
  as.data.table(.)
tmp[, hit_rate := sum(p_value <= 0.05)/.N, by = .(test, b, n, w)]
tmp <- unique(tmp[, .(test, n, b, w, hit_rate)])
tmp[, nn := paste('n =', n)]
tmp[, nn := factor(nn, levels = c('n = 2500', 'n = 5000', 'n = 10000'))]
tmp[, bb := paste('ntree =', b)]
tmp[, bb := factor(bb, levels = c('ntree = 100', 'ntree = 200', 
                                  'ntree = 500', 'ntree = 1000'))]
ggplot(tmp, aes(w, hit_rate, color = test, group = test)) + 
  geom_path() + 
  geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red') + 
  scale_color_d3() +
  theme_bw() + 
  facet_grid(bb ~ nn)

# Weird NA's on the t-test?
# t-test is the slight winner but Wilcox is basically equivalent
# binomial a distant third
# Interestingly: nothing seems to control type I error at this sample size
# Do we need larger n?

# Super unexpected: fewer trees is associated with better type I error control
# and has effectively ZERO impact on statistical power. So the fewer trees
# the better, basically? So weird!


ggplot(df, aes(w, p_value, color = method, group = method)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  geom_smooth(size = 0.5, method = 'loess') + 
  theme_bw()

df[, hit := ifelse(p_value <= 0.05, 1, 0)]
df[, hit_rate := sum(hit == 1)/.N, by = .(method, w)]
tmp <- unique(df[, .(method, w, hit_rate)])
ggplot(tmp, aes(w, hit_rate, color = method, group = method)) + 
  geom_point(size = 0.25, alpha = 0.5) + 
  geom_path(size = 0.5) + 
  geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red') + 
  theme_bw()





































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
      tmp_f <- randomForest(tmp_x, trn_y, ntree = 200)
      y_hat <- predict(tmp_f, newdata = tst_x)
      return(y_hat)
    })
  }
  eps_mat <- y_hat - tst_y
  mse <- colMeans(eps_mat^2)
  eps <- eps_mat[, which.min(mse)]
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
loco_test <- function(x, y, z, f) {
  eps_f0 <- eps_fn(trn_x = x[trn, ], trn_y = y[trn],
                   tst_x = x[tst, ], tst_y = y[tst], f)
  eps_f1 <- eps_fn(trn_x = cbind(x[trn, ], z[trn]), trn_y = y[trn],
                   tst_x = cbind(x[tst, ], z[tst]), tst_y = y[tst], f)
  delta <- abs(eps_f0) - abs(eps_f1)
  p.value <- wilcox.test(delta, alt = 'greater')$p.value
  return(p.value)
}



gcm_test <- function(x, y, z, trn, tst) {
  rf1 <- ranger(x = z[trn, ], y = x[trn], num.trees = 100, num.threads = 1)
  rf2 <- ranger(x = z[trn, ], y = y[trn], num.trees = 100, num.threads = 1)
  eps1 <- x[tst] - predict(rf1, z[tst, ], num.threads = 1)$predictions
  eps2 <- y[tst] - predict(rf2, z[tst, ], num.threads = 1)$predictions
  nn <- length(tst)
  R <- eps1 * eps2
  R.sq <- R^2
  meanR <- mean(R)
  z_score <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p.value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p.value)
}




# With GCM
f_x <- lm(x ~ z)
f_y <- lm(y ~ z)
for (j in seq_len(d_z)) {
  ### Rule 1 ###
  # Rule 1(i)
  f0_w <- lm(z[, j] ~ z[, -j])
  f0_y <- lm(y ~ z[, -j])
  p1.i <- gcm_test(residuals(f0_w), residuals(f0_y))
  rule1.i <- ifelse(p1.i <= alpha, TRUE, FALSE)
  # Rule 1(ii) ...Both must put nonzero weight on X tho?
  f1_w <- lm(z[, j] ~ zx[, -j])
  f1_y <- lm(y ~ zx[, -j])
  p1.ii <- gcm_test(residuals(f1_w), residuals(f1_y))
  rule1.ii <- ifelse(p1.ii <= alpha, FALSE, TRUE)
  # Therefore
  rule1 <- ifelse(rule1.i & rule1.ii, TRUE, FALSE)
  ### Rule 2 ###
  # Rule 2(i)
  p2.i <- gcm_test(residuals(f_x), residuals(f_y))
  rule2.i <- ifelse(p2.i <= alpha, FALSE, TRUE)
  # Rule 2(ii)
  f0_x <- lm(x ~ z[, -j])
  p2.ii <- gcm_test(residuals(f0_w), residuals(f0_x))
  rule2.ii <- ifelse(p2.ii <= alpha, TRUE, FALSE)
  # Rule 2(iii)
  rule2.iii <- !rule1.i
  # Therefore
  rule2 <- ifelse(rule2.i | (rule2.ii & rule2.iii), TRUE, FALSE)
}



# With LOCO
ci_outer <- function(h, f) {
  if (h == 'h0') y <- dat$y0 else y <- dat$y1
  p2.i <- loco_test(x = z, y = y, z = x, f)
  ci_inner <- function(j, f) {
    ### Rule 1 ###
    # Rule 1(i)
    p1.i <- loco_test(x = z[, -j], y = y, z = z[, j], f)
    # Rule 1(ii) ...Both must put nonzero weight on X tho?
    p1.ii <- loco_test(x = zx[, -j], y = y, z = zx[, j], f)
    ### Rule 2 ###
    # Rule 2(ii)
    p2.ii <- loco_test(x = z[, -j], y = x, z = z[, j], f)
    # Export
    out <- data.table(
      'j' = j, 'h' = h, 'f' = f,
      'idx' = c('1.i', '1.ii', '2.ii'),
      'p.value' = c(p1.i, p1.ii, p2.ii), 
    )
    return(out)
  }
  out <- foreach(a = seq_len(d_z), .combine = rbind) %dopar% ci_inner(a, f = f)
  out <- rbind(
    data.table(
      'j' = 0, 'h' = h, 'f' = f, 'idx' = '2.i', 'p.value' = p2.i
    ), out
  )
  return(out)
}
res <- foreach(a = c('h0', 'h1'), .combine = rbind) %do% ci_outer(a, f = 'lasso')
# Adjust p-values
res[, q.value := p.adjust(p.value, method = 'fdr'), by = h]
# Reprint 2.i for each j
res_2.i <- foreach(a = seq_len(d_z), .combine = rbind) %do%
  cbind(data.table('j' = a), res[j == 0, .(h, f, idx, p.value, q.value)])
res <- rbind(res[j != 0], res_2.i)
# Compute R1 and R2 for each feature
res[, r1 := ifelse(.SD[idx == '1.i', q.value <= alpha] & 
                     .SD[idx == '1.ii', q.value > alpha], TRUE, FALSE), by = .(j, h)]
res[, r2 := ifelse(.SD[idx == '2.i', q.value > alpha] | 
                     (.SD[idx == '2.ii', q.value <= alpha] & 
                        .SD[idx == '1.i', q.value > alpha]), TRUE, FALSE), by = .(j, h)]

