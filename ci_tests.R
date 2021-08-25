
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

