


library(gbm)

fit <- gbm.fit(x = z, y = x, distribution = 'gaussian', n.trees = 1000,
               interaction.depth = 2, nTrain = 0.8 * n)
n_tree <- gbm.perf(fit, plot.it = FALSE, method = 'test')
vimp <- summary(fit, plotit = FALSE, n.trees = n_tree, order = FALSE)

comp <- function(x, y, trn, tst) {
  d <- data.frame(x, y = y)
  fit <- gbm(y ~ ., data = d[trn, ], distribution = 'gaussian', n.trees = 1000,
             interaction.depth = 1, cv.folds = 5, n.cores = 5)
  n_tree <- gbm.perf(fit, plot.it = FALSE, method = 'cv')
  eps <- y[tst] - predict(fit, d[tst, ])
  mse <- mean(eps^2)
  
  f <- ranger(x = x[trn, ], y = y[trn])
  eps <- y[tst] - predict(f, x[tst, ])$predictions
  mse <- mean(eps^2)
}


# How many trees?
gbm_loop <- function(x, y) {
  d <- data.frame(y = y, x)
  fit <- gbm(y ~ ., data = d, distribution = 'gaussian', n.trees = 5000, 
             shrinkage = 0.1, interaction.depth = 1, train.fraction = 4/5, 
             n.cores = 1)
  n_tree <- gbm.perf(fit, plot.it = FALSE, method = 'test')
  return(n_tree)
}
outer_loop <- function(n, i) {
  # Simulate, extract
  sim <- sim_dat(n, d_z = 100, rho = 0.25, sp = 0.5, r2 = 2/3, 
                 lin_pr = 1/5, g = 'xy')
  dat <- sim$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- dat$x
  y <- dat$y
  zx <- cbind(z, x)
  zy <- cbind(z, y)
  d_z <- ncol(z)
  # Compute n_tree
  ntree_y0 <- gbm_loop(z, y)
  ntree_y1 <- gbm_loop(zx, y)
  ntree_x0 <- gbm_loop(z, x)
  ntree_x1 <- gbm_loop(zy, x)
  # Export
  out <- data.table(
    n, i, 
    model = c('y0', 'y1', 'x0', 'x1'), 
    ntree = c(ntree_y0, ntree_y1, ntree_x0, ntree_x1)
  )
  return(out)
}
df <- foreach(nn = c(2500, 5000, 1e4), .combine = rbind) %:%
  foreach(ii = 1:50, .combine = rbind) %dopar%
  outer_loop(nn, ii)

ggplot(df, aes(model, ntree, fill = model)) + 
  geom_boxplot() + 
  scale_fill_npg() + 
  theme_bw() + 
  facet_wrap(~ n)


# Here's what we learn: cap ntree at 2k for n = 1000, 2000
# Cap it at 3k for n = 4k. This is conservative, you could plausibly get 
# away with much less. All assumes shrinkage = 0.1. 










