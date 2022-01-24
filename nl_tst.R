
# Compare RF-RFE with XGB-RFE
rf_loop <- function(x, y, trn, tst) {
  fit <- randomForest(x[trn, ], y[trn], ntree = 200)
  vimp <- data.frame('feature' = colnames(x), 
                     'imp' = as.numeric(importance(fit))) %>%
    arrange(desc(imp))
  beta <- double(length = d_z)
  names(beta) <- paste0('z', seq_len(d_z))
  s <- subsets(m = 10, max_d = d_z, min_d = 5, decay = 2)
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
  betas <- cbind(betas, as.numeric(importance(fit))[seq_len(d_z)])
  out <- list(y_hat, betas)
  return(out)
}

xgb_loop <- function(x, y, trn, tst) {
  
}

  