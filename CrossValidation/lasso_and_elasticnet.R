
##=======================================================
##    cross validation for lasso and elastic-net model
##=======================================================

cvlasso<-function (x, y, lambda, blocks, sig) {
  N <- length(y)
  blocklist <- split(1:N, rep(1:blocks, length = N))
  blocklist1 <- split(which(y == 1), rep(1:blocks, length = length(which(y == 1))))
  blocklist2 <- split(which(y == 0), rep(1:blocks, length = length(which(y == 0))))
  for (i in 1:blocks) {
    blocklist[[i]] = c(blocklist1[[i]], blocklist2[[i]])
  }
  K <- length(blocklist)
  cutoff = seq(0.005, 0.5, 0.005)
  numcut = length(cutoff)
  gmeanmax <- array(0, c(K, length(lambda), numcut))
  avergmax <- array(0, c(K, length(lambda), numcut))
  overall <- array(0, c(K, length(lambda), numcut))
  for (i in seq(K)) {
    omit <- blocklist[[i]]
    x_train = x[-omit, ]
    y_train = y[-omit]
    x_test = x[omit, ]
    y_test = y[omit]
    fit <- glmnet(scale(x_train), y_train, lambda = lambda, family = "binomial")  
    # if this is for elastic-net model, then  fit <- glmnet(scale(x_train), y_train, lambda = lambda, family = "binomial",alpha=0.5)
    
    x_test = cbind(rep(1, length(y_test)), scale(x_test))
    coef = coef(fit)   ## coefficient matrix after glnmet function
    testpred = inv_logit_func(x_test %*% coef)
    for (j in 1:length(lambda)) {
      for (l in 1:numcut) {
        result = testpred[, j] >= cutoff[l]
        table = matrix(0, 2, 2)
        table[1, 1] = sum(result & (y_test == 1))
        table[1, 2] = sum(result & (y_test == 0))
        table[2, 1] = sum(!result & (y_test == 1))
        table[2, 2] = sum(!result & (y_test == 0))
        tpr = table[1, 1]/(table[1, 1] + table[2, 1])
        tnr = table[2, 2]/(table[1, 2] + table[2, 2])
        gmeanmax[i, j, l] = sqrt(tpr * tnr)
        avergmax[i, j, l] = (tpr + tnr)/2
        overall[i, j, l] = (table[1, 1] + table[2, 2])/sum(table)
      }
    }
  }
  gmeanmse <- gmeanmax[1, , ]
  avergmse <- avergmax[1, , ]
  overallmse <- overall[1, , ]
  for (i in 2:K) {
    gmeanmse <- gmeanmse + gmeanmax[i, , ]
    avergmse <- avergmse + avergmax[i, , ]
    overallmse <- overallmse + overall[i, , ]
  }
  result <- list(gmeanmse = gmeanmse, avergmse = avergmse, 
                 overallmse = overallmse, lambda = lambda)
}
