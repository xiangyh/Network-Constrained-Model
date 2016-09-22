
##=====================================================
##    cross validation for network-constrained model
##=====================================================

cvlasso<-function (x, y, lambda, blocks, lambda2, sig) {
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
    x_g = as.matrix(x_train)
    test_mat = t(x_g) %*% x_g
    if (det(test_mat) == 0) {
      x_g[row(x_g) == col(x_g)] = x_g[row(x_g) == col(x_g)] + 
        0.01
      x_train = x_g
    }
    mu.x <- apply(x_g, 2, mean)  ## the mean of each column of X
    x_g <- sweep(x_g, 2, mu.x)  ## centralize x_g
    decomp <- qr(x_g)  ## QR decomposition
    if (decomp$rank < dim(x_g)[2]) {
      x_train[1:dim(x_g)[2], ] = diag(dim(x_g)[2])
      x_train[-(1:dim(x_g)[2]), ] = 0
    }
    if (lambda2 == 0) {
      traindata = cbind(x_train, y_train)
      traindata = data.frame(traindata)
    }
    if (lambda2 != 0) {
      trainExtend = argumentMatrxIndependent(x_train, lambda2, sig = 0)
      x_train = 1/sqrt(1 + lambda2) * rbind(scale(x_train), trainExtend)
      y_train = c(y_train, rep(0, dim(trainExtend)[1]))
      traindata = cbind(x_train, y_train)
      traindata = data.frame(traindata)
    }
    fit <- glmnet(x_train, y_train, lambda = lambda, 
                  family = "binomial")
    x_test = cbind(rep(1, length(y_test)), scale(x_test))
    coef = 1/sqrt(1 + lambda2) * coef(fit)   ## coefficient matrix after glnmet function
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
