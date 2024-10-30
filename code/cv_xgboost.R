library(doParallel)
library(foreach)
library(parallel)
library(xgboost)
cv_xgboost<- function(x, y, nfold=10, seed=1, CPU=1) {
  
  set.seed(seed)
  n<-length(y)
  foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
  
  yp <- NULL
  yo <- NULL
  cl.cores <- detectCores()
  if (cl.cores <= 2 || CPU <= 1) {
    cl.cores <- 1
  }
  else if (cl.cores > 2) {
    if (cl.cores > 10) {
      cl.cores <- 10
    }
    else {
      cl.cores <- detectCores() - 1
    }
  }
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  
  res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", 
                 .packages = c("xgboost")) %dopar% {
                   id1 <- which(foldid != k)
                   id2 <- which(foldid == k)
                   x1 <- x[id1, ]
                   x2 <- x[id2, ]
                   y1 <- y[id1]
                   y2 <- y[id2]
                   xg <- xgboost(data = x1, label = y1,
                                 nrounds=206                         
                                 
                                 ,colsample_bynode=0.9
                                 ,colsample_bytree=0.4
                                 ,eta=0.03
                                 ,gamma=9
                                 ,reg_lambda=1
                                 
                                 
                                 ,max_depth=3
                                 # ,min_child_weight=24
                                 ,subsample=0.7
                                 ,set.seed(2),verbose = FALSE
                                 ,nthread=5
                   )
                   
                   yhat <- predict(xg, x2)
                   yp <- yhat
                   yo <- y2
                   data.frame(yp, yo)
                 }
  stopCluster(cl)
  r2 <- cor(res$yo, res$yp)^2
  
  print(r2)
}