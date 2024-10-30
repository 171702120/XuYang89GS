dir<-"E:/rice/data"
setwd(dir)
phe<-read.csv(file="phe.csv",row.names = 1)#Phenotypic data of parents
gen<-read.csv(file="gen.csv",row.names = 1)#Genotype of parents
hphe<-read.csv(file="hphe.csv",row.names = 1)#Phenotypic data of hybrids
met<-read.csv(file="met.csv",row.names = 1)#Metabolomic data of parents
met<-scale(met)

#Identify metabolites significantly associated with agronomic trait(metabolic markers)
library(lassopv)
y<-phe[,1]
x<-met
pv=lassopv(x,y)
p<-which(pv<0.05)

#
met1<-met[,p]  

#Inferring genotype of hybrid 
pn_g<-rownames(gen)
pf<-as.character(hphe[,1])
pm<-as.character(hphe[,2])
gen<-as.matrix(gen)
pf_gen<-gen[match(pf,pn_g),]
pm_gen<-gen[match(pm,pn_g),]
add<-(pf_gen+pm_gen)/2

#
pn_m<-rownames(met1)
pf<-as.character(hphe[,1])
pm<-as.character(hphe[,2])
met1<-as.matrix(met1)
pf_met<-met1[match(pf,pn_m),]
pm_met<-met1[match(pm,pn_m),]

#Additive coding matrix of metabolites
add_m<-(pf_met+pm_met)/2

#Dominance coding matrix of metabolites
dom_m<-abs(pf_met-pm_met)/2

#Calculate Kinship Matrix
m <- ncol(add)
n <- nrow(add)
ka<- matrix(0, n, n)
for (k in 1:m) {
  z <- add[, k, drop = F]
  ka<- ka+ z %*% t(z)
}
ka<- ka/m

#
m <- ncol(add_m)
n <- nrow(add_m)
kma<- matrix(0, n, n)
for (k in 1:m) {
  z <- add_m[, k, drop = F]
  kma<- kma+ z %*% t(z)
}
kma<- kma/m

#
m <- ncol(dom_m)
n <- nrow(dom_m)
kmd<- matrix(0, n, n)
for (k in 1:m) {
  z <- dom_m[, k, drop = F]
  kmd<- kmd+ z %*% t(z)
}
kmd<- kmd/m

#GBLUP for MM_GP
ka21 <- NULL
phe_name <- NULL

predparent_gen <- t(gen)
pred_name <- colnames(predparent_gen)
predparent_gen <- as.matrix(predparent_gen)
for (i in 1:(ncol(predparent_gen) - 1)) {
  ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
  ka2 <- tcrossprod(ha1, add)/ncol(add)
  ka21 <- rbind(ka21, ka2)
  test_name <- paste(pred_name[i], pred_name[-(1:i)],sep = "/")
  phe_name <- c(phe_name, test_name)
}

kma21<-NULL
kmd21<-NULL
predparent_m <- t(met1)
pred_name_m<- colnames(predparent_m)
predparent_m <- as.matrix(predparent_m)
for (i in 1:(ncol(predparent_m) - 1)) {
  hma1<- t((predparent_m[, i] + predparent_m[,-(1:i)])/2)
  kma2 <- tcrossprod(hma1,add_m)/ncol(add_m)
  kma21 <- rbind(kma21, kma2)
  
  hmd1 <- abs(t((predparent_m[, i] - predparent_m[, -(1:i)])/2))
  kmd2 <- tcrossprod(hmd1,dom_m)/ncol(dom_m)
  kmd21 <- rbind(kmd21, kmd2)
}

#
mixed<- function(fix = NULL, y, kk) {
  n <- length(y)
  y <- as.matrix(y)
  if (is.null(fix)) {
    fix <- matrix(1, n, 1)
  } else {
    fix <- as.matrix(fix)
  }
  ##Negative nloglikelihood Function
  g <- length(kk)
  nloglik_REML <- function(pm) {
    v_phi <- 0
    for (p in 1:g) {
      v_phi <- v_phi + kk[[p]] * pm[p]
    }
    v_sigma <- diag(n) * pm[g + 1]
    v <- v_phi + v_sigma + diag(1e-09, n)
    v_i <- solve(v, tol = -50)
    beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
    nloglik <- 0.5 * (unlist(determinant(v))[[1]] + unlist(determinant(t(fix) %*%
                                                                         v_i %*% fix))[[1]] + t(y - fix %*% beta) %*% v_i %*% (y - fix %*% beta))
    return(nloglik)
  }
  parm0 <- rep(1, g + 1)
  parm <- optim(par = parm0, fn = nloglik_REML, method = "L-BFGS-B", hessian = FALSE,
                lower = 0)
  v_phi <- 0
  for (p in 1:g) {
    v_phi <- v_phi + kk[[p]] * parm$par[p]
  }
  v_sigma <- diag(n) * parm$par[g + 1]
  v <- v_phi + v_sigma
  v_i <- solve(v, tol = -50)
  beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
  res <- list(v_i = v_i, var = parm$par[-(g + 1)], ve = parm$par[g + 1], beta = beta)
  return(res)
}

#Predict the Performance of Hybrids
y <- hphe[,3]
parm<-mixed(y = y, kk = list(ka,kma,kmd))
v_i<-parm$v_i
beta<-parm$beta
va<-parm$var[1]
vpa<-parm$var[2]
vpd<-parm$var[3]
ve<-parm$ve

fix <- matrix(1, n, 1)
n1 <- nrow(ka21)
fixnew <- matrix(1, n1, 1)
G21 <- ka21*va+kma21*vpa+kmd21*vpd

pred_phe <- fixnew %*% beta + G21 %*% v_i %*% (y - fix %*% beta)
row.names(pred_phe) <- phe_name

#"all"
Results_select <- as.data.frame(pred_phe)
colnames(Results_select) <- paste("all_", nrow(pred_phe), sep = "")
write.csv(Results_select,"E:/rice/predict/GBLUP/yield_all.csv")

#select == "top"
Results_select <- as.data.frame(sort(pred_phe[, 1], decreasing = T)[c(1:150)])
names(Results_select) <- paste("top_", 150, sep = "")
write.csv(Results_select,"E:/rice/predict/GBLUP/yield_top.csv")

# select == "bottom"
Results_select <- as.data.frame(sort(pred_phe[, 1], decreasing = F)[c(1:150)])
colnames(Results_select) <- paste("bottom_", 150,sep="")
write.csv(Results_select,"E:/rice/predict/GBLUP/yield_bottom.csv")


#XGBoost for MM_GP
ha1_xgb<-NULL
phe_name<- NULL
predparent_gen <- t(gen)
pred_name <- colnames(predparent_gen)
predparent_gen <- as.matrix(predparent_gen)
for (i in 1:(ncol(predparent_gen) - 1)) {
  ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
  ha1_xgb<- rbind(ha1_xgb, ha1)
  test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
  phe_name<- c(phe_name, test_name)
}

#
hma1_xgb<-NULL
hmd1_xgb<-NULL

predparent_m <- t(met1)
pred_name_m<- colnames(predparent_m)
predparent_m <- as.matrix(predparent_m)
for (i in 1:(ncol(predparent_m) - 1)) {
  hma1<- t((predparent_m[, i] + predparent_m[, -(1:i)])/2)
  hma1_xgb <- rbind(hma1_xgb, hma1)
  
  hmd1 <- abs(t((predparent_m[, i] - predparent_m[, -(1:i)])/2))
  hmd1_xgb <- rbind(hmd1_xgb, hmd1)
}

##Predict the Performance of Hybrids
x1 <- cbind(add,add_m,dom_m)
y1 <- hphe[,3]

library(xgboost)
xg <- xgboost(data = x1, label = y1,
               nrounds=206                         
              ,colsample_bynode=0.9
              ,colsample_bytree=0.4
              ,eta=0.03
              ,gamma=9
              ,reg_lambda=1
              ,max_depth=3
              #,min_child_weight=24
              ,subsample=0.7
              ,set.seed(2)
              ,verbose = FALSE
              ,nthread=5
)

x2 <- cbind(ha1_xgb,hma1_xgb,hmd1_xgb)
yhat <- predict(xg, x2)
yhat<-as.matrix(yhat)
row.names(yhat)<-phe_name

#"all"
Results_select <- as.data.frame(yhat)
colnames(Results_select) <- paste("all_", nrow(yhat), sep = "")
write.csv(Results_select,"E:/rice/predict/XGBoost/yield_all.csv")

#select == "top"
Results_select <- as.data.frame(sort(yhat[, 1],decreasing = T)[c(1:150)])
names(Results_select) <- paste("top_", 150, sep = "")
write.csv(Results_select,"E:/rice/predict/XGBoost/yield_top.csv")

# select == "bottom"
Results_select <- as.data.frame(sort(yhat[, 1], decreasing = F)[c(1:150)])
colnames(Results_select) <- paste("bottom_", 150,sep="")
write.csv(Results_select,"E:/rice/predict/XGBoost/yield_bottom.csv")