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

m <- ncol(add_m)
n <- nrow(add_m)
kma<- matrix(0, n, n)
for (k in 1:m) {
  z <- add_m[, k, drop = F]
  kma<- kma+ z %*% t(z)
}
kma<- kma/m

m <- ncol(dom_m)
n <- nrow(dom_m)
kmd<- matrix(0, n, n)
for (k in 1:m) {
  z <- dom_m[, k, drop = F]
  kmd<- kmd+ z %*% t(z)
}
kmd<- kmd/m

#GBLUP for MM_GP
dir<-"E:/rice/code"
setwd(dir)
source(file="cv_gblup.R")
cv_gblup(y=hphe[,3],kk=list(ka,kma,kmd))

#XGBoost for MM_GP
dir<-"E:/rice/code"
setwd(dir)
source(file="cv_xgboost.R")
x<-cbind(add,add_m,dom_m)
cv_xgboost(x,y=hphe[,3])
