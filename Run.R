library(dplyr)
library(SNPknock)

load("Input_dataset.RData")
#load("MCEM_Single-SNP.RData")
load("MCEM_Multiple-SNP.RData")

# args <- commandArgs()
# n <- as.numeric(args[6])
# print(n)

m=4
K=10
rep=10
pp=matrix(data=NA,ncol=rep,nrow=20)
for (i in 1:rep) {
  #set.seed(2000*n+i+K)
  sim_data=simKnockoffGenotypes(n=200,Beta=log(rf(K,m,m)))
  pp[,i]=profilelkhd(data=sim_data,mvals=c(1:20),N=1000)
}

# save(pp,file=paste0("K=50_rep",n,".RData"))







