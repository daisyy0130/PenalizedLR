library(dplyr)
library(SNPknock)

load("Input_dataset.RData")
#load("MCEM_Multiple-SNP.RData")
load("MCEM_Single-SNP.RData")

args <- commandArgs()
n <- as.numeric(args[6])
print(n)

m=4
K=10
set.seed(2020*n+K)
#sim_data=simUnmatched(n=200,Beta=log(rf(K,m,m)),p=0)     # continuous data
sim_data=simKnockoffGenotypes(n=200,Beta=log(rf(K,m,m)))  # SNP data
pp=profilelkhd(data=sim_data,mvals=c(1:20),N=1000)

save(pp,file=paste0("K=10_rep",n,".RData"))










