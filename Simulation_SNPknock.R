devtools::install_github("msesia/snpknock")
library(SNPknock)
library(dplyr)

## Fit the Hidden Markov Model on Genotype data with fastPHASE
# Load input genotype data on 179 CNs for 69 SNPs in NEDD9
load("/Users/daisyyu/Desktop/Input_dataset.RData") 
Xinp_file=writeXtoInp(X)
fp_path="/Users/daisyyu/Desktop/Simulation/fastPHASE" 
fp_outPath=runFastPhase(fp_path, Xinp_file, K=12, numit=30)
r_file=paste(fp_outPath,"_rhat.txt",sep="")
alpha_file=paste(fp_outPath,"_alphahat.txt",sep="")
theta_file=paste(fp_outPath,"_thetahat.txt",sep="")
char_file=paste(fp_outPath,"_origchars",sep="")
hmm=loadHMM(r_file,alpha_file,theta_file,char_file)


simKnockoffGenotypes=function(n,Beta=log(rf(K,m,m)),hmm) {
  # Input: 
  # - n is the sample size
  # - Beta is the value of log-OR simulated from log-F(m,m) distribution
  # - hmm is the hidden Markov model contructed by fastPHASE
  # Output: sim_data
  
  # Simulate knockoff genotypes 
  seed=sample(10000:50000,60,replace=F)
  data=as.data.frame(NULL)
  for (i in 1:60) {
    Xk=knockoffGenotypes(X,hmm$r,hmm$alpha,hmm$theta,seed=seed[i])
    data=rbind(data,Xk)
  }
  
  # Calculate disease probability for each individual
  # choose a subset of SNPs to be casually associated with the disease
  cov_ind=sample(1:69,K,replace=F) 
  data$case=NULL
  for (i in 1:dim(data)[1]) {
    s=sum(data[i,cov_ind]*Beta)
    prob=exp(-5+s)/(1+exp(-5+s))
    if (prob >= 0.5) {data$case[i]=1}
    else {data$case[i]=0}
  }
  
  # Sample case/control status
  # case-control ratio: 1:4
  no_case=n/5*1
  no_con=n/5*4
  colnames(data)[1:69]=seq(1:69)
  data_con=data %>% filter(case==0) %>% select(cov_ind,case)
  data_case=data %>% filter(case==1) %>% select(cov_ind,case)
  sim_data=rbind(data_con[sample(nrow(data_con),no_con,replace=F),],
                 data_case[sample(nrow(data_case),no_case,replace=F),]) 
  colnames(sim_data)=c(paste0("x",1:K),"case");rownames(sim_data)=NULL
  return (sim_data)
}

m=2
K=30
n=200
sim_data=simKnockoffGenotypes(n,Beta=log(rf(K,m,m)),hmm)
