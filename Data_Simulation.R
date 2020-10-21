simMAF=function(n,Beta,p) {
  # Input:
  # - n is total sample size
  # - beta is value of parameter of interest
  # - p is the pre-determined MAF
  # Output: sim_data
  
  ncov=length(Beta)
  data=NULL
  for(i in 1:ncov) {
    data=cbind(data,rbinom(10000,2,p))
  }
  data=as.data.frame(data)
  data$case=NULL
  
  for (i in 1:dim(data)[1]) {
    s=sum(data[i,]*Beta)
    prob=exp(s)/(1+exp(s))
    data$case[i]=rbinom(1,1,prob)
  }
  
  # Sample case/control status
  # case-control ratio: 1:1
  no_case=n/2
  no_con=n/2
  data_con=data %>% filter(case==0)
  data_case=data %>% filter(case==1)
  sim_data=rbind(data_con[sample(nrow(data_con),no_con,replace=F),],
                 data_case[sample(nrow(data_case),no_case,replace=F),])
  colnames(sim_data)=c(paste0("X",1:K),"case");rownames(sim_data)=NULL
  return(sim_data)
}


simKnockoffGenotypes=function(n,Beta=log(rf(K,m,m))) {
  # Input:
  # - n is the sample size
  # - Beta is the value of log-OR simulated from log-F(m,m) distribution
  # Output: sim_data
  
  Xinp_file=writeXtoInp(X)
  fp_path="/scratch/yya188/fastPHASE"
  #fp_path="/Users/daisyyu/Desktop/Simulation/fastPHASE"
  fp_outPath=runFastPhase(fp_path, Xinp_file, K=12, numit=25)
  r_file=paste(fp_outPath,"_rhat.txt",sep="")
  alpha_file=paste(fp_outPath,"_alphahat.txt",sep="")
  theta_file=paste(fp_outPath,"_thetahat.txt",sep="")
  char_file=paste(fp_outPath,"_origchars",sep="")
  hmm=loadHMM(r_file,alpha_file,theta_file,char_file)
  
  # Simulate knockoff genotypes
  seed=sample(10000:50000,20,replace=F)
  data=as.data.frame(NULL)
  for (i in 1:20) {
    Xk=knockoffGenotypes(X,hmm$r,hmm$alpha,hmm$theta,seed=seed[i])
    data=rbind(data,Xk)
  }
  
  # Calculate disease probability for each individual
  # choose a subset of K SNPs to be casually associated with the disease
  cov_ind=sample(1:dim(data)[2],K,replace=F)
  data=data[,cov_ind]
  data$case=NULL
  
  for (i in 1:dim(data)[1]) {
    s=sum(data[i,]*Beta)
    prob=exp(s)/(1+exp(s))
    data$case[i]=rbinom(1,1,prob)
  }
  
  # Sample case/control status
  # case-control ratio: 1:1
  no_case=n/2
  no_con=n/2
  data_con=data %>% filter(case==0)
  data_case=data %>% filter(case==1)
  sim_data=rbind(data_con[sample(nrow(data_con),no_con,replace=F),],
                 data_case[sample(nrow(data_case),no_case,replace=F),])
  colnames(sim_data)=c(paste0("X",1:K),"case");rownames(sim_data)=NULL
  return (sim_data)
}


simUnmatched=function(n,Beta,p,scale=FALSE) {
  # Input:
  # - n is total sample size
  # - beta is value of parameter of interest
  # - p is number of nuisance covariates
  # Output: sim_data

  ncase=n/2; ncon=n/2  # assuming 1:1 con:case ratio
  Beta=c(Beta,rep(1,p))
  ncov=p+length(Beta)
  conX=caseX = NULL
  for(i in 1:ncov) {
    conX=cbind(conX,rnorm(ncon,mean=0,sd=1))
    caseX=cbind(caseX,rnorm(ncase,mean=Beta[i],sd=1))
  }
  X=rbind(caseX,conX)
  if(scale) X = round(scale(X))
  colnames(X)=paste0("X",1:ncov);rownames(X) = NULL
  case=c(rep(1,ncase),rep(0,ncon))
  return (data.frame(X,case))
}
