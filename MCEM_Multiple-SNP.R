simKnockoffGenotypes=function(n,Beta=log(rf(K,m,m))) {
  # Input: 
  # - n is the sample size
  # - Beta is the value of log-OR simulated from log-F(m,m) distribution
  # Output: sim_data
  
  # Fit the Hidden Markov Model on Genotype data with fastPHASE
  # load input genotype data on 179 CNs for 69 SNPs in NEDD9
  Xinp_file=writeXtoInp(X)
  fp_path="/scratch/yya188/fastPHASE" 
  #fp_path="/Users/daisyyu/Desktop/Simulation/fastPHASE"
  fp_outPath=runFastPhase(fp_path, Xinp_file, K=12, numit=30)
  r_file=paste(fp_outPath,"_rhat.txt",sep="")
  alpha_file=paste(fp_outPath,"_alphahat.txt",sep="")
  theta_file=paste(fp_outPath,"_thetahat.txt",sep="")
  char_file=paste(fp_outPath,"_origchars",sep="")
  hmm=loadHMM(r_file,alpha_file,theta_file,char_file)
  
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
  data=data[,cov_ind]
  data$case=NULL
  
  for (i in 1:dim(data)[1]) {
    s=sum(data[i,]*Beta)
    prob=exp(s)/(1+exp(s))
    data$case[i]=rbinom(1,1,prob)
  }
  
  # Sample case/control status
  # case-control ratio: 1:4
  no_case=n/5*1
  no_con=n/5*4
  data_con=data %>% filter(case==0) 
  data_case=data %>% filter(case==1) 
  sim_data=rbind(data_con[sample(nrow(data_con),no_con,replace=F),],
                 data_case[sample(nrow(data_case),no_case,replace=F),]) 
  colnames(sim_data)=c(paste0("x",1:K),"case");rownames(sim_data)=NULL
  return (sim_data)
}


MCEM=function(m,data,N) {
  # Input: 
  # - m is the value of m
  # - data is the simulated case-control data obtained by simUnmatched()
  # - N is number of Monte Carlo replicates
  # Output: Alpha_star
  
  # model=glm(case~.,data=data,family=binomial(link="logit")) 
  # alpha_star_initial=model$coefficients[1]
  
  AlphaStar=numeric()
  AlphaStar[1]=0
  AlphaStar[2]=-5
  
  p=2
  threshold=1E-04
  data=data.matrix(data)
  
  Weight=function(beta) {
    prod(exp(data[,K+1]*(AlphaStar[p]+data[,1:K]%*%beta))/(1+exp(AlphaStar[p]+data[,1:K]%*%beta)))
  }
  
  Y=rep(data[,K+1],times=N)
  betas=matrix(data=log(rf(K*N,m,m)),nrow=K,ncol=N)
  
  O=numeric() # offsets
  for (j in 1:N) {
    O=c(O,data[,1:K]%*%betas[,j])
  }
  
  while(abs(AlphaStar[p]-AlphaStar[p-1])>=threshold) { 
    W_t=numeric() # weights
    for (j in 1:N) {
      W_t[j]=Weight(betas[,j])
    }
    W=rep(W_t,each=dim(data)[1])
    
    g=glm(Y~offset(O),weights=W,family=binomial(link="logit"))
    #cat("EM iteration",p-1,":",g$coefficients,"\n")
    p=p+1
    AlphaStar[p]=g$coefficients
  }
  return(last(AlphaStar))
}


lkhdk=function(alpha_star,data,m,N) {
  # Input: 
  # - alpha_star is the output of MCEM 
  # - data is the simulated case-control data obtained by simKnockoffGenotypes()
  # - m is the value of m
  # - N is number of Monte Carlo replicates
  # Output: Monte Carlo estimate of the profile likelihood
  
  data=data.matrix(data)
  betas=matrix(data=log(rf(K*N,m,m)),nrow=K,ncol=N)
  lvec=rep(NA,N)
  for(j in 1:N) {
    lvec[j]=prod(exp(data[,K+1]*(alpha_star+data[,1:K]%*%betas[,j]))/(1+exp(alpha_star+data[,1:K]%*%betas[,j])))
  }
  return(log(mean(lvec)))
}  


profilelkhd=function(data,mvals,N) {
  # Input: 
  # - data is the simulated case-control data obtained by simKnockoffGenotypes()
  # - mvals is a set of values of m
  # - N is number of Monte Carlo replicates
  # Output: profile likelihood of m
  
  ll=rep(NA,length(mvals))
  for(m in 1:length(mvals)) {
    cat("Estimating profile log-likelihood for m =",m,"\n")
    alpha_star=MCEM(m,data,N)
    ll[m]=lkhdk(alpha_star,data,m,N)
  }
  return(ll)
}


#save(simKnockoffGenotypes,MCEM,lkhdk,profilelkhd,file="MCEM_Multiple-SNP.RData")




