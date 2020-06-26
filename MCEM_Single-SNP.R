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
  data=data[,cov_ind]
  data$case=NULL
  for (i in 1:dim(data)[1]) {
    s=sum(data[i,]*Beta)
    prob=exp(-2+s)/(1+exp(-2+s))
    if (prob >= 0.5) {data$case[i]=1}
    else {data$case[i]=0}
  }
  
  # Sample case/control status
  # case-control ratio = 1:4
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
  
  model=glm(case~.,data=data,family=binomial(link="logit")) 
  alpha_star_initial=model$coefficients[1]
  
  AlphaStar=numeric()
  AlphaStar[1]=0
  AlphaStar[2]=alpha_star_initial
  
  p=2
  threshold=1E-04
  
  # Weight=function(beta,con,case) { 
  #   sum(-log(1+exp(AlphaStar[p]+con*beta)))+sum(AlphaStar[p]+case*beta-log(1+exp(AlphaStar[p]+case*beta)))
  # }
  Weight=function(beta) { 
    prod(exp(data_k$case*(AlphaStar[p]+data_k$X*beta))/(1+exp(AlphaStar[p]+data_k$X*beta)))
  }
  
  # caseX=data %>% filter(case=="1") %>% select(-case)
  # conX=data %>% filter(case=="0") %>% select(-case)
  Y=rep(data$case,times=N)
  betas=log(rf(N,m,m))
  
  O=numeric() # offset
  for (j in 1:N) {
    O=c(O,data$X*betas[j])
  }
  
  while(abs(AlphaStar[p]-AlphaStar[p-1])>=threshold) { 
    W_t=numeric() # weight
    for (j in 1:N) {
      W_t[j]=Weight(betas[j])  
    }
    W=rep(W_t,each=dim(data)[1])
    
    g=glm(Y~offset(O),weights=W,family=binomial(link="logit"))
    cat("EM iteration",p-1,":",g$coefficients,"\n")
    p=p+1
    AlphaStar[p]=g$coefficients
  }
  return(last(AlphaStar))
}


lkhdk=function(alpha_k,data,m,N) {
  # Input: 
  # - alpha_k is the output of MCEM 
  # - data is the simulated case-control data obtained by simUnmatched()
  # - m is the value of m
  # - N is number of Monte Carlo replicates
  # Output: Monte Carlo estimate of the profile likelihood
  
  betas=log(rf(N,m,m))
  lvec=rep(NA,N)
  for(j in 1:N) {
    lvec[j]=prod(exp(data$case*(alpha_k+data$X*betas[j]))/(1+exp(alpha_k+data$X*betas[j])))
  }
  return(log(mean(lvec)))
} 


profilelkhd=function(data,mvals,N) {
  # Input: 
  # - data is the simulated case-control data obtained by simUnmatched()
  # mvals is a set of values of m
  # - N is number of Monte Carlo replicates
  # Output: profile likelihood of m
  
  ll=rep(NA,length(mvals))
  for(m in mvals) {
    cat("Estimating profile log-likelihood for m =",m,"\n")
    ll[m]=0
    for(k in 1:K) {
      data_k=data[,c(K+1,k)] # first K covars, then the case
      names(data_k)=c("case","X")
      alpha_k=MCEM(m,data_k,N) # estimating alpha_k by MCEM
      ll[m]=ll[m]+lkhdk(alpha_k,data_k,m,N)
    }
  }
  return(ll)
}


save(simKnockoffGenotypes,MCEM,lkhdk,profilelkhd,file="MCEM_Single-SNP.RData")

