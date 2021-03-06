
#' @title chat_sim
#' @param n number of SNPs
#' @param d number of subgroups
#' @param betasd the vriance of true effect
#' @param esd the standard error of E in chat=mu+c+E, i.e., E~N(0,diag(esd^2))
#' @export

chat_sim=function(n=1000,d=3,betasd=1,esd=0.1,K=10){
  J=0.10*n
  
#   configs = matrix((rnorm(d*K)),byrow=T,ncol=d) # A matrix of K classes (patterns) across R subgroups 
#   F=as.matrix(configs);ftf=lapply(seq(1:K),function(k){F[k,]%*%t(F[k,])}) ## each entry of F is the the factor of decomposition of covariance of effect sizes
#   cormat=lapply(seq(1:K),function(k){A=matrix(runif(d*d),nrow=d);cormat=A+t(A);diag(cormat)=rep(1,d);return(cormat)})
#   k = nrow(F) # number of factors
#   
#   possible_loadings = diag(k) #set possible loadings to be "sparse" (loaded on one factor each)
  z = sample(K,J,replace=TRUE) # randomly sample factor to be loaded on for each real snp
 
#   covmat=lapply(seq(1:K),function(k){
#     A=ftf[[k]]*cormat[[k]]##multiply by pairwise correlations;
#     A%*%t(A)##ensure positive semidefinite
#     A/max(diag(A))##scale so max diag is 1
#     })

  #RNGkind("Mersenne-Twister")##make sure not to mess with seed
    covmat=lapply(seq(1:K),function(k){
   A=genPositiveDefMat("eigen",dim=d)$Sigma
    A/max(diag(A))##scale so max diag is 1
    })
                    
  mus=rnorm(n)  ###generate a list of n mus
  mumat=matrix(rep(mus,d),ncol=d) ##generate a matrix of mus for each gene
   beta=t(sapply(seq(1:J),function(j){
    k=z[j]
    omega=abs(rnorm(1,mean=0,sd=betasd))##effect size variance can be big or small
    #mvrnorm(1,mu=rep(0,d),Sigma=omega*covmat[[k]])
    rmvnorm(1,mean = rep(0,d),sigma=omega*covmat[[k]])
  }))
  
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  c=beta+mumat
  sj=abs(matrix(rnorm(n*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
t=chat/sj
 return(list(beta=beta,chat=chat,covmat=covmat,components=z,t=t,mumat=mumat,shat=sj,error=e,ceff=c))
}



#' @title chat_sim_fact
#' @param n number of SNPs
#' @param d number of subgroups
#' @param betasd the vriance of true effect
#' @param esd the standard error of E in chat=mu+c+E, i.e., E~N(0,diag(esd^2))
#' @return omega the size of the effects
#' @export

chat_sim_fact=function(n=1000,d=44,betasd=1,esd=0.1,K=10){
  library("MASS")
  library("mvtnorm")
  J=0.10*n
  
  configs = matrix((rnorm(d*K)),byrow=T,ncol=d) # A matrix of K classes (patterns) across R subgroups 
  F=as.matrix(configs);
  covmat=lapply(seq(1:K),function(k){
    A=F[k,]%*%t(F[k,]);
    A/max(diag(A))})
  ## each entry of F is the the factor of decomposition of covariance of effect sizes
  z = sample(K,J,replace=TRUE) # randomly sample factor to be loaded on for each real snp
  
  
  mus=rnorm(n)  ###generate a list of n mus
  mumat=matrix(rep(mus,d),ncol=d)##generate a matrix of mus for each gene
  omega=abs(rnorm(J,mean=0,sd=betasd))##effect size variance can be big or small
  beta=t(sapply(seq(1:J),function(j){
    k=z[j]
    mvrnorm(1,mu=rep(0,d),Sigma=omega[j]*covmat[[k]])
    #rmvnorm(1,mean = rep(0,d),sigma=omega*covmat[[k]])
  }))
  
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  c=beta+mumat
  sj=abs(matrix(rnorm(n*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
  t=chat/sj
  return(list(beta=beta,chat=chat,covmat=covmat,components=z,t=t,mumat=mumat,shat=sj,error=e,ceff=c,F=F,omega=omega))
}



#' @title chat_sim_fsimple
#' @param n number of SNPs
#' @param d number of subgroups
#' @param betasd the vriance of true effect
#' @param esd the standard error of E in chat=mu+c+E, i.e., E~N(0,diag(esd^2))
#' @return omega the size of the effects
#' @return f the factors from which they were simulated
#' @export

chat_sim_fsimple=function(n=1000,d=8,betasd=1,esd=0.1,K=10){
  library("MASS")
  library("mvtnorm")
  J=0.10*n
  temp=rep(list(c(0,1)),d)
  configs = expand.grid(temp) # all possible 2^d combinations
  S=sample(seq(1:nrow(configs)),size = K,replace = FALSE)##which factors will be used
  F=as.matrix(configs[S,])
  covmat=lapply(seq(1:K),function(k){
    A=F[k,]%*%t(F[k,]);
    A/max(diag(A))})
  ## each entry of F is the the factor of decomposition of covariance of effect sizes
  z = sample(K,J,replace=TRUE) # randomly sample factor to be loaded on for each real snp
  
  
  mus=rnorm(n)  ###generate a list of n mus
  mumat=matrix(rep(mus,d),ncol=d)##generate a matrix of mus for each gene
  omega=abs(rnorm(J,mean=0,sd=betasd))##effect size variance can be big or small
  beta=t(sapply(seq(1:J),function(j){
    k=z[j]
    mvrnorm(1,mu=rep(0,d),Sigma=omega[j]*covmat[[k]])
    #rmvnorm(1,mean = rep(0,d),sigma=omega*covmat[[k]])
  }))
  
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  c=beta+mumat
  sj=abs(matrix(rnorm(n*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
  t=chat/sj
  return(list(beta=beta,chat=chat,covmat=covmat,components=z,factors=F,t=t,mumat=mumat,shat=sj,error=e,ceff=c,omega=omega))
}

#' @title chat_sim_fsingle
#' @param n number of SNPs
#' @param d number of subgroups
#' @param betasd the vriance of true effect
#' @param esd the standard error of E in chat=mu+c+E, i.e., E~N(0,diag(esd^2))
#' @return omega the size of the effects
#' @return f the factors from which they were simulated
#' @export

chat_sim_fsingle=function(n=1000,d=8,betasd=1,esd=0.1,K=10){
  library("MASS")
  library("mvtnorm")
  J=0.10*n
  config=rep(0,d)
  config[d]=1
  config[d-1]=1
  mus=rnorm(n)  ###generate a list of n mus
  covmat=(config)%*%t(config)##only active in tissues d and d-1
  mumat=matrix(rep(mus,d),ncol=d)##generate a matrix of mus for each gene
  omega=abs(rnorm(J,mean=0,sd=betasd))##effect size variance can be big or small
  beta=t(sapply(seq(1:J),function(j){
    
    mvrnorm(1,mu=rep(0,d),Sigma=omega[j]*covmat)}))
  
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  c=beta+mumat
  sj=abs(matrix(rnorm(n*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
  t=chat/sj
  return(list(beta=beta,chat=chat,covmat=covmat,t=t,mumat=mumat,shat=sj,error=e,ceff=c,omega=omega))
}


#' @title chat_sim_fsingle_fixedomega
#' @param n number of SNPs
#' @param d number of subgroups
#' @param betasd the vriance of true effect
#' @param esd the standard error of E in chat=mu+c+E, i.e., E~N(0,diag(esd^2))
#' @return omega the size of the effects
#' @return f the factors from which they were simulated
#' @export

chat_sim_fsingle_fixedomega=function(n=1000,d=8,omega=2,esd=0.1,mu=0,reals=0.10){
  library("MASS")
  library("mvtnorm")
  J=reals*n
  config=rep(0,d)
  config[d]=1
  config[d-1]=1
if (mu==FALSE) {
    mumat=matrix(rep(0,n*d),ncol=d)
  } else {
    mus=rnorm(n,mean = mu) 
    mumat=matrix(rep(mus,d),ncol=d)
  }
   ###generate a list of n mus
  covmat=(config)%*%t(config)##only active in tissues d and d-1
  ##generate a matrix of mus for each gene
  beta=t(sapply(seq(1:J),function(j){
    
    mvrnorm(1,mu=rep(0,d),Sigma=omega*covmat)}))
  
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  c=beta+mumat
  sj=matrix(rep(esd,n*d),ncol=d)##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
  t=chat/sj
  return(list(beta=beta,chat=chat,covmat=covmat,t=t,mumat=mumat,shat=sj,error=e,ceff=c,omega=omega))
}