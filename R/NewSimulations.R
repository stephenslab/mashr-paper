#' @details Generate J beta,betahats, T statistics across R tissues from a set of covariance matrices
#' @param Creates factors with magnitude based on betasd where each betaj is loaded on a minimal number of factors by simulating from rmvnorm
#' @param J Number of Gene-SNP Pairs and by definition, only 0.8 percent of them will be true eQTL; pi0 is 0.20 and the number of genes is J/100)
#' @param d Number of subgroups
#' @param n number of real associations
#' @param betasd size of covariance of true effects
#' @param esd standard error
#' @details to simulate strong sharing, I use the gtex covariance matrices

#' @title factor.sim.new
#' @export

factor_sim_new=function(J,d=44,betasd=1,esd=0.11,tspec=0,n=400,R=44){
  #n=trunc(0.008*J,units = 0)##number of significant gene-snp Pairs, so there are 100 snps in cis of a gene and one causal snp
  n=n
  #F=t(sapply(seq(1:K),function(x){rnorm(d,mean=0,sd=betasd)})) 
  #covmat=readRDS("~/matrix_ash/inst/simdata/covmatforsimulation.rds")[2:9]
  covmat=readRDS(system.file('simdata/covmatforsimulation.rds', package = 'mashr'))[2:9]
  covmat=lapply(seq(1:length(covmat)),function(x){covmat[[x]]/max(diag(covmat[[x]]))})
 
  if(tspec!=0){
    configs=matrix(0,nrow=R,ncol=R)
    for(r in 1:R){
      configs[r,r]=1}
    specs=sample(seq(1:d),tspec,replace=FALSE)##generate tissue specific configs if required
    singleton.mat=list()
    for(t in 1:length(specs)){
      tissue=specs[t]
      singleton.mat[[t]]=(configs[tissue,]%*%t(configs[tissue,]))}
    covmat=append(covmat,singleton.mat)
  }
  
  library("mvtnorm")
  library("MASS")
  K=length(covmat)
  possible_loadings = diag(K) #set possible loadings to be "sparse" (loaded on one factor each)
  
  z = sample(K,n,replace=TRUE)
  
  beta=t(sapply(seq(1:n),function(j){
    k=z[j]
    omega=abs(rnorm(1,mean=0,sd=betasd))##effect size variance can be big or small
    mvrnorm(1,mu=rep(0,d),Sigma=omega*covmat[[k]])
  }))
  beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))
  sj=abs(matrix(rnorm(J*d,0.11,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  load.e=matrix(rnorm(J*K,0,0.001),ncol=K)##to make sure all loadings are covered
  #sign.index=matrix(rbinom(n*d,size=1,prob=0.5)+1,ncol=d)
  #sign.choice=c(-1,1)
  #sign.mat=matrix(sign.choice[sign.index],ncol=d)
  #e=sign.mat*sebetahat##make some negative and some positive
  betahat = rbind(beta + e)
  #betahat=rbind(betahat,)
  tstat=betahat/abs(sj)
  lambda=rbind(possible_loadings[z,],matrix(rep(0,(J-n)*K),ncol=K))+load.e
  return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=sj,tstat=tstat,component.id=z))
}


#' @title independent.simulation
#' @export

independent.data=function(J,d=44,betasd=1,esd=0.1,tspec=0){
  n=trunc(0.008*J,units = 0)##number of significant gene-snp Pairs, so there are 100 snps in cis of a gene and one causal snp
  beta=matrix(rnorm(d*n,mean=0,sd=betasd),ncol=d,nrow=n) 
  beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))
  sj=abs(matrix(rnorm(J*d,0.11,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  library("mvtnorm")
  library("MASS")
  betahat = rbind(beta + e)
  
  tstat=betahat/abs(sj)
  
  return(list(beta=beta,betahat=betahat,sebetahat=sj,tstat=tstat))
}


#' @title independent.from.omega
#' @details simulate shared effects with randomly chosen omegas
#' @export

independent.from.omega=function(J,d=44,betasd=1,esd=0.1,tspec=0){
  n=trunc(0.008*J,units = 0)##number of significant gene-snp Pairs, so there are 100 snps in cis of a gene and one causal snp
  betasd=c(0.1,0.5,0.75,1)
  beta=matrix(rnorm(d*n,mean=0,sd=sample(betasd,d*n,replace = T)),ncol=d,nrow=n) 
  beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))
  sj=abs(matrix(rnorm(J*d,0.11,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  library("mvtnorm")
  library("MASS")
  betahat = (beta + e)
  
  tstat=betahat/abs(sj)
  
  return(list(beta=beta,betahat=betahat,sebetahat=sj,tstat=tstat))
}

#' @title independent.totally
#' @export

independent.totally=function(J,d=44,betasd=1,esd=1,tspec=0){
  n=trunc(0.02*J,units = 0)##number of significant gene-snp Pairs, so there are 100 snps in cis of a gene and one causal snp
  betasd=c(0.1,0.5,0.75,1)
  beta=matrix(rep(0,J*d),ncol=d)
  sim.sd=matrix(rep(0,J*d),ncol=d)
for(r in 1:d){
    true.effects=sample(seq(1:J),n,replace=FALSE)
    bsd=sample(betasd,n,replace = T)
    sim.b=rnorm(n,mean=0,sd=bsd)
    beta[true.effects,r]=sim.b
    sim.sd[true.effects,r]=bsd
}
  
  sj=abs(matrix(rnorm(J*d,0.11,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rnorm(d,mean=0,sd=x)}))
  library("mvtnorm")
  library("MASS")
  betahat = beta + e
  
  tstat=betahat/sj
  
  return(list(beta=beta,betahat=betahat,sebetahat=sj,tstat=tstat))
}

#' @title sim.with.error
#' @export

sim.with.error=function(J,d=44,betasd=1,esd=0.11,n=400,rho=0.8){
   n=n
  covmat=readRDS(system.file('simdata/covmatforsimulation.rds', package = 'mash'))[2:9]
  covmat=lapply(seq(1:length(covmat)),function(x){covmat[[x]]/max(diag(covmat[[x]]))})
  
  
  
  library("mvtnorm")
  library("MASS")
  K=length(covmat)
  
  
  if(n!=0){
  z = sample(K,n,replace=TRUE)
  omega=abs(rnorm(n,mean=0,sd=betasd))##effect size variance can be big or small
  beta=t(sapply(seq(1:n),function(j){
    k=z[j]
    o=omega[j]
    mvrnorm(1,mu=rep(0,d),Sigma=o*covmat[[k]])
  }))
  beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))}
  if(n==0){
    beta=matrix(rep(0,(J-n)*d),ncol=d)
  }
  
  s.j.r=abs(rnorm(d,esd,0.001))##simulate with the same standard error for every J
  cormat=rep(1,d)%*%t(rep(1,d))*rho+(1-rho)*diag(1,d)#make the errors correlated
  v.mat=diag(s.j.r)%*%cormat%*%diag(s.j.r)##now v.j.r will be the same for every J
  cormat=cov2cor(v.mat)
  e=rmvnorm(J,mean=rep(0,d),sigma=v.mat)
  betahat = beta + e
  s.j=matrix(rep(as.matrix(s.j.r)),nrow(betahat),byrow=T,ncol=d)
  t.stat=betahat/abs(s.j)
  if(n!=0){
    return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=s.j,t.stat=t.stat,component.id=z,error=e,var.mat=v.mat,omega=omega,cormat=cormat))
  }
  if(n==0){
    return(list(beta=beta,betahat=betahat,sebetahat=s.j,t.stat=t.stat,error=e,var.mat=v.mat))
  }
}


#' @title gtexchatsim
#' @export

gtexchatsim=function(J,d=44,betasd=1,esd=0.11,tspec=0,n=400){
  #n=trunc(0.008*J,units = 0)##number of significant gene-snp Pairs, so there are 100 snps in cis of a gene and one causal snp
  n=n
  #F=t(sapply(seq(1:K),function(x){rnorm(d,mean=0,sd=betasd)})) 
  #covmat=readRDS("~/matrix_ash/inst/simdata/covmatforsimulation.rds")[2:9]
  covmat=readRDS(system.file('simdata/covmatforsimulation.rds', package = 'mash'))[2:9]
  covmat=lapply(seq(1:length(covmat)),function(x){covmat[[x]]/max(diag(covmat[[x]]))})
  mus=rnorm(J)  ###generate a list of n mus
  mumat=matrix(rep(mus,d),ncol=d) ##generate a matrix of mus for each gene
  
  if(tspec!=0){
    configs=matrix(0,nrow=R,ncol=R)
    for(r in 1:R){
      configs[r,r]=1}
    specs=sample(seq(1:d),tspec,replace=FALSE)##generate tissue specific configs if required
    singleton.mat=list()
    for(t in 1:length(specs)){
      tissue=specs[t]
      singleton.mat[[t]]=(configs[tissue,]%*%t(configs[tissue,]))}
    covmat=append(covmat,singleton.mat)
  }
  
  library("mvtnorm")
  library("MASS")
  K=length(covmat)
  
  z = sample(K,n,replace=TRUE)
  
  beta=t(sapply(seq(1:n),function(j){
    k=z[j]
    omega=abs(rnorm(1,mean=0,sd=betasd))##effect size variance can be big or small
    mvrnorm(1,mu=rep(0,d),Sigma=omega*covmat[[k]])
  }))
  beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))
  sj=abs(matrix(rnorm(J*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  mus=rnorm(J)  ###generate a list of n mus
  mumat=matrix(rep(mus,d),ncol=d) ##generate a matrix of mus for each gene
  c=beta+mumat

  chat=c+e
  t=chat/sj
  #betahat=rbind(betahat,)
  tstat=chat/abs(sj)
  return(list(beta=beta,chat=chat,covmat=covmat,components=z,t=t,mumat=mumat,shat=sj,error=e,ceff=c))
  #return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=sj,tstat=tstat,component.id=z))
}






