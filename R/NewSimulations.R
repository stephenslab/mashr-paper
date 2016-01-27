#'@title factor.sim.new
#'@description Generate J beta,betahats, T statistics across R tissues from a set of covariance matrices
#'@description Creates factors with magnitude based on betasd - each betaj is loaded on a minimal number of factors by simulating from rmvnorm
#'@param J Number of Gene-SNP Pairs (by definition, only 0.8% of them will be true eQTL; pi0 is 0.20 and the number of genes is 1/100 J)
#'@param d Number of subgroups
#'@param betasd size of covariance of true effects
#'@param esd standard error
#'@details: to simulate strong sharing, I use the gtex covariance matrices
#'@export

factor_sim_new=function(J,d=44,betasd=1,esd=0.1,tspec=0){
  n=trunc(0.008*J,units = 0)##number of significant gene-snp Pairs, so there are 100 snps in cis of a gene and one causal snp
  #F=t(sapply(seq(1:K),function(x){rnorm(d,mean=0,sd=betasd)})) 
  covmat=readRDS("~/Dropbox/Aug12/covmatAug13withED.rds")[2:9]
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
    omega=abs(rnorm(1,mean=0,sd=betasd^2))##effect size can be big or small
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
  tstat=betahat/abs(sgp)
  lambda=rbind(possible_loadings[z,],matrix(rep(0,(J-n)*K),ncol=K))+load.e
  return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=sgp,tstat=tstat,component.id=z))
}

