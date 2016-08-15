
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
  library('clusterGeneration')
    covmat=lapply(seq(1:K),function(k){
   A=genPositiveDefMat("eigen",dim=d)$Sigma
    A/max(diag(A))##scale so max diag is 1
    })
                    
  mus=rnorm(n)  ###generate a list of n mus
  mumat=matrix(rep(mus,d),ncol=d) ##generate a matrix of mus for each gene
   beta=t(sapply(seq(1:J),function(j){
    k=z[j]
    omega=abs(rnorm(1,mean=0,sd=betasd))##effect size variance can be big or small
    mvrnorm(1,mu=rep(0,d),Sigma=omega*covmat[[k]])
  }))
  
  beta=rbind(beta,matrix(rep(0,(n-J)*d),ncol=d))
  c=beta+mumat
  sj=abs(matrix(rnorm(n*d,esd,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sj,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  chat=c+e
t=chat/sj
 return(list(beta=beta,chat=chat,covmat=covmat,components=z,t=t,mumat=mumat,shat=sj,error=e,ceff=c))
}



