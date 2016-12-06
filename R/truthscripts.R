#'@title compute.covmat.using.truth
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standarderrors
#'@export

compute.covmat.using.truth = function(b.gp.hat,sebetahat,A,zero=FALSE,covlist,maxp=1,minp=1,power=1){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat,maxp,minp)
  omega.table=data.frame(omega^power)
  
  
  U.0=get.prior.covar.with.truth(R = ncol(b.gp.hat),omega.table = omega.table,covlist)
  covmat=unlist(U.0,recursive=F)
  if(zero==TRUE){covmat=c(covmat,list(matrix(rep(0,R*R),ncol=R,nrow=R)))}
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}




#' @title get.prior.covar.with.truth
#' @param R number of tissues
#' @param omega.table L vector grid weights
#' @return a L x K list of covariance matrices
#' @export
get.prior.covar.with.truth <- function(R,omega.table,covlist)  {
  test=list()
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    
    test[[l]][[1]]=omega*diag(1,R)
    for(c in 1:length(covlist)){
      test[[l]][[c+1]]=omega*covlist[[c]]}}
  
  return(U.0kl=test)
}



#'@title compute.covmat.using.truth.fixed.omega
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standarderrors
#'@export

compute.covmat.using.truth.fixed.omega = function(covlist,A,omega,zero=TRUE){
  
  test=list()
  R=dim(covlist[[1]])[1]
  test[[1]]=diag(1,R)
  for(c in 1:length(covlist)){
test[[c+1]]=omega*covlist[[c]]}
if(zero==TRUE){covmat=c(test,list(matrix(rep(0,R*R),ncol=R,nrow=R)))}
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}


#'@title compute.covmat.using.fixed.omega.withQ
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standarderrors
#'@export

compute.covmat.using.fixed.omega.withQ=function(max.step,Q,bma=TRUE,omega=2,zero=TRUE,A,power=1)  {
  test=list()
  test[[1]]=omega*diag(1,R)
    data.prox=max.step$true.covs[1,,]
    d.norm=data.prox/max(diag(data.prox))
    
    
    test[[2]]=omega*d.norm
    
    
    
    cov.pc=max.step$true.covs[3,,]
    
    
    cov.pc.norm=cov.pc/max(diag(cov.pc))
    
    
    
    test[[3]]=omega*(cov.pc.norm)
    if(Q!=0){for(q in 1:Q){
      
      
      rank.prox=max.step$true.covs[3+q,,]
      a=t(rank.prox)%*% rank.prox
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[q+3]]=omega*a.norm
    }}
    full.rank=max.step$true.covs[2,,]
    b.norm=full.rank/max(diag(full.rank))
    test[[Q+4]]=omega*b.norm
    
    if(bma==TRUE){
      configs=matrix(0,nrow=R,ncol=R)
      for(r in 1:R){
        configs[r,r]=1}
      
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        
        mat=(configs[c,]%*%t(configs[c,]))
        
        test[[Q+4+c]]=omega*mat}}
  
  covmat=test
  
  if(zero==TRUE){
    z=matrix(rep(0,R*R),ncol=R,nrow=R)
    covmat=c(covmat,list(z))
  }
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(list(covmat=covmat,omega))

}

#'@title compute.hm.mash.covmat.all.max.step.fixed.omega
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standarderrors
#'@export

compute.hm.mash.covmat.all.max.step.fixed.omega=function(bma=TRUE,omega,Q,lambda.mat,A,factor.mat,max.step,zero=FALSE,maxp=1,minp=1,power=1){
  
  
  
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  test=list()
  R=dim(max.step$true.covs)[3]
  test=list()
    
  test[[1]]=omega*diag(1,R)
    data.prox=max.step$true.covs[1,,]
    d.norm=data.prox/max(diag(data.prox))
    
    
    test[[2]]=omega*d.norm
    
    
    
    cov.pc=max.step$true.covs[3,,]
    
    
    cov.pc.norm=cov.pc/max(diag(cov.pc))
    
    
    
    test[[3]]=omega*(cov.pc.norm)
    if(Q!=0){for(q in 1:Q){
      
      load=as.matrix(lambda.mat[,q])
      fact=as.matrix(factor.mat[q,])
      rank.prox=load%*%t(fact)
      a=(t(rank.prox)%*% rank.prox)
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[q+3]]=omega*a.norm
    }}
    full.rank=max.step$true.covs[2,,]
    b.norm=full.rank/max(diag(full.rank))
    test[[Q+4]]=omega*b.norm
    
    if(bma==TRUE){
      configs=matrix(0,nrow=R,ncol=R)
      
      R=ncol(factor.mat)
      for(r in 1:R){
        configs[r,r]=1}
      
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        
        mat=(configs[c,]%*%t(configs[c,]))
        
        test[[Q+4+c]]=omega*mat}}
  covmat=test
  
  if(zero==TRUE){
    z=matrix(rep(0,R*R),ncol=R,nrow=R)
    covmat=c(covmat,list(z))
  }
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(list(covmat=covmat,omega))}




#'@title compute.covmat.using.oneK.fixed.omega
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standarderrors
#'@export

compute.covmat.using.oneK.fixed.omega = function(covlist,A,omega,zero=TRUE){
  
  test=list()
  R=dim(covlist[[1]])[1]
  test[[1]]=omega*diag(1,R)
  for(c in 1:length(covlist)){
    xnorm=covlist[[c]]/max(covlist[[c]])
    test[[c+1]]=omega*xnorm}
  if(zero==TRUE){covmat=c(test,list(matrix(rep(0,R*R),ncol=R,nrow=R)))}
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}


