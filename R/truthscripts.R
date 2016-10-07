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




