#' @title autoselect.omega
#' @param betahat = column of betahats (or zs)
#' @param sebetahat= column of standard errs
#' @param minparam = what the standard error is divided by (want to be small)
#' @param mult = density parameter (intervals between points)
#' @param maxparam = how big the grid needs to be
#' @concept this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2 
#' @return vector of variances for one tissue

autoselect.omega=function(betahat,sebetahat,minparam,mult,maxparam){
  ##convert autoselected sigmas to reasonable choice of variances
  sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
  sigmaamin = (min(sebetahat)/minparam)^2 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = maxparam*(max(betahat^2-sebetahat^2)) 
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}


#' @title mult.tissue.omega
#' @param betahat = JxR matrix of betahats (or zs)
#' @param sebetahat= JxR matrix of standard errs
#' @param minparam = what the standard error is divided by (want to be small)
#' @param multone = density parameter (intervals between points for within tissue variances, default 2 so that the ratio between sd is sqrt(2))
#' @param maxparam = how big the grid needs to be
#' @concept this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2 
#' @return vector of variances for all tissues


mult.tissue.omega=function(mult,betahat,sebetahat,minparam,multone,maxparam){
  
R=ncol(betahat);
mix.weights=unlist(sapply(seq(1:ncol(betahat)),function(r){
  
  autoselect.omega(betahat = betahat[,r],sebetahat = sebetahat[,r],minparam=minparam,mult=multone,maxparam=maxparam)}))
sigmaamin=min(mix.weights);sigmaamax=max(mix.weights);
npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult));
omega=mult^((-npoint):0) * sigmaamax;return(omega)}



#' @title compute.hm.covmat.all.max.step.omega
#' @param betahat = JxR matrix of betahats (or zs)
#' @param sebetahat= JxR matrix of standard errs
#' @param minparam = what the standard error is divided by (want to be small)
#' @param multone = density parameter (intervals between points for within tissue variances, default 2 so that the ratio between sd is sqrt(2))
#' @param maxparam = how big the grid needs to be
#' @param Q number of factors to consider
#' @param t.stat matrix of strong t statistics (needs to be R in columns but not necessarily j in rows)
#' @param factor mat a fxR matrix of factors
#' @param BMA=true wheter or not to use singleton and shared config
#' @param max.step = object with denoised covariance matrices
#' @return A list of covariance matrices







compute.hm.covmat.all.max.step.omega=function(mult=sqrt(2),minparam=4,b.hat,se.hat,t.stat,v.j,Q,lambda.mat,A,factor.mat,max.step,multone,maxparam)
  {
  X.real=as.matrix(t.stat)
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  R=ncol(X.c)
  omega=mult.tissue.omega(mult,b.hat,se.hat,minparam,multone,maxparam)
  omega.table=data.frame(omega)
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  U.0kl=get.prior.covar.with.all.max.step(X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
  covmat=unlist(U.0kl,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  return(covmat)
}