
autoselect.omega=function(betahat,sebetahat,mult){
  sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
  sigmaamin = (min(sebetahat)/4)^2 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = (max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}



mult.tissue.omega=function(mult=sqrt(2),betahat,sebetahat){
  
R=ncol(betahat);
mix.weights=unlist(sapply(seq(1:ncol(betahat)),function(r){
  autoselect.omega(betahat = betahat[,r],sebetahat = sebetahat[,r],mult)}))
sigmaamin=min(mix.weights);sigmaamax=max(mix.weights);
npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult));
omega=mult^((-npoint):0) * sigmaamax;return(omega)}

compute.hm.covmat.all.max.step.omega=function(mult=sqrt(2),b.hat,se.hat,t.stat,v.j,Q,lambda.mat,A,factor.mat,max.step)
  {
  X.real=as.matrix(t.stat)
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  R=ncol(X.c)
  omega=mult.tissue.omega(mult=sqrt(2),b.hat,se.hat)
  omega.table=data.frame(omega)
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  U.0kl=get.prior.covar.with.all.max.step(X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
  covmat=unlist(U.0kl,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  return(covmat)
}