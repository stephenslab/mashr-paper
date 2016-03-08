#'@title mixEM
##'   Fits a k component mixture model \deqn{f(x|\pi) = \sum_k \pi_k f_k(x)} to independent
##' and identically distributed data \eqn{x_1,\dots,x_n}. 
##' Estimates posterior on mixture proportions \eqn{\pi} by Variational Bayes, 
##' with a Dirichlet prior on \eqn{\pi}. 
##' Algorithm adapted from Bishop (2009), Pattern Recognition and Machine Learning, Chapter 10.
##' 
#' @param matrix_lik a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi.init the initial value of the posterior parameters. If not specified defaults to the prior parameters.
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 
#' 
#' @return A list, whose components include point estimates (pihat), 
#' the parameters of the fitted posterior on \eqn{\pi} (pipost),
#' the bound on the log likelihood for each iteration (B)
#' and a flag to indicate convergence (converged).
#'  
#' @export

mixEM = function(matrix_lik,prior,pi.init=NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  k=dim(matrix_lik)[2]
  if(is.null(pi.init)){
    pi.init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi.init,fixptfn=fixpoint, objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, 
              niter = res$iter, converged=res$convergence))
}

# helper functions used by mixEM
normalize = function(x){return(x/sum(x))}

fixpoint = function(pi, matrix_lik, prior){  
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
}



negpenloglik = function(pi,matrix_lik,prior){return(-penloglik(pi,matrix_lik,prior))}

penloglik = function(pi, matrix_lik, prior){
    pi = normalize(pmax(0,pi))
    m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik = sum(log(m.rowsum))
    subset = (prior != 1.0)
    priordens = sum((prior-1)[subset]*log(pi[subset]))
    return(loglik+priordens)
}


#The kth element of this vector is the derivative 
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with respect to $\pi_0$ at $\pi_0=1$.
gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - colSums(matrix_lik/matrix_lik[,1]) 
  return(grad)
}

##for when bma is used and some of the likelihoods are 0

mixEMbmaonly= function(matrix_lik,prior,pi.init=NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  k=dim(matrix_lik)[2]
  if(is.null(pi.init)){
    pi.init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi.init,fixptfn=fixpoint.play,objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, 
              niter = res$iter, converged=res$convergence))
}



#' @title compute.hm.train.log.lik
#' @description = takes a matrix of training sumamry statistics and their standard errors and computes the likelihood matrix according to a list of covariance matrices, using the exponent of the log likelihood - max (llj)
#' @param train.b =  JxR matrix of training beta hats
#' @param se.train = JxR matrix of training standard errors
#' @param covmat = LxK dimenstional (unlisted list) of prior covariance matrices
#' @param  A  output file name
#' @return An object containing pis and model fit from the EM, and a pdf barplot
#' @export



compute.hm.train.log.lik=function(train.b,se.train,covmat,A,pen=FALSE){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){log.lik.func(b.mle=train.b[x,],V.gp.hat=diag(se.train[x,])^2,covmat)}))##be sure to use log lik function
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  #train=lik.mat###but this matrix is the loglik matrix
  train=t(apply(lik.mat,1,function(x){
    e=exp(x-max(x))
    return(e)}
  ))
  #pis=mixEM.normlik(matrix_lik=train,prior=rep(1,ncol(train)))##here the matrix_lik is log normalized
  if(pen==TRUE){pis=mixEM(matrix_lik=train,prior=c(rep(1,ncol(train)-1),10))}
  pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}


##Here the matrix is normal likelihoods
mixEM.normlik= function(matrix_lik,prior,pi.init=NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  k=dim(matrix_lik)[2]
  if(is.null(pi.init)){
    pi.init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi.init,fixptfn=fixpoint.norm.lik,objfn=negpenloglik.normnew,matrix_lik=matrix_lik, prior=prior, control=controlinput)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, 
              niter = res$iter, converged=res$convergence))
}

fixpoint.norm.lik = function(pi, matrix_lik, prior){  
  pi = normalize(pmax(0,pi))
  log.pi=log(pi)#avoid occasional problems with negative pis due to rounding
  matrix.log.lik=matrix_lik
  
  
  
  
  classprob=t(apply(matrix.log.lik,1,function(x){
    log.lik.minus.max=x-max(x)
    s=log.lik.minus.max+log.pi
    exp.vec=exp(s)
    classprob=exp.vec/sum(exp.vec)
    return(classprob)}))
  
#   log.lik.minus.max=t(apply(matrix.log.lik,1,function(x){x-max(x)}))
#   log.pi=log(pi)
#   s=log.lik.minus.max+log.pi
#   exp.vec=(exp(s))
#   classprob=exp.vec/rowSums(exp.vec)
  
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
}

negpenloglik.normnew = function(pi,matrix_lik,prior){return(-penloglik.normnew(pi,matrix_lik,prior))}

penloglik.normnew = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi))
  log.pi=log(pi)
  matrix.log.lik=matrix_lik
  
  max.log.lik=apply(matrix.log.lik,1,function(x){max(x)})
  m.rowsum.logmethod=t(apply(matrix.log.lik,1,function(x){
    max.log.lik=max(x)
    log.lik.minus.max=x-max.log.lik
    s=log.lik.minus.max+log.pi
    exp.vec=exp(s)
    m.rowsum.logmethod=exp(max.log.lik)*sum(exp.vec)
    return(m.rowsum.logmethod)}))
  loglik.logmethod = sum(log(m.rowsum.logmethod))
  
  
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik.logmethod+priordens)
}