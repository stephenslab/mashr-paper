#'@title mixEM
#'  @details Fits a k component mixture model \deqn{f(x|\pi) = \sum_k \pi_k f_k(x)} to independent
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



#The kth element of this vector is the derivative 
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with respect to $\pi_0$ at $\pi_0=1$.
gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - colSums(matrix_lik/matrix_lik[,1]) 
  return(grad)
}
