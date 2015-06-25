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


#' @title post.b.jk.cov
#' @return RxR posterior covariance matrix for a given prior covariance matrix
#' @export
post.b.jk.cov=function(V.j.hat.inv, U.0k){
  U.j1k <- U.0k %*% solve(V.j.hat.inv %*% U.0k + diag(nrow(U.0k)))
  return(U.j1k)
}


#' @title post.b.jk.mean
#' @param U.0k.l let U.0k.l represent a specific matrix in U.0kl (.e.g, U.0kl[[l]][[k]])
#' @return return a 1 * R vector of posterior means for a given prior covariance matrix
#' @export
post.b.jk.mean = function(b.mle, V.j.hat.inv, U.1jk){
  mu.j1k <- U.j1k %*% V.j.hat.inv %*% b.mle
  return(mu.j1k)}



tinv=solve(U.k+V.j.hat)

post.b.jk.ed.mean = function(b.mle, tinv,U.k){
b.jk=U.k%*%tinv%*%b.mle
  return(b.jk)}


post.b.jk.ed.cov = function(b.mle, tinv,U.k){
  B.jk=U.k-U.k%*%tinv%*%U.k
  return(B.jk)}

#' @title tarray
#' @export
tarray <- function(x) aperm(x, rev(seq_along(dim(x))))



#' @title em.array.generator
#' @param b.j.hat = 1xR vector of MLEs
#' @param se.j.jat=1xR vector of their standard errors
#' @param covmat = K dimensional list
#' @param lik.mat JxK matrix of likelihoods
#' @return list of JxKxR conditional posterior means, JxKxRxR covariance matries and JxK normalized likelihoods
#' @export



em.array.generator=function(b.j.hat,J,se.j.hat,true.covs,pi){
  R=ncol(b.j.hat)
  K=dim(true.covs)[1]
  post.means=array(NA,dim=c(J,K,R))
  post.covs=array(NA,dim=c(J,K,R,R))
  q.mat=array(NA,dim=c(J,K))
  
  for(j in 1:J){
    b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
    V.j.hat=diag(se.j.hat[j,]^2)
    #V.j.hat.inv <- diag(se.j.hat[j,]^-2)##to avoid having to 'solve' since we know that it is simply diag(1/s^2)
    lik=sapply(seq(1:K),function(k){dmvnorm(x=b.mle, sigma=true.covs[k,,] + V.j.hat)})##compute K element likeilihood for each idndiviual
    a=(lapply(seq(1:K),function(k){post.b.jk.ed.cov(b.mle,tinv=solve(true.covs[k,,]+V.j.hat),true.covs[k,,])}))##create a K dimensional list of covariance matrices 
    
    post.covs[j,,,] <- tarray(array(unlist(a), c( 44, 44,K))) ###store a K dimensional list of posterior covariances for each J (KxRxR)
    
    pm=(lapply(seq(1:K),function(k){t(post.b.jk.ed.mean(b.mle,tinv=solve(true.covs[k,,]+V.j.hat),true.covs[k,,]))}))##compute a K dimensional list of posterior means for each J
    post.means[j,,]=matrix(unlist(pm),ncol=44,byrow=TRUE) ##store as KxR matrix for each indiviudal
    q.mat[j,]=pi*lik/sum(pi*lik)##compute a k dimensional normalixed likelihood for each individual 
#     for(k in 1:K){
#       
#       U.k=true.covs[k,,]
#       tinv=solve(true.covs[k,,]+V.j.hat)
#       U.j1k <- post.b.jk.ed.cov(b.mle,tinv,U.k)
#       mu.j1k <- post.b.jk.ed.mean(b.mle,tinv,U.k)
#       post.means[j,k,]=mu.j1k
#       post.covs[j,k,,]=U.j1k ##critically, now store the actual matrix rather than just its diagonal
#       q.mat[j,]=pi*lik/sum(pi*lik)
#     }
  }
  return(list(post.means=post.means,post.covs=post.covs,q.mat=q.mat))
}


#' @title max.step.func 
#' @param post.means JxKxR matrix of posterior means
#' @param post.covs JxKxR matrix of posterior covariances
#' @param q.mat JxK matrix of normalixed likelihoods e.g., p(K|D)
#' @return K x R xR array of 'true covariance matrices' and Kx1 vector of pis
#' @export
max.step.func = function(post.means,post.covs,q.mat){
  J=dim(post.means)[1]
  #true.means=array(NA,dim=c(K,R)) but we force to be 0
  true.covs=array(NA,dim=c(K,R,R))
  q=colSums(q.mat) ##compute the sum of the posterior weights across individuals
  pis=q/J ##normalise
  d=array(NA,dim=c(J,R,R))
  for(k in 1:K){
    if(q[k]==0){
    #true.means[k,]=rep(0,R)
    true.covs[k,,]=array(rep(0,R*R),dim=c(R,R))}
  else{
  #true.means[k,]=1/q[k]*(q.mat[,k]*post.means[,k,])
  for(j in seq(1:J)){
      d[j,,]=q.mat[j,k]*(-post.means[j,k,]%*%t(-post.means[j,k,])+post.covs[j,k,,])##produce a RxR matrix of weighted 'truth' for each individual
  }
  true.covs[k,,]=1/q[k]*apply(d, MARGIN=c(2, 3), sum)}##get the latent identifier weighted sum of all J individuals to produce 1 RxR matrix of truth
  }
  return(list(true.covs=true.covs,pis=pis))
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
