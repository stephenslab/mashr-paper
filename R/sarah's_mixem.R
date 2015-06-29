library("mvtnorm")
###Recall that T is the variance of the marginal distribution of $hat{b}$ (i.e., integrating over the uncertainty in B)
#tinv=solve(U.k+V.j.hat)


#' @title post.b.jk.ed.mean
#' @param b.mle=Rx1 vector of mles
#' @param tinv = variance of marginal distirbution of bhat, RxR matrix (e.g.,tinv=solve(U.k+V.j.hat)) )
#' @param U.k = RxR prior covariance matrix
#' @return component specific vector of posterior means for the jth gene snp pair

post.b.jk.ed.mean = function(b.mle, tinv,U.k){
  b.jk=U.k%*%tinv%*%b.mle
  return(b.jk)}

#' @title lik.func.em
#' @param b.mle=Rx1 vector of mles
#' @param B.j.hat = matrix of squared standard erros
#' @param true.covs = KxRxR array of prior covariance matrices
#' @return 1xK vector of likelihoods for the jth gene snp pair across all K covariance matrices

lik.func.em=function(true.covs,b.mle,V.j.hat){sapply(seq(1:K),function(k){dmvnorm(x=b.mle, sigma=true.covs[k,,] + V.j.hat)})}

#' @title post.b.jk.ed.cov
#' @param tinv = variance of marginal distirbution of bhat, RxR matrix (e.g.,tinv=solve(U.k+V.j.hat)) )
#' @param U.k = prior covariance matrix (RxR)
#' @return RxR posterior covariance matrix for the Jth gene snp pair at the kth componenet

post.b.jk.ed.cov = function(tinv,U.k){
  B.jk=U.k-U.k%*%tinv%*%U.k
  return(B.jk)}

#' @title tarray
#' @details ensures that the unlisted k dimensional lists are properly sotred
#' @export
tarray <- function(x) aperm(x, rev(seq_along(dim(x))))


#' @title dim.true.cov.fun
#' @details outputs a list of dimensions for both the true covariance matrices and the vector of pis
#' @export


dim.true.cov.fun=function(b.j.hat,se.j,hat,max.step.unlist){
  L=length(max.step.unlist)
  K=L/(R^2+1)
  dim.true.covs=c(K,R,R)
  pi.length=K
  return(list(dim.true.covs=dim.true.covs,pi.length=pi.length))
}

#' @title init.covmat
#' @param t.stat matrix of strong t statistics (MxR) from which to derive covariance matrices
#' @param factor.mat KxR matrix of factors from SFA
#' @param lambda.mat KxR matric of loadings, also from SFA
#' @param K number of components to fit (i.e., the first dimension of the array)
#' @param P rank of SVD approximation
#' @return KxRxR matric of prior covariance matrices to initialize the EM
#' @export

init.covmat=function(t.stat=t.stat,factor.mat=factors,lambda.mat=lambda,K=3,P=2){
  K=K
  R=ncol(t.stat)#number of tissues
  true.covs=array(NA,dim=c(K,R,R))
  
  X.t=as.matrix(t.stat)
  X.real=X.t
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  
  M=nrow(X.c)
  data.prox=((t(X.c)%*%X.c)/M)
  true.covs[1,,]=data.prox
  full.rank=lambda.mat%*%factor.mat
  sfa.prox=(t(full.rank)%*%full.rank)/(M)
  true.covs[2,,]=sfa.prox
  
  svd.X=svd(X.c)
  v=svd.X$v;u=svd.X$u;d=svd.X$d
  cov.pc=1/M*v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P])%*%t(v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P]))
  true.covs[3,,]=cov.pc
  return(true.covs)
}



#' @title em.array.generator
#' @param b.j.hat = 1xR vector of MLEs
#' @param se.j.hat=1xR vector of their standard errors
#' @param max.step=list containing two arrays: true.covs= KxRxR arrays of prior covariance matrices
#' @param pi K vector of prior weights
#' @return list of 3 arrays: JxKxR conditional posterior means, JxKxRxR posterior covariance matries and JxK normalized likelihoods (P(K|Data))
#' @export



em.array.generator=function(max.step,b.j.hat,se.j.hat){
  true.covs=max.step$true.covs
  pi=max.step$pi
  R=ncol(b.j.hat)
  K=dim(true.covs)[1]
  J=nrow(b.j.hat)
  
  post.means=array(NA,dim=c(J,K,R))
  post.covs=array(NA,dim=c(J,K,R,R))
  q.mat=array(NA,dim=c(J,K))
  
  
  for(j in 1:J){
    b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
    V.j.hat=diag(se.j.hat[j,]^2)
    lik=lik.func.em(true.covs,b.mle,V.j.hat)
    ##compute a K dimensional list of posterior covariance for each J, this is faster than for looping over the K, but in the main function we have to do within R also
    B.j.=lapply(seq(1:K),function(k){
      
      tinv=solve(true.covs[k,,]+V.j.hat)##covariance matrix of the marginal distribution
      post.b.jk.ed.cov(tinv=tinv,true.covs[k,,])
        
    }
    )##create a K dimensional list of covariance matrices 
    post.covs[j,,,] <- tarray(array(unlist(B.j.), c( R, R,K))) ###store a K dimensional list of posterior covariances for each J (JxKxRxR) in the post.means array
    ##this is different than in my maincode, where I forloop through each K, because there i need to compute the tail probabilities for each r within a k and so I would have to forloop inside K which is ugly
    
    ##compute a K dimensional list of posterior means for each J, again this is faster than for looping over the K, but in the main function we have to do within R also
   b.j.=(lapply(seq(1:K),function(k){
      tinv=solve(true.covs[k,,]+V.j.hat)##covariance matrix of the marginal distribution
      b=post.b.jk.ed.mean(b.mle,tinv=tinv,true.covs[k,,])##for each component, compute posterior mean
      
    }
    ))
    post.means[j,,]=matrix(unlist(b.j.),ncol=R,byrow=TRUE) ##store as KxR matrix for each indiviudal
    q.mat[j,]=pi*lik/sum(pi*lik)##compute a k dimensional normalized likelihood for each individual 

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
  K=dim(post.means)[2]
  R=dim(post.means)[3]
  #true.means=array(NA,dim=c(K,R)) but we force to be 0
  true.covs=array(NA,dim=c(K,R,R))
  q=colSums(q.mat) ##compute the sum of the posterior weights across individuals
  pi=q/J ## average across individuals (represents relative freqnecy in data se)
  pi=normalize(pmax(0,pi))## to prevent rounding errors
  d=array(NA,dim=c(J,R,R))##create an array to store the componenet specific true covariance for each gene snp pair which will be summed
  for(k in 1:K){##loop through to compute a true covariance matrix for each Q
    if(q[k]==0){##ask about this
      #true.means[k,]=rep(0,R)
      true.covs[k,,]=array(rep(0,R*R),dim=c(R,R))}
    else{
      #true.means[k,]=1/q[k]*(q.mat[,k]*post.means[,k,])
      for(j in seq(1:J)){
        d[j,,]=q.mat[j,k]*(-post.means[j,k,]%*%t(-post.means[j,k,])+post.covs[j,k,,])##produce a RxR matrix of weighted 'truth' for each individual
      }
      true.covs[k,,]=1/q[k]*apply(d, MARGIN=c(2, 3), sum)}##get the latent identifier weighted sum of all J individuals to produce 1 RxR matrix of truth
  }
  return(list(true.covs=true.covs,pi=pi))
}


#' @title fixpoint.cov
#' @details combines e and m step
#' @param max.step a list of K pis and the true.covs arrays KxRxR which will be parsed by the e.step
#' @param dim.true.covs=kxrxr
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export
fixpoint.cov = function(max.step.unlist,b.j.hat,se.j.hat){  
  L=length(max.step.unlist)
  R=ncol(b.j.hat)
  r2=R^2
  K=L/(r2+1)
  dim.true.covs=c(K,R,R)
  pi.length=K
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  e.step=em.array.generator(max.step=max.step,b.j.hat = b.j.hat,se.j.hat = se.j.hat)
  
  post.means=e.step$post.means;
  post.covs=e.step$post.covs;
  q.mat=e.step$q.mat
  
  max.step=max.step.func(post.means=post.means,post.covs = post.covs,q.mat = q.mat)
  max.step.unlist=unlist(max.step)
  return(max.step.unlist)
}

normalize = function(x){return(x/sum(x))}

lik.func.em=function(true.covs,b.mle,V.j.hat){sapply(seq(1:K),function(k){dmvnorm(x=b.mle, sigma=true.covs[k,,] + V.j.hat)})}

negpenlogliksarah = function(max.step.unlist,b.j.hat,se.j.hat){return(-penlogliksarah(max.step.unlist,b.j.hat,se.j.hat))}

penlogliksarah = function(max.step.unlist,b.j.hat,se.j.hat){
  
  
  L=length(max.step.unlist)
  R=ncol(b.j.hat)
  r2=R^2
  K=L/(r2+1)
  dim.true.covs=c(K,R,R)

  pi.length=K
  
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  pi=max.step$pi
  true.covs=max.step$true.covs
  matrix_lik=t(sapply(seq(1:J),function(x){lik.func.em(true.covs,b.mle=b.j.hat[x,],V.j.hat=diag(se.j.hat[x,])^2)}))
  pi = (normalize(pmax(0,pi)))
  m  = pi*matrix_lik # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
    return(loglik)
}

