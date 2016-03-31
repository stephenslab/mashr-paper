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
#' @export dim.true.cov.fun


dim.true.cov.fun=function(max.step.unlist){
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

init.covmat=function(t.stat=t.stat,factor.mat=factors,lambda.mat=lambda,K,P){
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
  if(K>1){
  sfa.prox=(t(full.rank)%*%full.rank)/(M)
  true.covs[2,,]=sfa.prox
  
  svd.X=svd(X.c)
  v=svd.X$v;u=svd.X$u;d=svd.X$d
  cov.pc=1/M*v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P])%*%t(v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P]))
  true.covs[3,,]=cov.pc}
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
    lik=lik.func.em(true.covs,b.mle,V.j.hat,K)
    tinv=lapply(seq(1:K),function(k){solve(true.covs[k,,]+V.j.hat)})##create a K list of intverted Ts for each J, covariance matrix of the marginal distribution
    
    
    ##compute a K dimensional list of posterior covariance for each J, this is faster than for looping over the K, but in the main function we have to do within R also
    B.j.=lapply(seq(1:K),function(k){
      post.b.jk.ed.cov(tinv=tinv[[k]],true.covs[k,,])
    }
    )##create a K dimensional list of covariance matrices 
    post.covs[j,,,] <- tarray(array(unlist(B.j.), c( R, R,K))) ###store a K dimensional list of posterior covariances for each J (JxKxRxR) in the post.means array
    ##this is different than in my maincode, where I forloop through each K, because there i need to compute the tail probabilities for each r within a k and so I would have to forloop inside K which is ugly
    
    ##compute a K dimensional list of posterior means for each J, again this is faster than for looping over the K, but in the main function we have to do within R also
   b.j.=(lapply(seq(1:K),function(k){
      
      b=post.b.jk.ed.mean(b.mle,tinv=tinv[[k]],true.covs[k,,])##for each component, compute posterior mean
      
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
#' @export max.step.func 
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

lik.func.em=function(true.covs,b.mle,V.j.hat,K){sapply(seq(1:K),function(k){dmvnorm(x=b.mle, sigma=true.covs[k,,] + V.j.hat)})}

negpenlogliksarah = function(max.step.unlist,b.j.hat,se.j.hat){return(-penlogliksarah(max.step.unlist,b.j.hat,se.j.hat))}

penlogliksarah = function(max.step.unlist,b.j.hat,se.j.hat){
  
  
  L=length(max.step.unlist)
  R=ncol(b.j.hat)
  r2=R^2
  K=L/(r2+1)
  dim.true.covs=c(K,R,R)
  J=nrow(b.j.hat)
  pi.length=K
  
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  pi=max.step$pi
  true.covs=max.step$true.covs
  matrix_lik=t(sapply(seq(1:J),function(x){lik.func.em(true.covs,b.mle=b.j.hat[x,],V.j.hat=diag(se.j.hat[x,])^2,K)})) ##correct problem if K=1
  if(K==1){
    matrix_lik=t(matrix_lik)
  }
  pi = (normalize(pmax(0,pi)))
  m  = t(pi * t(matrix_lik))# matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
    return(loglik)
}


#' @title deconvolution.em
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

deconvolution.em <- function(t.stat,factor.mat,lambda.mat,K,P,permsnp=1000){
  init.cov=init.covmat(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K=K,P=P)
  pi=rep(1/K,K)
  R=ncol(t.stat)
  
  par.init=list(true.covs=init.cov,pi=rep(1/K,K))
  par.init.unlist=unlist(par.init)
  
  t.stat=data.frame(t.stat)
  maxes=apply(t.stat,1,function(x){mean(abs(x))})##takes the strongest t statistics
  a=cbind(t.stat,maxes)
  f=ncol(a)
  t=a[order(a$maxes,decreasing=TRUE),-f]
  t.strong=t[1:permsnp,]
  v.strong=matrix(rep(1,R*nrow(t.strong)),nrow=nrow(t.strong))

  s=squarem(par=par.init.unlist,b.j.hat=t.strong,se.j.hat=v.strong,fixptfn=fixpoint.cov, objfn=negpenlogliksarah)
  max.step.unlist=s$par
  dim.true.covs=c(K,R,R)
  pi.length=length(pi)
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  return(max.step)
}


#' @title get.prior.covar.with.max.step
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

get.prior.covar.with.max.step <- function(P, X.c,max.step,lambda.mat, Q, factor.mat,omega.table,bma=TRUE)  {
  test=list()
  R=ncol(X.c)
  M=nrow(X.c)
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=omega*diag(1,R)
    data.prox=max.step$true.covs[1,,]
    d.norm=data.prox/max(diag(data.prox))
    
    
    test[[l]][[2]]=omega*d.norm
    
    
    svd.X=svd(X.c)
    
    v=svd.X$v;u=svd.X$u;d=svd.X$d
    
    cov.pc=1/M*v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P])%*%t(v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P]))
    
    
    cov.pc.norm=cov.pc/max(diag(cov.pc))
    
    
    
    test[[l]][[3]]=omega*(cov.pc.norm)
    if(Q!=0){for(q in 1:Q){
      
      load=as.matrix(lambda.mat[,q])
      fact=as.matrix(factor.mat[q,])
      rank.prox=load%*%t(fact)
      a=(1/M*(t(rank.prox)%*% rank.prox))
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[l]][[q+3]]=omega*a.norm
    }}
    full.rank=as.matrix(lambda.mat)%*%as.matrix(factor.mat)
    b=(1/M*(t(full.rank)%*%full.rank))
    b[is.nan(b)]=0
    b.norm=b/max(diag(b))
    test[[l]][[Q+4]]=omega*b.norm
    
    if(bma==TRUE){
      configs=matrix(0,nrow=R,ncol=R)
      
      R=ncol(factor.mat)
      for(r in 1:R){
        configs[r,r]=1}
      
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        
        mat=(configs[c,]%*%t(configs[c,]))
        
        test[[l]][[Q+4+c]]=omega*mat}}}
  return(U.0kl=test)
}


#' @title compute.hm.covmat
#' @details use a prespecified X`X from the output of the 'deconvoluting EM'
#' @param max.step: the list ouput of deconvolution EM (list of length 2, contains the denoised covariance  matrices and vector of pi)
#' @param Q number of single rank factors to include in the set of covariance matrices
#' @param P PC approximation
#' @param BMA Whether or not to include singleton and full configurations
#' @return a list of covariance matrices
 
compute.hm.covmat = function(t.stat,v.j,Q,lambda.mat,P,A,factor.mat,max.step){
X.real=as.matrix(t.stat)
X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
R=ncol(X.c)
omega=mult.tissue.grid(mult=sqrt(2),t.stat,v.j)
omega.table=data.frame(omega)
lambda.mat=lambda.mat
A=A
factor.mat=factor.mat
U.0kl=get.prior.covar.with.max.step(P = P,X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
covmat=unlist(U.0kl,recursive=F)
saveRDS(covmat,paste0("covmat",A,".rds"))

return(covmat)}

#' @title test.funct
#' @details shows that the e step properly stores the posterior mean for a given gene snp pair and componenet combination
#' @param j gene snp pair of interest
#' @param k componenet of interest
#' @param R number of tissues

test.funct=function(j,max.step.unlist,k,R){
  
  (L=length(max.step.unlist))
  K=L/(R^2+1)
  dim.true.covs=c(K,R,R)
  pi.length=K
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  true.covs=max.step$true.covs
  pi=max.step$pi
  
  b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
  V.j.hat=diag(se.j.hat[j,]^2)
  tinv=lapply(seq(1:K),function(k){solve(true.covs[k,,]+V.j.hat)})##create a K list of intverted Ts for each J
  lik=sapply(seq(1:K),function(k){dmvnorm(x=b.mle, sigma=true.covs[k,,] + V.j.hat)})##compute K element likeilihood for each idndiviual
  B.j.=(lapply(seq(1:K),function(k){
    #tinv=solve(true.covs[k,,]+V.j.hat)##covariance matrix of the marginal distribution
    post.b.jk.ed.cov(tinv=tinv[[k]],true.covs[k,,])
  }
  ))##create a K dimensional list of covariance matrices 
  
  ##compute a K dimensional list of posterior means for each J
  b.j.=(lapply(seq(1:K),function(k){
    #tinv=solve(true.covs[k,,]+V.j.hat)##covariance matrix of the marginal distribution
    post.b.jk.ed.mean(b.mle,tinv=tinv[[k]],true.covs[k,,])##for each component, compute posterior mean
  }))
  
  
  test=em.array.generator(max.step = max.step,b.j.hat = b.j.hat,se.j.hat = se.j.hat)### show that the output of the posterior means in the E step matches the actual computation above
  pm=test$post.means;pc=test$post.covs;q.mat=test$q.mat
  
  par(mfrow=c(1,2))
  plot(pm[j,k,],b.j.[[k]])##test to make sure posterior mean is stored properly
  plot(diag(pc[j,k,,]),diag(B.j.[[k]]))##test to make sure posterior covariance is stored properly
}




#' @title get.prior.covar.with.all.max.step
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

get.prior.covar.with.all.max.step <- function(X.c,max.step,lambda.mat, Q, factor.mat,omega.table,bma=TRUE)  {
  test=list()
  R=ncol(X.c)
  M=nrow(X.c)
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=omega*diag(1,R)
    data.prox=max.step$true.covs[1,,]
    d.norm=data.prox/max(diag(data.prox))
    
    
    test[[l]][[2]]=omega*d.norm
    
  
    
    cov.pc=max.step$true.covs[3,,]
    
    
    cov.pc.norm=cov.pc/max(diag(cov.pc))
    
    
    
    test[[l]][[3]]=omega*(cov.pc.norm)
    if(Q!=0){for(q in 1:Q){
      
      load=as.matrix(lambda.mat[,q])
      fact=as.matrix(factor.mat[q,])
      rank.prox=load%*%t(fact)
      a=(1/M*(t(rank.prox)%*% rank.prox))
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[l]][[q+3]]=omega*a.norm
    }}
    full.rank=max.step$true.covs[2,,]
    b.norm=full.rank/max(diag(full.rank))
    test[[l]][[Q+4]]=omega*b.norm
    
    if(bma==TRUE){
      configs=matrix(0,nrow=R,ncol=R)
      
      R=ncol(factor.mat)
      for(r in 1:R){
        configs[r,r]=1}
      
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        
        mat=(configs[c,]%*%t(configs[c,]))
        
        test[[l]][[Q+4+c]]=omega*mat}}}
  return(U.0kl=test)
}

  
#' @title compute.hm.covmat.all.max.step
#' @details use a prespecified X`X from the output of the 'deconvoluting EM'
#' @param max.step: the list ouput of deconvolution EM (list of length 2, contains the denoised covariance  matrices and vector of pi)
#' @param Q number of single rank factors to include in the set of covariance matrices
#' @param P PC approximation
#' @param b.hat,se.hat matrices to help choose grid weights (JxR)
#' @param v.j a JxR matrix of standard errors for scaling (in the zstat case this is all 1s)
#' @param BMA Whether or not to include singleton and full configurations
#' @return a list of covariance matrices

compute.hm.covmat.all.max.step = function(b.hat,se.hat,t.stat,v.j,Q,lambda.mat,A,factor.mat,max.step,zero=FALSE){
  X.real=as.matrix(t.stat)
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  R=ncol(X.c)
  omega=mult.tissue.grid(mult=sqrt(2),b.hat,se.hat)
  omega.table=data.frame(omega)
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  U.0kl=get.prior.covar.with.all.max.step(X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
  covmat=unlist(U.0kl,recursive=F)
  if(zero==TRUE){
    z=matrix(rep(0,R*R),ncol=R,nrow=R)
    covmat=c(covmat,list(z))
  }
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}
  

#' @title get.prior.covar.with.all.max.step
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

get.prior.covar.with.all.max.step.no.pc <- function(X.c,max.step,lambda.mat, Q, factor.mat,omega.table,bma=TRUE)  {
  test=list()
  R=ncol(X.c)
  M=nrow(X.c)
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=omega*diag(1,R)
    data.prox=max.step$true.covs[1,,]
    d.norm=data.prox/max(diag(data.prox))
    
    
    test[[l]][[2]]=omega*d.norm
    
    
    

    if(Q!=0){for(q in 1:Q){
      
      load=as.matrix(lambda.mat[,q])
      fact=as.matrix(factor.mat[q,])
      rank.prox=load%*%t(fact)
      a=(1/M*(t(rank.prox)%*% rank.prox))
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[l]][[q+2]]=omega*a.norm
    }}
    full.rank=max.step$true.covs[2,,]
    b.norm=full.rank/max(diag(full.rank))
    test[[l]][[Q+3]]=omega*b.norm
    
    if(bma==TRUE){
      configs=matrix(0,nrow=R,ncol=R)
      
      R=ncol(factor.mat)
      for(r in 1:R){
        configs[r,r]=1}
      
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        
        mat=(configs[c,]%*%t(configs[c,]))
        
        test[[l]][[Q+3+c]]=omega*mat}}}
  return(U.0kl=test)
}



compute.hm.covmat.all.max.step.no.pc = function(t.stat,v.j,Q,lambda.mat,A,factor.mat,max.step){
  X.real=as.matrix(t.stat)
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  R=ncol(X.c)
  omega=mult.tissue.grid(mult=sqrt(2),t.stat,v.j)
  omega.table=data.frame(omega)
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  U.0kl=get.prior.covar.with.all.max.step.no.pc(X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
  covmat=unlist(U.0kl,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}

#' @title deconvolution.em.with.max.iter
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

deconvolution.em.with.max.iter <- function(t.stat,factor.mat,lambda.mat,K,P,permsnp=1000,maxiter=1000){
  init.cov=init.covmat(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K=K,P=P)
  pi=rep(1/K,K)
  R=ncol(t.stat)
  
  par.init=list(true.covs=init.cov,pi=rep(1/K,K))
  par.init.unlist=unlist(par.init)
  
  t.stat=data.frame(t.stat)
  maxes=apply(t.stat,1,function(x){mean(abs(x))})##takes the strongest t statistics
  a=cbind(t.stat,maxes)
  f=ncol(a)
  t=a[order(a$maxes,decreasing=TRUE),-f]
  t.strong=t[1:permsnp,]
  v.strong=matrix(rep(1,R*nrow(t.strong)),nrow=nrow(t.strong))
  maxiter=maxiter
  s=squarem(par=par.init.unlist,b.j.hat=t.strong,se.j.hat=v.strong,fixptfn=fixpoint.cov, control=list(maxiter=maxiter))#objfn=negpenlogliksarah)
  max.step.unlist=s$par
  dim.true.covs=c(K,R,R)
  pi.length=length(pi)
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  return(max.step)
}


#' @title deconvolution.em.with.loop
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

deconvolution.em.with.loop <- function(t.stat,factor.mat,lambda.mat,K,P,permsnp=1000,maxiter=1000){
  init.cov=init.covmat(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K=K,P=P)
  pi=rep(1/K,K)
  R=ncol(t.stat)
  
  par.init=list(true.covs=init.cov,pi=rep(1/K,K))
  par.init.unlist=unlist(par.init)
  
  t.stat=data.frame(t.stat)
  maxes=apply(t.stat,1,function(x){mean(abs(x))})##takes the strongest t statistics
  a=cbind(t.stat,maxes)
  f=ncol(a)
  t=a[order(a$maxes,decreasing=TRUE),-f]
  t.strong=t[1:permsnp,]
  v.strong=matrix(rep(1,R*nrow(t.strong)),nrow=nrow(t.strong))
  maxiter=maxiter
  max.step.unlist=par.init.unlist
  for(i in 1:(maxiter-1)){
    m=fixpoint.cov(max.step.unlist,b.j.hat = t.strong,se.j.hat = v.strong)
    max.step.unlist=m
  }
 
  nb=negpenlogliksarah(max.step.unlist = max.step.unlist,b.j.hat = t.strong,se.j.hat = v.strong)
  m=fixpoint.cov(max.step.unlist,b.j.hat = t.strong,se.j.hat = v.strong)
  max.step.unlist=m
  negpen=negpenlogliksarah(max.step.unlist = max.step.unlist,b.j.hat = t.strong,se.j.hat = v.strong)-nb
  
  dim.true.covs=c(K,R,R)
  pi.length=length(pi)
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  return(list(max.step=max.step,negpen=negpen))
}



#' @title deconvolution.em.bovy
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat
#' @param P = rank of PC approxiatmion
#' @param Q = rank of SFA approximation
#' @param permsnp = number of strong stats you want to train on (default is 1000)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

deconvolution.em.with.bovy=function(t.stat,factor.mat,v.j,lambda.mat,K,P){
R=ncol(t.stat)
init.cov=init.covmat(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K=K,P=P)
init.cov.list=list()
for(i in 1:K){init.cov.list[[i]]=init.cov[i,,]}
head(init.cov.list)
mean.mat=matrix(rep(0,ncol(t.stat)*nrow(t.stat)),ncol=ncol(t.stat),nrow=nrow(t.stat))

ydata=  t.stat
xamp= rep(1/K,K)
xcovar= init.cov.list
fixmean= TRUE     
ycovar=  v.j     
xmean=   mean.mat   
projection= list();for(l in 1:nrow(t.stat)){projection[[l]]=diag(1,R)}

e=extreme_deconvolution(ydata=ydata,ycovar=ycovar,xamp=xamp,xmean=xmean,xcovar=init.cov.list,fixmean=T,projection=projection)

true.covs=array(dim=c(K,R,R))
for(i in 1:K){true.covs[i,,]=e$xcovar[[i]]}
pi=e$xamp
max.step=list(true.covs=true.covs,pi=pi)
return(max.step)}



init.covmat.single=function(t.stat=t.stat,factor.mat=factors,lambda.mat=lambda,K,P,Q){
  K=K
  R=ncol(t.stat)#number of tissues
  true.covs=array(NA,dim=c(K+Q,R,R))
  
  X.t=as.matrix(t.stat)
  X.real=X.t
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  
  M=nrow(X.c)
  data.prox=((t(X.c)%*%X.c)/M)
  true.covs[1,,]=data.prox
  full.rank=lambda.mat%*%factor.mat
  if(K>1){
  sfa.prox=(t(full.rank)%*%full.rank)/(M)
  true.covs[2,,]=sfa.prox
  
  svd.X=svd(X.c)
  v=svd.X$v;u=svd.X$u;d=svd.X$d
  cov.pc=1/M*v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P])%*%t(v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P]))
  true.covs[3,,]=cov.pc
  
  
  
   if(Q!=0){for(q in 1:Q){
      
      load=as.matrix(lambda.mat[,q])
      fact=as.matrix(factor.mat[q,])
      rank.prox=load%*%t(fact)
      a=(1/M*(t(rank.prox)%*% rank.prox))
      
    true.covs[3+q,,]=a
      }}}
  return(true.covs)
}






compute.hm.covmat.all.max.step.withQ = function(b.hat,se.hat,t.stat,v.j,lambda.mat,A,factor.mat,max.step,Q){
  X.real=as.matrix(t.stat)
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  R=ncol(X.c)
  omega=mult.tissue.grid(mult=sqrt(2),b.hat,se.hat)
  omega.table=data.frame(omega)
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  U.0kl=get.prior.covar.with.all.max.step.withQ(X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
  covmat=unlist(U.0kl,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}


get.prior.covar.with.all.max.step.withQ <- function(X.c,max.step,lambda.mat, Q, factor.mat,omega.table,bma=TRUE)  {
  test=list()
  R=ncol(X.c)
  M=nrow(X.c)
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=omega*diag(1,R)
    data.prox=max.step$true.covs[1,,]
    d.norm=data.prox/max(diag(data.prox))
    
    
    test[[l]][[2]]=omega*d.norm
    
    
    
    cov.pc=max.step$true.covs[3,,]
    
    
    cov.pc.norm=cov.pc/max(diag(cov.pc))
    
    
    
    test[[l]][[3]]=omega*(cov.pc.norm)
    if(Q!=0){for(q in 1:Q){
      
      
      rank.prox=max.step$true.covs[3+q,,]
      a=(1/M*(t(rank.prox)%*% rank.prox))
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[l]][[q+3]]=omega*a.norm
    }}
    full.rank=max.step$true.covs[2,,]
    b.norm=full.rank/max(diag(full.rank))
    test[[l]][[Q+4]]=omega*b.norm
    
    if(bma==TRUE){
      configs=matrix(0,nrow=R,ncol=R)
      
      R=ncol(factor.mat)
      for(r in 1:R){
        configs[r,r]=1}
      
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        
        mat=(configs[c,]%*%t(configs[c,]))
        
        test[[l]][[Q+4+c]]=omega*mat}}}
  return(U.0kl=test)
}


compute.hm.covmat.all.max.step.squared = function(b.hat,se.hat,t.stat,v.j,Q,lambda.mat,A,factor.mat,max.step){
  X.real=as.matrix(t.stat)
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  R=ncol(X.c)
  omega=mult.tissue.grid(mult=sqrt(2),b.hat,se.hat)
  omega.table=data.frame(omega)
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  U.0kl=get.prior.covar.with.all.max.step(X.c,max.step = max.step,lambda.mat = lambda.mat,Q = Q,factor.mat = factor.mat,omega.table=omega.table,bma = TRUE)
  covmat=unlist(U.0kl,recursive=F)
  covmat=lapply(covmat,function(x){x^2})
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}