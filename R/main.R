#' @title post.b.gpkl.cov
#' @param U.0k.l let U.0k.l represent a specific matrix in U.0kl (.e.g, U.0kl[[l]][[k]])
#' @return post.b.gpkl.cov ## returns an R*R posterior covariance matrix for jth gene snp pair
#' @return post.b.gpkl.mean return a 1 * R vector of posterior means for a given prior covariance matrix
#' @export
post.b.gpkl.cov <- function(V.gp.hat.inv, U.0k.l){
    U.gp1kl <- U.0k.l %*% solve(V.gp.hat.inv %*% U.0k.l + diag(nrow(U.0k.l)))
    return(U.gp1kl)
}


#' @title post.b.gpkl.mean
#' @param U.0k.l let U.0k.l represent a specific matrix in U.0kl (.e.g, U.0kl[[l]][[k]])
#' @return post.b.gpkl.mean return a 1 * R vector of posterior means for a given prior covariance matrix
#' @export
post.b.gpkl.mean <- function(b.mle, V.gp.hat.inv, U.gp1kl){
    mu.gp1kl <- U.gp1kl %*% V.gp.hat.inv %*% b.mle
    return(mu.gp1kl)
}

#' @title get.prior.covar.Ukl.with.rho 
#' @description Compute the prior covariance matrix for each gene SNP pair in the rows and componenets in the columns
#' @param P integer number of PCs
#' @param omega.table table of grid weights
#' @param rhos correlation factor for morphed matrices
#' @param X.c = column centered matrix of strong t statistics from which to derive covariance matrices
#' @return U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export
get.prior.covar.Ukl.with.rho <- function(P, lambda.mat,X.c, Q, factor.mat,omega.table,rhos)  {
    test=list()
    id=diag(1,R)
    M=nrow(X.c)
    for(l in 1:nrow(omega.table)){
        test[[l]]=list()
        omega=omega.table[l,]
        for(r in 1:length(rhos)){
          test[[l]][[r]]=list()
          rho=rhos[r]
        
        test[[l]][[r]][[1]]=rho*omega*id+(1-rho)*omega*id
        data.prox=((t(X.c)%*%X.c)/M)
        d.norm=data.prox/max(diag(data.prox))
        
        
        test[[l]][[r]][[2]]=rho*omega*id+ (1-rho)*omega*d.norm
        
        
        svd.X=svd(X.c)
        
        v=svd.X$v;u=svd.X$u;d=svd.X$d
        
        cov.pc=1/M*v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P])%*%t(v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P]))
        
        
        cov.pc.norm=cov.pc/max(diag(cov.pc))
        
        
        
        test[[l]][[r]][[3]]=rho*id+(1-rho)*omega*(cov.pc.norm)
        if(Q!=0){for(q in 1:Q){
            fact=factor.mat
            load=as.matrix(lambda[,q])
            fact=t((as.matrix(fact[q,])))
            rank.prox=load%*%fact
            a=(1/M*(t(rank.prox)%*% rank.prox))
            a[is.nan(a)] = 0
            a.norm=a/max(diag(a))
            test[[l]][[r]][[q+3]]=rho*id+(1-rho)*omega*a.norm
        }
        full.rank=as.matrix(lambda)%*%as.matrix(factor.mat)
        b=(1/M*(t(full.rank)%*%full.rank))
        b[is.nan(b)]=0
        b.norm=b/max(diag(b))
            test[[l]][[r]][[Q+4]]=rho*id+(1-rho)*omega*b.norm}}}
    return(U.0kl=test)
}


#' @title get.prior.covar.Ukl
#' @param P number of pcs (must be greater than 1)
#' @param lambda.mat matrix of JxQ lambdas
#' @param Q number of single rank factor matrices (if 0, then only the rank 5 approximation wil be used)
#' @param factor.mat nfactors x R tissue matrix of factors
#' @param omega.table Lx2 table of grid weights
#' @param bma whether or not to use singleton and shared configurations
#' @param X.c = column centered matrix of strong t statistics from which to derive covariance matrices
#' @return U.0kl
#' @export
get.prior.covar.Ukl <- function(P, X.c,lambda.mat, Q, factor.mat,omega.table,bma=TRUE)  {
  test=list()
  M=nrow(X.c)
  R=ncol(X.c)
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=omega*diag(1,R)
    data.prox=((t(X.c)%*%X.c)/M)
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

#' @title get.prior.covar.bma.only
#' @param R number of tissues
#' @param omega.table L vector grid weights
#' @return a L x K list of covariance matrices
#' @export
get.prior.covar.bma.only <- function(R,omega.table)  {
  test=list()
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
  
    configs=matrix(0,nrow=R,ncol=R)
    for(r in 1:R){
        configs[r,r]=1}
      configs=rbind(configs,rep(1,R))
      for(c in 1:nrow(configs)) {
        mat=(configs[c,]%*%t(configs[c,]))
      test[[l]][[c]]=omega*mat}}
  return(U.0kl=test)
}

#' @title autoselect.mix.sd 
#' @param betahat is an J x R matrix of betahats (or t statistics)
#' @param sebetahat is a J x R matrix of their standard errors (or 1s, if T staitsitcs)
#' @return mult^((-npoint):0) * sigmaamax)
#' @export


autoselect.mix.sd = function(betahat,sebetahat,mult){
    sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
    sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
    if(all(betahat^2<=sebetahat^2)){
        sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
    }else{
        sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
    }
    if(mult==0){
        return(c(0,sigmaamax/2))
    }else{
        npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
        return(mult^((-npoint):0) * sigmaamax)
    }
}


#' @title mult.tissue.grid 
#' @param mult is density
#' @param betahat is matrix of betahats
#' @param sebetahat is matrix of sebetahats
#' @return omega
#' @export
mult.tissue.grid = function(mult=sqrt(2),betahat,sebetahat,maxp,minp){ 
  R=ncol(betahat);
mix.weights=unlist(sapply(seq(1:ncol(betahat)),function(r){
  autoselect.mix.sd(betahat = betahat[,r],sebetahat = sebetahat[,r],mult=2)}))
    ;sigmaamin=minp*min(mix.weights);sigmaamax=maxp*max(mix.weights);npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult));
    omega=mult^((-npoint):0) * sigmaamax;return(omega)}
##we tried with omega^2 as output, to be consistent with production of covariance matrices (not coSD matrices) but better likelihood with denser narrower grid

#' @title lik.func 
#' @details computes likelihood for each betahat
#' @param b.mle Rx1 vector of mles
#' @param V.gp.hat RxR matrix of standard errors
#' @param U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export

lik.func=function(b.mle,V.gp.hat,covmat)
{ sapply(seq(1:length(covmat)),function(x){dmvnorm(x=b.mle, sigma=covmat[[x]] + V.gp.hat)})
}


#' @title log.lik.func 
#' @details computes log.likelihood for each betahat
#' @param b.mle Rx1 vector of mles
#' @param V.gp.hat RxR matrix of standard errors
#' @param U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export log.lik.func

log.lik.func=function(b.mle,V.gp.hat,covmat)
{ sapply(seq(1:length(covmat)),function(x){dmvnorm(x=b.mle, sigma=covmat[[x]] + V.gp.hat,log=TRUE)})
}

#' @title hm.weight.gen
#' @param lik.mat = JxK matrix of marginal likelihoods
#' @param Jtrain = size of train set
#' @param Jtest=size of test set
#' @return compute HM weights on training and test set for compoarsion
#' @export
hm.weight.gen=function(lik.mat,jtrain,jtest){
    
    train=lik.mat[1:jtrain,]
    test=lik.mat[jtrain:nrow(lik.mat),]
    pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
    return(pis$pihat)
}


#' @title total.lik.func
#' @param test = J x R matrix of likielihoods from test data set
#' @param pis HM weights
#' @return likelihood of whole data set
#' @export

total.lik.func=function(test,pis){
    
    sum(log(test%*%pis))}



#' @title post.weight.func 
#' @details converts the matrix of likelihood for each gene snp pairs to matrix of posterior weights 
#' @param pis = object from EM output with prior weight P(Z=K) as computed from
#' @param lik.mat = a JxK matrix of likelihoods (may be training set) for the P(D|Z=K)
#' @return a 1xK vector of posterior weight for each gene snp pait
#' @export

post.weight.func=function(pis,lik.mat){d=t(apply(lik.mat,1,function(x){x*pis}))
    
    marg=rowSums(d)
    return(d/marg)}


#' @title post.array.generator
#' @param b.mle =  Rx1 vector of beta.hats
#' @param V.gp.hat = Rx1 vector of standard errors
#' @param covmat = LxK dimenstional (unlisted list) of prior covariance matrices
#' @return post.means JxKxR array of posterior means, correspodning to the posterior mean ##' for the Jth individual in the Kth compoenent across all R tissues
#' @return post.covs JxKxR array of posterior vars, correspodning to the posterior vars ##' for the Jth individual in the Kth compoenent across all R tissues
#' @return post.ups JxKxR array of posterior tailprobs, corresponding to the marginal
#' @return upper tail probability for the Jth individual in the Kth compoenent across all R
#' @return post.ups JxKxR array of posterior tailprobs, corresponding to the marginal
#' @return upper  tail probability for the Jth individual in the Kth component across all R
#' @return post.nulls JxKxR array of posterior nullprobs, corresponding to the marginal
#' @return "null probability"" for the Jth individual in the Kth component across all R
#' @export

post.array.generator=function(b.gp.hat,J,A,se.gp.hat,covmat){
    
    R=ncol(b.gp.hat)
    
    K=length(covmat)
    post.means=array(NA,dim=c(J,K,R))
    post.covs=array(NA,dim=c(J,K,R))
    post.ups=array(NA,dim=c(J,K,R))
    post.nulls=array(NA,dim=c(J,K,R))
    post.downs=array(NA,dim=c(J,K,R))
    
    
    
    for(j in 1:J){
        b.mle=as.vector(t(b.gp.hat[j,]))##turn i into a R x 1 vector
        V.gp.hat=diag(se.gp.hat[j,])^2
        V.gp.hat.inv <- solve(V.gp.hat)
        
        for(k in 1:K){
            
            
            U.gp1kl <- (post.b.gpkl.cov(V.gp.hat.inv, covmat[[k]]))
            mu.gp1kl <- as.array(post.b.gpkl.mean(b.mle, V.gp.hat.inv, U.gp1kl))
            post.means[j,k,]=mu.gp1kl
            post.covs[j,k,]=as.array(diag(U.gp1kl))
            for(r in 1:R){##` marginal calculation for each tissue in each component
                if(post.covs[j,k,r]==0){###if the covariance matrix has entry 0, then p(null=1)
                    post.ups[j,k,r]=0#need to figure out how to make the list have individual entries
                    post.downs[j,k,r]=0
                    post.nulls[j,k,r]=1}
                else{
                    post.ups[j,k,r]=pnorm(0,mean=mu.gp1kl[r],sd=sqrt(diag(U.gp1kl)[r]),lower.tail=F)
                    post.downs[j,k,r]=pnorm(0,mean=mu.gp1kl[r],sd=sqrt(diag(U.gp1kl)[r]),lower.tail=T)
                    post.nulls[j,k,r]=0}
            }
        }
    }
    
    return(list(post.means=post.means,post.nulls=post.nulls,post.covs=post.covs,post.ups=post.ups,post.downs=post.downs))
}




#' @title total.mean
#' @description generate a K x R matrix of post.weighted quantieis for each gene snp pair and sum them to get total weighted
#' @param post.means J x K x R arrays of posterior means for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#' @return get vector of total weighted quanties for each gene SNP pair
#' @export

total.mean=function(j,post.means,post.weights){weightedmeans=(as.matrix(post.means[j,,1:R]*post.weights[j,]))
    colSums(weightedmeans)}

#' @title total.up
#' @description generate a K x R matrix of post.weighted upper tail probabilites and sum them to get total weighted
#' @param post.ups J x K x R arrays of posterior upper tail probabilities for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#' @return get vector of total weighted quanties for each gene SNP pair
#' @export
total.up=function(j,post.ups,post.weights){colSums(as.matrix(post.ups[j,,1:R]*post.weights[j,]))}

#' @title total.down
#' @description generate a K x R matrix of post.weighted lower tail probabilites and sum them to get total weighted
#' @param post.downs J x K x R arrays of posterior lower tail probabilities for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#' @export

total.down=function(j,post.downs,post.weights){colSums(as.matrix(post.downs[j,,1:R]*post.weights[j,]))}

#' @title total.null
#' @description generate a K x R matrix of post.weighted null probabilites and sum them to get total weighted
#' @param post.nulls J x K x R arrays of posterior upper tail probabilities for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#' @export

total.null=function(j,post.nulls,post.weights){colSums(as.matrix(post.nulls[j,,1:R]*post.weights[j,]))}

#' @title total.covs.partone
#' @description generate a K x R matrix of post.weighted marginal variances and sum them to get total weighted
#' @param post.covs J x K x R arrays of posterior marginal variance (diagonals of posteior covmat) for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#' @export

total.covs.partone=function(j,post.means,post.covs,post.weights){colSums(as.matrix((post.covs[j,,1:R]+post.means[j,,1:R]^2)*post.weights[j,]))}




#' @title compare.func
#' @details compares upper and lower tail probability at each tissue to determine larger
#' @param all.upper and all.lower are JxR matrices of weighted upper and lower tail probabilities for all gene pairs
#' @return 1 x R vector of lfsr at each tissue
#' @export
compare.func=function(j,all.upper,all.lower){
    as.matrix(apply(rbind(all.upper[j,],all.lower[j,]),2,function(j){1-max(j)}))}







#' @title compute.covmat.with.rho
#' @export

compute.covmat.with.rho = function(b.gp.hat,sebetahat,Q,lambda.mat,P,A,factor.mat,t.stat){
  

  R=ncol(t.stat)#number of tissues
  
  
  X.t=as.matrix(t.stat)
  X.real=X.t
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  
  M=nrow(X.c)
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
    
    omega.table=data.frame(omega)
    rhos=seq(0.1,1,by=0.1)
    lambda=lambda
    A=A
    factor.mat=factor.mat
    X.c=X.c
    Q=Q
    
 
    U.0kl=get.prior.covar.Ukl.with.rho(P=2,lambda=lambda,Q=Q,factor.mat=factor.mat, X.c,omega.table=omega.table,rhos=rhos)
    
    covmat1=unlist(U.0kl,recursive=F)
    covmat=unlist(covmat1,recursive=F)
    saveRDS(covmat,paste0("covmat",A,".rds"))
    
    return(covmat)}


#' @title compute.covmat
#' @param b.gp.hat a JxR matrix of betahats - necessary for grid computation
#' @param sebetahat a JxR matrix of their standard errors - necessary for grid computation
#' @param Q number of factors to consider
#' @param t.stat matrix of strong t statistics (needs to be R in columns but not necessarily j in rows)
#' @param factor mat a fxR matrix of factors
#' @param BMA=true wheter or not to use singleton and shared config
#' @param zero=true whether or not to include zero matrix
#' @return A list of covariance matrices
#' @export

compute.covmat = function(b.gp.hat,sebetahat,Q,t.stat,lambda.mat,P,A,factor.mat,bma=TRUE,zero=FALSE,maxp=1,minp=1,power=1){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat,maxp,minp)
  
  omega.table=data.frame(omega^power)
 
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat

  Q=Q
  R=ncol(b.gp.hat)
  R=ncol(t.stat)#number of tissues
  
  X.t=as.matrix(t.stat)
  X.real=X.t
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  
  M=nrow(X.c)
  
  
  
  
  
U.0kl=get.prior.covar.Ukl(P,X.c,lambda.mat=lambda.mat,Q=Q,factor.mat=factor.mat, omega.table=omega.table,bma)
    
   
    covmat=unlist(U.0kl,recursive=F)
if(zero==TRUE){
  z=matrix(rep(0,R*R),ncol=R,nrow=R)
  covmat=c(covmat,list(z))
}
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}


#' @title compute.covmat.bma.only
#' @param b.gp.hat a JxR matrix of betahats
#' @param sebetahat a JxR matrix of their standard errors
#' @return A list of covariance matrices
#' @export


compute.covmat.bma.only = function(b.gp.hat,sebetahat,A){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
  
  omega.table=data.frame(omega)
  
    
  U.0=get.prior.covar.bma.only(R = ncol(b.gp.hat),omega.table = omega.table)
  
  
  covmat=unlist(U.0,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}




#' @title compute.mixture.dist
#' @description compute a list of JxKxR arrays with componenet specific quanitites and a JxK matrix of posterior weights
#' @param b.gp.hat JxR matrix of MLE betahats
#' @param J number of gene snp pairs to consider
#' @param se.gp.hat JxR matrix of standard errors
#' @param A file name to save things in
#' @param save=FALSE whether to save array binaries (can be large memory constraint)
#' @return all.arrays a list of JxKxR arrays of posterior quantities
#' @export
compute.mixture.dist=function(b.gp.hat,J,se.gp.hat,covmat,A,save=FALSE){
    
    J=J
    R=ncol(b.gp.hat)
    
    if(file.exists(paste0("likelihood",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=b.gp.hat[x,],V.gp.hat=diag(se.gp.hat[x,])^2,covmat)}))
    
    saveRDS(lik.mat,paste0("likelihood",A,".rds"))}
    
    else(lik.mat=readRDS(paste0("likelihood",A,".rds")))
    
    
    pis=hm.weight.gen(lik.mat,J/2,J/2)
    write.table(pis,paste0("pis",A,".txt"))
    test=lik.mat[((J/2):(J)),]
    write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
    
    post.weights=as.matrix(post.weight.func(pis,lik.mat))
    
    saveRDS(post.weights,paste0("post.weight.",A,".rds"))
    
    rm(post.weights) ## to conserve memory
    rm(lik.mat) ## to conserve memory
    
    all.arrays=post.array.generator(b.gp.hat,A=A,J=J,se.gp.hat,covmat)
    
 
    
    
    if(save==TRUE){
    
    post.means=all.arrays$post.means
    post.covs=all.arrays$post.covs
    post.ups=all.arrays$post.ups
    post.downs=all.arrays$post.downs
    post.nulls=all.arrays$post.nulls
    
    # return(list(post.means,...))
    saveRDS(post.means,paste0("post.means",A,".rds"))
    saveRDS(post.covs,paste0("post.covs",A,".rds"))
    saveRDS(post.ups,paste0("post.ups",A,".rds"))
    saveRDS(post.downs,paste0("post.downs",A,".rds"))
    saveRDS(post.nulls,paste0("post.nulls",A,".rds"))}
    else(return(all.arrays))
  }

#'@title compute.total.quant
#'@param all.arrays list of JxKxR arrays of posterior componenet-specific quantities
#'@param J number of gene snp pairs to consider
#'@param A filename to save weighted quantities in
#'@param If file present, read array file; else use list in memory
#'@export
compute.total.quant=function(A,J,all.arrays){
   
      if(missing(all.arrays)){
    post.means=readRDS(file=paste0("post.means",A,".rds"))
    post.covs=readRDS(file=paste0("post.covs",A,".rds"))
    post.ups=readRDS(file=paste0("post.ups",A,".rds"))
    post.nulls=readRDS(file=paste0("post.nulls",A,".rds"))
    post.downs=readRDS(file=paste0("post.downs",A,".rds"))
    post.weights=readRDS(file=paste0("post.weight.",A,".rds"))
    }
    
     
  
  else{
      post.means=all.arrays$post.means
      post.covs=all.arrays$post.covs
      post.ups=all.arrays$post.ups
      post.downs=all.arrays$post.downs
      post.nulls=all.arrays$post.nulls
      post.weights=readRDS(file=paste0("post.weight.",A,".rds"))
  }
    
    all.covs.partone=t(sapply(seq(1:J),function(j){total.covs.partone(j,post.means,post.covs,post.weights)}))

    all.mus=t(sapply(seq(1:J),function(j){total.mean(j,post.means,post.weights)}))
    all.upper=t(sapply(seq(1:J),function(j){total.up(j,post.ups,post.weights)}))
    all.lower=t(sapply(seq(1:J),function(j){total.down(j,post.downs,post.weights)}))
    all.nuller=t(sapply(seq(1:J),function(j){total.null(j,post.nulls,post.weights)}))
    
    lfsr.mat=t(sapply(seq(1:J),function(j){compare.func(j,all.upper,all.lower)}))
    marginal.var=all.covs.partone-all.mus^2
    
    
    write.table(marginal.var,paste0("marginal.var.",A,".txt"))
    write.table(all.mus,paste0("post.mean.",A,".txt"))
    write.table(all.upper,paste0("post.up.",A,".txt"))
    write.table(all.lower,paste0("post.low.",A,".txt"))
    write.table(all.nuller,paste0("post.null.",A,".txt"))
    write.table(lfsr.mat,paste0("lfsr.",A,".txt"))
}




#'@title checkfunc
#'@export

checkfunc=function(j,b.gp.hat,se.gp.hat,A,k ) {
    post.means=readRDS(file=paste0("post.means",A,".rds"))
    post.covs=readRDS(file=paste0("post.covs",A,".rds"))
    post.ups=readRDS(file=paste0("post.ups",A,".rds"))
    post.nulls=readRDS(file=paste0("post.nulls",A,".rds"))
    post.downs=readRDS(file=paste0("post.downs",A,".rds"))
    posterior.means=read.table(file=paste0("post.mean.",A,".txt"))
    post.weights=readRDS(paste0("post.weight.",A,".rds"))
    
    b.mle=t(as.vector(b.gp.hat[j,]))
    V.gp.hat=diag(se.gp.hat[j,])^2
    V.gp.hat.inv <- solve(V.gp.hat)
    
    U.gp1kl <- post.b.gpkl.cov(V.gp.hat.inv, covmat[[k]])
    
    pdf(paste0("matchingvariance",A,".pdf"))
    plot(diag(U.gp1kl),post.covs[j,k,1:R])
    dev.off()
    
    pdf(paste0("matchingmus",A,".pdf"))
    mu.gp1kl <- as.array(post.b.gpkl.mean(b.mle, V.gp.hat.inv, U.gp1kl))
    plot(mu.gp1kl,post.means[j,k,])
    dev.off()
    
    pdf(paste0("postmeancheck",A,".pdf"))
    plot((post.weights[j,]%*%post.means[j,,]),posterior.means[j,])
    dev.off()
    
}

#'@title post.array.per.snp
#'@export

post.array.per.snp=function(j,covmat,b.gp.hat,se.gp.hat){

  
  R=ncol(b.gp.hat)
  K=length(covmat)
  post.means=array(NA,dim=c(K,R))
  post.covs=array(NA,dim=c(K,R))
  post.ups=array(NA,dim=c(K,R))
  post.nulls=array(NA,dim=c(K,R))
  post.downs=array(NA,dim=c(K,R))

  b.mle=as.vector(t(b.gp.hat[j,]))##turn i into a R x 1 vector
  V.gp.hat=diag(se.gp.hat[j,])^2
  V.gp.hat.inv <- solve(V.gp.hat)
  
  for(k in 1:K){
    
    
    U.gp1kl <- (post.b.gpkl.cov(V.gp.hat.inv, covmat[[k]]))
    mu.gp1kl <- as.array(post.b.gpkl.mean(b.mle, V.gp.hat.inv, U.gp1kl))
    post.means[k,]=mu.gp1kl
    post.covs[k,]=as.array(diag(U.gp1kl))
    for(r in 1:R){##` marginal calculation for each tissue in each component
      if(post.covs[k,r]==0){###if the covariance matrix has entry 0, then p(null=1)
        post.ups[k,r]=0#need to figure out how to make the list have individual entries
        post.downs[k,r]=0
        post.nulls[k,r]=1}
      else{
        post.ups[k,r]=pnorm(0,mean=mu.gp1kl[r],sd=sqrt(diag(U.gp1kl)[r]),lower.tail=F)
        post.downs[k,r]=pnorm(0,mean=mu.gp1kl[r],sd=sqrt(diag(U.gp1kl)[r]),lower.tail=T)
        post.nulls[k,r]=0}
    }
  }
  return(list(post.means=post.means,post.ups=post.ups,post.downs=post.downs,post.covs=post.covs,post.nulls=post.nulls))
}


#'@title total.mean.per.snp
#'@export

total.mean.per.snp=function(post.weights,post.means){
  post.weights%*%post.means
}

#'@title total.null.per.snp
#'@export


total.null.per.snp=function(post.weights,post.nulls){
  post.weights%*%post.nulls
}

#'@title total.up.per.snp
#'@export

total.up.per.snp=function(post.weights,post.ups){
  post.weights%*%post.ups
}

#'@title total.down.per.snp
#'@export

total.down.per.snp=function(post.weights,post.downs){
  post.weights%*%post.downs
}

#'@title total.covs.partone.persnp
#'@export

total.covs.partone.persnp=function(post.means,post.covs,post.weights){
  post.weights%*%(post.covs+post.means^2)
}

#'@title lfsr.per.snp
#'@export

lfsr.per.snp=function(all.ups,all.downs){
  tailprob=rbind(all.ups,all.downs)
  as.matrix(apply(tailprob,2,function(r){1-max(r)}))
}



#'@title total.quant.per.snp
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standard errs
#'@param pis Kx1 matrix of prior weights
#'@return writes the posterior weighted quantities to a file
#'@export

total.quant.per.snp=function(j,covmat,b.gp.hat,se.gp.hat,pis,A,checkpoint=FALSE){
  gene.snp.name=rownames(b.gp.hat)[j]
  V.gp.hat=diag(se.gp.hat[j,])^2
  b.mle=b.gp.hat[j,]
  R=ncol(b.mle)
#   lik.snp=lik.func(b.mle,V.gp.hat,covmat)
#   if(sum(lik.snp*pis)==0){post.weights=rep(1/length(pis),length(pis))}
#   else{post.weights=t(lik.snp*pis/sum(lik.snp*pis))}

log.lik.snp=log.lik.func(b.mle,V.gp.hat,covmat)
log.lik.minus.max=log.lik.snp-max(log.lik.snp)
#log.pi=log(pis)
#s=log.lik.minus.max+log.pi
exp.vec=exp(log.lik.minus.max)
post.weights=t(exp.vec*pis/sum(exp.vec*pis))
                            
  all.arrays=post.array.per.snp(j,covmat,b.gp.hat,se.gp.hat)
  post.means=all.arrays$post.means
  post.ups=all.arrays$post.ups
  post.downs=all.arrays$post.downs
  post.covs=all.arrays$post.covs
  post.nulls=all.arrays$post.nulls
  
  all.nulls=total.null.per.snp(post.weights,post.nulls)
  all.mus=total.mean.per.snp(post.weights,post.means)
  all.ups=total.up.per.snp(post.weights,post.ups)
  all.downs=total.down.per.snp(post.weights,post.downs)
  lfsr=t(lfsr.per.snp(all.ups,all.downs))
  all.covs.partone=total.covs.partone.persnp(post.means,post.covs,post.weights)
  marginal.var=all.covs.partone-all.mus^2
  rownames(all.mus)=gene.snp.name
  rownames(all.ups)=gene.snp.name
  rownames(all.downs)=gene.snp.name
  rownames(lfsr)=gene.snp.name
  rownames(marginal.var)=gene.snp.name
  rownames(all.nulls)=gene.snp.name

  if(checkpoint==FALSE){
  write.table(all.mus,paste0(A,"posterior.means.txt"),append=TRUE,col.names=FALSE)
  write.table(all.ups,paste0(A,"posterior.ups.txt"),append=TRUE,col.names=FALSE)
  write.table(all.downs,paste0(A,"posterior.downs.txt"),append=TRUE,col.names=FALSE)
  write.table(all.nulls,paste0(A,"posterior.nulls.txt"),append=TRUE,col.names=FALSE)
  write.table(marginal.var,paste0(A,"marginal.var.txt"),append=TRUE,col.names=FALSE)
  write.table(lfsr,paste0(A,"lfsr.txt"),append=TRUE,col.names=FALSE)
  write.table(post.weights,paste0(A,"post.weights.txt"),append=TRUE,col.names=FALSE)}
  else{return(list(posterior.means=all.mus,posterior.nulls=all.nulls,posterior.downs=all.downs,posterior.ups=all.ups,lfsr=lfsr,marginal.var=marginal.var,post.weights=post.weights))}
}


#'@title plotting.func
#'@export

plotting.func=function(j,posterior.means,lfsr.mat,marginal.var,genesnpnames,tissue.names,mle.mat){
  
  R=ncol(posterior.means)
  posterior.means=as.matrix(posterior.means)
  col.mat=NULL
  
    for(r in 1:R){
      
      if (lfsr.mat[j,r]<=0.10) {
        col.mat[r]=1
      } else if (lfsr.mat[j,r]<0.5) {
        col.mat[r]=2
      } else if (lfsr.mat[j,r]>=0.50) {
        col.mat[r]=3
      } 
    }
  
  
  pdf(paste0("posteriormeans",genesnpnames[j],".pdf"))
  b=barplot((posterior.means[j,]),main=paste0("PostTissueMeans,",genesnpnames[j]),names=tissue.names,col=col.mat,las=2,ylim=c((min(posterior.means[j,])-0.05),(max(posterior.means[j,])+0.10)),cex.names=0.5,ylab="PosteriorMean")
  lfsr=lfsr.mat[j,]
  mean=as.numeric(posterior.means[j,])
  sd=as.numeric(sqrt(marginal.var[j,]))
  segments(b, mean - sd, b, mean + sd, lwd=2)
  dev.off()
  
  mle=mle.mat[j,]
  a.m=(mle)
  a.p=(posterior.means[j,])
  low=min(a.m-0.01,a.p-0.01)
  high=max(a.m+0.01,a.p+0.01)
  
  pdf(paste0("mlevsposteriormeans",genesnpnames[j],".pdf"))
  plot(as.matrix(mle),as.matrix(posterior.means[j,]),main=paste0("MLEvsPostMeans,",genesnpnames[j]),ylim=c(low,high),xlim=c(low,high))
  abline(0,1,col="red")
  dev.off()
  
  pdf("ShrinkagePlots.pdf")
  plot(seq(1:R),lapply(1:44,function(x){mean(abs(mle.mat[,x])>=posterior.means[,x])}),ylim=c(0,1),main="|MLE|>|Posteriormean|")
  dev.off()
}


#'@title checkfunc
#'@export
checkfunc.html=function(j,b.gp.hat,se.gp.hat,A,k ) {
  post.means=readRDS(file=paste0("post.means",A,".rds"))
  post.covs=readRDS(file=paste0("post.covs",A,".rds"))
  post.ups=readRDS(file=paste0("post.ups",A,".rds"))
  post.nulls=readRDS(file=paste0("post.nulls",A,".rds"))
  post.downs=readRDS(file=paste0("post.downs",A,".rds"))
  posterior.means=read.table(file=paste0("post.mean.",A,".txt"))
  post.weights=readRDS(paste0("post.weight.",A,".rds"))
  
  b.mle=t(as.vector(b.gp.hat[j,]))
  V.gp.hat=diag(se.gp.hat[j,])^2
  V.gp.hat.inv <- solve(V.gp.hat)
  
  U.gp1kl <- post.b.gpkl.cov(V.gp.hat.inv, covmat[[k]])
  
  par(mfrow=c(1,3))
  plot(diag(U.gp1kl),post.covs[j,k,1:R])
  
  
  
  mu.gp1kl <- as.array(post.b.gpkl.mean(b.mle, V.gp.hat.inv, U.gp1kl))
  plot(mu.gp1kl,post.means[j,k,])


  plot((post.weights[j,]%*%post.means[j,,]),posterior.means[j,])

  
}


#'@title plotting.func.html
#'@export
plotting.func.html=function(j,posterior.means,lfsr.mat,marginal.var,genesnpnames,tissue.names,mle.mat){
  
  R=ncol(posterior.means)
  posterior.means=as.matrix(posterior.means)
  col.mat=NULL
  
  for(r in 1:R){
    
    if (lfsr.mat[j,r]<=0.10) {
      col.mat[r]=1
    } else if (lfsr.mat[j,r]<0.5) {
      col.mat[r]=2
    } else if (lfsr.mat[j,r]>=0.50) {
      col.mat[r]=3
    } 
  }
  
  #par(mfrow=c(2,1))
  
  b=barplot((posterior.means[j,]),main=paste0("PostTissueMeans,",genesnpnames[j]),names=tissue.names,col=col.mat,las=2,ylim=c((min(posterior.means[j,])-0.05),(max(posterior.means[j,])+0.10)),cex.names=0.5,ylab="PosteriorMean")
  lfsr=lfsr.mat[j,]
  mean=as.numeric(posterior.means[j,])
  sd=as.numeric(sqrt(marginal.var[j,]))
  segments(b, mean - sd, b, mean + sd, lwd=2)
  t=as.numeric(levels(as.factor(col.mat)))
  key=c("<0.10","0.10<x<0.50",">0.5")
  legend("topright",legend=key[t],col=t,pch=1,title="lfsr")

  
  mle=mle.mat[j,]
  a.m=(mle)
  a.p=(posterior.means[j,])
  low=min(a.m-0.01,a.p-0.01)
  high=max(a.m+0.01,a.p+0.01)
  
  
  #plot(as.matrix(mle),as.matrix(posterior.means[j,]),main=paste0("MLEvsPostMeans,",genesnpnames[j]),ylim=c(low,high),xlim=c(low,high))
  #abline(0,1,col="red")

  
  
  #plot(seq(1:R),lapply(1:44,function(x){mean(abs(mle.mat[,x])>=posterior.means[,x])}),ylim=c(0,1),main="|MLE|>|Posteriormean|")
 
}


#'@title plotting.func.html.neg
#'@export
plotting.func.html.neg=function(j,posterior.means,lfsr.mat,marginal.var,genesnpnames,tissue.names,mle.mat){
  
  R=ncol(posterior.means)
  posterior.means=as.matrix(posterior.means)
  col.mat=NULL
  
  for(r in 1:R){
    
    if (lfsr.mat[j,r]<=0.10) {
      col.mat[r]=1
    } else if (lfsr.mat[j,r]<0.5) {
      col.mat[r]=2
    } else if (lfsr.mat[j,r]>=0.50) {
      col.mat[r]=3
    } 
  }
  

  
  b=barplot((posterior.means[j,]),main=paste0("PostTissueMeans,",genesnpnames[j]),names=tissue.names,col=col.mat,las=2,ylim=c((min(posterior.means[j,])-0.05),0),cex.names=0.5,ylab="PosteriorMean")
  lfsr=lfsr.mat[j,]
  mean=as.numeric(posterior.means[j,])
  sd=as.numeric(sqrt(marginal.var[j,]))
  segments(b, mean - sd, b, mean + sd, lwd=2)
  t=as.numeric(levels(as.factor(col.mat)))
  key=c("<0.10","0.10<x<0.50",">0.5")
  legend("topright",legend=key[t],col=t,pch=1,title="lfsr")
  
#   
#   mle=mle.mat[j,]
#   a.m=(mle)
#   a.p=(posterior.means[j,])
#   low=min(a.m-0.01,a.p-0.01)
#   high=max(a.m+0.01,a.p+0.01)
#   
#   
#   plot(as.matrix(mle),as.matrix(posterior.means[j,]),main=paste0("MLEvsPostMeans,",genesnpnames[j]),ylim=c(low,high),xlim=c(low,high))
#   abline(0,1,col="red")
#   
  
  
 # plot(seq(1:R),lapply(1:44,function(x){mean(abs(mle.mat[,x])>=posterior.means[,x])}),ylim=c(0,1),main="|MLE|>|Posteriormean|")
  
}


#'@title compute.hm.train
#'@export
compute.hm.train=function(train.b,se.train,covmat,A){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=train.b[x,],V.gp.hat=diag(se.train[x,])^2,covmat)}))
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  train=lik.mat
  pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}

#'@title compute.lik.test
#'@param b.test JxR matrix with the mles for test gene-snp pairs
#'@param J number of gene snp pairs to consider
#'@param se.test JxR matrix of standard errors of the test set
#'@param covmat K list of covariance matrices
#'@param pis K matrix of HM weights form compute.hm.train
#'@export
compute.lik.test=function(b.test,J,se.test,covmat,A,pis){
  
  J=J
  R=ncol(b.test)
  
  if(file.exists(paste0("liketest",A,".rds"))==FALSE){
    #lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=b.test[x,],V.gp.hat=diag(se.test[x,])^2,covmat)}))
    lik.mat=t(sapply(seq(1:J),function(x){log.lik.func(b.mle=b.test[x,],V.gp.hat=diag(se.test[x,])^2,covmat)}))
    saveRDS(lik.mat,paste0("liketest",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketest",A,".rds")))
  
  

  test=exp(lik.mat)
  write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
  #post.weights=as.matrix(post.weight.func(pis,lik.mat))
  #saveRDS(post.weights,paste0("post.weight.",A,".rds"))
  #rm(post.weights) ## to conserve memory
  rm(lik.mat) ## to conserve memory
  
  
} 



#'@title compute.mix.test
#'@export
compute.mix.test=function(b.gp.hat,J,se.gp.hat,covmat,A,save=FALSE){

  

  all.arrays=post.array.generator(b.gp.hat,A=A,J=J,se.gp.hat,covmat)
  
  
  
  
  if(save==TRUE){
    
    post.means=all.arrays$post.means
    post.covs=all.arrays$post.covs
    post.ups=all.arrays$post.ups
    post.downs=all.arrays$post.downs
    post.nulls=all.arrays$post.nulls
    
    # return(list(post.means,...))
    saveRDS(post.means,paste0("post.means",A,".rds"))
    saveRDS(post.covs,paste0("post.covs",A,".rds"))
    saveRDS(post.ups,paste0("post.ups",A,".rds"))
    saveRDS(post.downs,paste0("post.downs",A,".rds"))
    saveRDS(post.nulls,paste0("post.nulls",A,".rds"))}
  else(return(all.arrays))
}

#'@title test.quant
#'@export
test.quant=function(A,all.arrays){
  
  
  if(missing(all.arrays)){
    post.means=readRDS(file=paste0("post.means",A,".rds"))
    post.covs=readRDS(file=paste0("post.covs",A,".rds"))
    post.ups=readRDS(file=paste0("post.ups",A,".rds"))
    post.nulls=readRDS(file=paste0("post.nulls",A,".rds"))
    post.downs=readRDS(file=paste0("post.downs",A,".rds"))
    post.weights=readRDS(file=paste0("post.weight.",A,".rds"))
  }
  
  
  
  else{
    post.means=all.arrays$post.means
    post.covs=all.arrays$post.covs
    post.ups=all.arrays$post.ups
    post.downs=all.arrays$post.downs
    post.nulls=all.arrays$post.nulls
    post.weights=readRDS(file=paste0("post.weight.",A,".rds"))
  }
  
  J=nrow(post.means)
  all.covs.partone=t(sapply(seq(1:J),function(j){total.covs.partone(j,post.means,post.covs,post.weights)}))
  
  all.mus=t(sapply(seq(1:J),function(j){total.mean(j,post.means,post.weights)}))
  all.upper=t(sapply(seq(1:J),function(j){total.up(j,post.ups,post.weights)}))
  all.lower=t(sapply(seq(1:J),function(j){total.down(j,post.downs,post.weights)}))
  all.nuller=t(sapply(seq(1:J),function(j){total.null(j,post.nulls,post.weights)}))
  
  lfsr.mat=t(sapply(seq(1:J),function(j){compare.func(j,all.upper,all.lower)}))
  marginal.var=all.covs.partone-all.mus^2
  
  
  write.table(marginal.var,paste0("marginal.var.",A,".txt"))
  write.table(all.mus,paste0("post.mean.",A,".txt"))
  write.table(all.upper,paste0("post.up.",A,".txt"))
  write.table(all.lower,paste0("post.low.",A,".txt"))
  write.table(all.nuller,paste0("post.null.",A,".txt"))
  write.table(lfsr.mat,paste0("lfsr.",A,".txt"))
}

#'@title factor.sim
# @description Generate J beta,betahats, T statistics across R tissues from a set of covariance matrices
# @description Creates factors with magnitude based on betasd - each betaj is loaded on a minimal number of factors by simulating from rmvnorm
# @param J Number of Gene-SNP Pairs (by definition, only 0.08% of them will be true eQTL; pi0 is 0.20 and the number of genes is 1/10 J)
# @param d Number of subgroups
# @param betasd size of covariance of true effects
# @param esd standard error
#'@export

factor_sim=function(J,d=3,betasd=1,esd=0.1,K=10){
  n=trunc(0.008*J,units = 0)##number of significant gene snp paris
  F=t(sapply(seq(1:K),function(x){rnorm(d,mean=0,sd=betasd)})) 
  
  covmat=lapply(seq(1:K),function(f){t=as.matrix((F[f,]))%*%F[f,]})
  library("mvtnorm")
  library("MASS")
  possible_loadings = diag(K) #set possible loadings to be "sparse" (loaded on one factor each)
  
  z = sample(nrow(F),n,replace=TRUE)
  
  beta=t(sapply(seq(1:n),function(j){
    k=z[j]
    mvrnorm(1,mu=rep(0,d),Sigma=covmat[[k]])
  }))
  beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))
  sgp=abs(matrix(rnorm(J*d,0.11,0.001),ncol=d))##use uniform to simulate 'shrunken'
  e=t(apply(sgp,1,function(x){rmvnorm(1,mean=rep(0,d),sigma=diag(x)^2)}))
  load.e=matrix(rnorm(J*K,0,0.001),ncol=K)##to make sure all loadings are covered
  #sign.index=matrix(rbinom(n*d,size=1,prob=0.5)+1,ncol=d)
  #sign.choice=c(-1,1)
  #sign.mat=matrix(sign.choice[sign.index],ncol=d)
  #e=sign.mat*sebetahat##make some negative and some positive
  betahat = rbind(beta + e)
  #betahat=rbind(betahat,)
  tstat=betahat/abs(sgp)
  lambda=rbind(possible_loadings[z,],matrix(rep(0,(J-n)*K),ncol=K))+load.e
  return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=sgp,factors=F,lambda=lambda,tstat=tstat))
}

#' @title compute.determinant
#' @description compute product of inverse sebhats for transofrmations of likielihood
#' @export 

compute.determinant=function(se.gp.hat){
  apply(se.gp.hat,1,function(x){prod(x^-1)})
}

#' @title convert.liks
#' @description converttestlikelihoods between models
#' @export 

convert.liks=function(se.gp.hat,test.lik,pis){
  test.liks=test.lik%*%pis
  conversion.factor=compute.determinant(se.gp.hat)
  sum(log(test.liks))+sum(log(conversion.factor))
  
}

#' @title get.prior.covar.nopcs
#' @description get.prior.covar.nopcs
#' @export 
get.prior.covar.nopcs <- function(lambda.mat, Q, X.c,simat,omega.table,bma=TRUE)  {
  test=list()
  R=ncol(X.c)
  M=nrow(X.c)
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=omega*diag(1,R)
    data.prox=((t(X.c)%*%X.c)/M)
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
    full.rank=as.matrix(lambda.mat)%*%as.matrix(factor.mat)
    b=(1/M*(t(full.rank)%*%full.rank))
    b[is.nan(b)]=0
    b.norm=b/max(diag(b))
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


#' @title get.prior.covar.nopcs
#' @description get.prior.covar.nopcs
#' @export 

compute.covmat.nopc = function(b.gp.hat,sebetahat,Q,t.stat,lambda.mat,A,factor.mat,bma=TRUE){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
  
  omega.table=data.frame(omega)
  
  lambda.mat=lambda.mat
  A=A
  factor.mat=factor.mat
  R=ncol(t.stat)#number of tissues
  
  X.t=as.matrix(t.stat)
  X.real=X.t
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  
  M=nrow(X.c)
  
  U.0kl=get.prior.covar.nopcs(lambda.mat=lambda.mat,Q=Q,X.c,factor.mat=factor.mat, omega.table=omega.table,bma)
  
  
  covmat=unlist(U.0kl,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}




post.array.per.snp.list=function(j,covmat,b.gp.hat,se.gp.hat){
  
  
  R=ncol(b.gp.hat)
  K=length(covmat)
  
  b.mle=as.vector(t(b.gp.hat[j,]))
  V.j.hat=diag(se.gp.hat[j,])^2
  lik=lik.func(b.mle,V.j.hat,covmat)
  tinv=lapply(seq(1:K),function(k){solve(covmat[[k]]+V.j.hat)})
  B.j.=lapply(seq(1:K),function(k){
    diag(post.b.jk.ed.cov(tinv=tinv[[k]],covmat[[k]]))
  }
  )
  post.covs <- matrix(unlist(B.j.),ncol=R,byrow=TRUE)
  b.j.=(lapply(seq(1:K),function(k){
   b=post.b.jk.ed.mean(b.mle,tinv=tinv[[k]],covmat[[k]])##for each component, compute posterior mean
  }
  ))
  post.means=matrix(unlist(b.j.),ncol=R,byrow=TRUE) ##store as KxR matrix for each indiviudal
  
  u.t.=(lapply(seq(1:K),function(k){
    mean=b.j.[[k]];sigma=sqrt((B.j.[[k]]))
    sapply(sigma,function(x){if(x==0){0}else(pnorm(0,mean[x],x,lower.tail=FALSE))})}))
    
post.ups=matrix(unlist(u.t.),ncol=R,byrow=TRUE) ##store as KxR matrix for each indiviudal
    
l.t.=(lapply(seq(1:K),function(k){
      mean=b.j.[[k]];sigma=sqrt((B.j.[[k]]))
      sapply(sigma,function(x){if(x==0){0}else(pnorm(0,mean[x],x,lower.tail=TRUE))})}))

      post.downs=matrix(unlist(l.t.),ncol=R,byrow=TRUE) ##store as KxR matrix for each indiviudal

      null.t.=(lapply(seq(1:K),function(k){
        sigma=sqrt((B.j.[[k]]))
        sapply(sigma,function(x){if(x==0){1}else(0)})
      }))
post.nulls=matrix(unlist(null.t.),ncol=R,byrow=TRUE)

return(list(post.means=post.means,post.covs=post.covs,post.ups=post.ups,post.downs=post.downs,post.nulls=post.nulls))} ##store as KxR matrix for each indiviudal
        
      

#'@title compute.lik.test.semat
#'@param b.test JxR matrix with the mles for test gene-snp pairs
#'@param J number of gene snp pairs to consider
#'@param v.mat RxR matrix of residual variance matrix, estimated as the same for all genes
#'@param covmat K list of covariance matrices
#'@param pis K matrix of HM weights form compute.hm.train
#'@export
compute.lik.test.semat=function(b.test,J,v.mat,covmat,A,pis){
  
  J=J
  R=ncol(b.test)
  
  if(file.exists(paste0("liketest",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=b.test[x,],V.gp.hat=v.mat,covmat)}))
    saveRDS(lik.mat,paste0("liketest",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketest",A,".rds")))
  
  
  
  test=lik.mat
  write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
  post.weights=as.matrix(post.weight.func(pis,lik.mat))
  saveRDS(post.weights,paste0("post.weight.",A,".rds"))
  rm(post.weights) ## to conserve memory
  rm(lik.mat) ## to conserve memory
  
  
} 



#'@title compute.lik.test.loglik.vmat
#'@param b.test JxR matrix with the mles for test gene-snp pairs
#'@param J number of gene snp pairs to consider
#'@param v.mat RxR matrix of residual variance matrix, estimated as the same for all genes
#'@param covmat K list of covariance matrices
#'@param pis K matrix of HM weights form compute.hm.train
#'@export

compute.lik.test.loglik.vmat=function(b.test,J,v.mat,covmat,A,pis){
  
  J=J
  R=ncol(b.test)
  
  if(file.exists(paste0("liketest",A,".rds"))==FALSE){
    #lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=b.test[x,],V.gp.hat=diag(se.test[x,])^2,covmat)}))
    lik.mat=t(sapply(seq(1:J),function(x){log.lik.func(b.mle=b.test[x,],V.gp.hat=v.mat,covmat)}))
    saveRDS(lik.mat,paste0("liketest",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketest",A,".rds")))
  
  
  
  test=exp(lik.mat)
  write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
  #post.weights=as.matrix(post.weight.func(pis,lik.mat))
  #saveRDS(post.weights,paste0("post.weight.",A,".rds"))
  #rm(post.weights) ## to conserve memory
  rm(lik.mat) ## to conserve memory
  
  
} 









#'@title compute.hm.train.semat
#'@export
compute.hm.train.semat=function(train.b,v.mat,covmat,A){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=train.b[x,],V.gp.hat=v.mat,covmat)}))
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  train=lik.mat
  pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}


#'@title compute.hm.train.log.lik.vmat
#'@details here, use the matrix of residual variation that may have non-diagonal elements
#'@details v.mat is the input (as opposed to the vector of standard error) and is the same for all J genes
#'@export


compute.hm.train.log.lik.vmat=function(train.b,v.mat,covmat,A){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){log.lik.func(b.mle=train.b[x,],V.gp.hat=v.mat,covmat)}))##be sure to use log lik function
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  #train=lik.mat###but this matrix is the loglik matrix
  train=t(apply(lik.mat,1,function(x){
    e=exp(x-max(x))
    return(e)}
  ))
  #pis=mixEM.normlik(matrix_lik=train,prior=rep(1,ncol(train)))##here the matrix_lik is log normalized
  pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}






#'@title total.quant.per.snp.with.var
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standard errs
#'@param pis Kx1 matrix of prior weights
#'@return writes the posterior weighted quantities to a file
#'@export

total.quant.per.snp.with.var=function(j,covmat,b.gp.hat,se.mat,pis,A,checkpoint=FALSE){
  gene.snp.name=rownames(b.gp.hat)[j]
  b.mle=b.gp.hat[j,]
  R=ncol(b.mle)
  lik.snp=lik.func(b.mle,se.mat,covmat)
  post.weights=t(lik.snp*pis/sum(lik.snp*pis))
  
  all.arrays=post.array.per.snp.with.mat(j,covmat,b.gp.hat,se.mat)
  post.means=all.arrays$post.means
  post.ups=all.arrays$post.ups
  post.downs=all.arrays$post.downs
  post.covs=all.arrays$post.covs
  post.nulls=all.arrays$post.nulls
  
  all.mus=total.mean.per.snp(post.weights,post.means)
  all.ups=total.up.per.snp(post.weights,post.ups)
  all.downs=total.down.per.snp(post.weights,post.downs )
  lfsr=t(lfsr.per.snp(all.ups,all.downs))
  all.covs.partone=total.covs.partone.persnp(post.means,post.covs,post.weights)
  marginal.var=all.covs.partone-all.mus^2
  rownames(all.mus)=gene.snp.name
  rownames(all.ups)=gene.snp.name
  rownames(all.downs)=gene.snp.name
  rownames(lfsr)=gene.snp.name
  rownames(marginal.var)=gene.snp.name
  
  if(checkpoint==FALSE){
    write.table(all.mus,paste0(A,"posterior.means.txt"),append=TRUE,col.names=FALSE)
    write.table(all.ups,paste0(A,"posterior.ups.txt"),append=TRUE,col.names=FALSE)
    write.table(all.downs,paste0(A,"posterior.downs.txt"),append=TRUE,col.names=FALSE)
    write.table(marginal.var,paste0(A,"marginal.var.txt"),append=TRUE,col.names=FALSE)
    write.table(lfsr,paste0(A,"lfsr.txt"),append=TRUE,col.names=FALSE)
    write.table(post.weights,paste0(A,"post.weights.txt"),append=TRUE,col.names=FALSE)}
  else{return(list(posterior.means=all.mus,posterior.downs=all.downs,posterior.ups=all.ups,lfsr=lfsr,marginal.var=marginal.var))}
}


#'@title post.array.per.snp.with.mat
#'@export

post.array.per.snp.with.mat=function(j,covmat,b.gp.hat,se.mat){
  
  
  R=ncol(b.gp.hat)
  K=length(covmat)
  post.means=array(NA,dim=c(K,R))
  post.covs=array(NA,dim=c(K,R))
  post.ups=array(NA,dim=c(K,R))
  post.nulls=array(NA,dim=c(K,R))
  post.downs=array(NA,dim=c(K,R))
  
  b.mle=as.vector(t(b.gp.hat[j,]))##turn i into a R x 1 vector
  V.gp.hat=se.mat^2
  V.gp.hat.inv <- solve(V.gp.hat)
  
  for(k in 1:K){
    
    
    U.gp1kl <- (post.b.gpkl.cov(V.gp.hat.inv, covmat[[k]]))
    mu.gp1kl <- as.array(post.b.gpkl.mean(b.mle, V.gp.hat.inv, U.gp1kl))
    post.means[k,]=mu.gp1kl
    post.covs[k,]=as.array(diag(U.gp1kl))
    for(r in 1:R){##` marginal calculation for each tissue in each component
      if(post.covs[k,r]==0){###if the covariance matrix has entry 0, then p(null=1)
        post.ups[k,r]=0#need to figure out how to make the list have individual entries
        post.downs[k,r]=0
        post.nulls[k,r]=1}
      else{
        post.ups[k,r]=pnorm(0,mean=mu.gp1kl[r],sd=sqrt(diag(U.gp1kl)[r]),lower.tail=F)
        post.downs[k,r]=pnorm(0,mean=mu.gp1kl[r],sd=sqrt(diag(U.gp1kl)[r]),lower.tail=T)
        post.nulls[k,r]=0}
    }
  }
  return(list(post.means=post.means,post.ups=post.ups,post.downs=post.downs,post.covs=post.covs,post.nulls=post.nulls))
}


#' @title get.prior.covar.bmafull.only
#' @param R number of tissues
#' @param omega.table L vector grid weights
#' @return a L x K list of covariance matrices
#' @export
get.prior.covar.bmafull.only <- function(R,omega.table)  {
  test=list()
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    
      
    temp=rep(list(c(0,1)),R)
    configs = as.matrix(expand.grid(temp)[-1,])
    colnames(configs)=NULL
   
    for(c in 1:nrow(configs)) {
      mat=(configs[c,]%*%t(configs[c,]))
      test[[l]][[c]]=omega*mat}
  }
  return(U.0kl=test)
}


#' @title compute.covmat.bmafull.only
#' @param b.gp.hat a JxR matrix of betahats
#' @param sebetahat a JxR matrix of their standard errors
#' @return A list of covariance matrices
#' @export


compute.covmat.bmafull.only = function(b.gp.hat,sebetahat,A){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
  
  omega.table=data.frame(omega)
  
  
  U.0=get.prior.covar.bmafull.only(R = ncol(b.gp.hat),omega.table = omega.table)
  
  
  covmat=unlist(U.0,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}



#' @title lik.func.with.tol 
#' @details computes likelihood for each betahat
#' @param b.mle Rx1 vector of mles
#' @param V.gp.hat RxR matrix of standard errors
#' @param U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export

lik.func.with.tol=function(b.mle,R,V.gp.hat,covmat)
{ sapply(seq(1:length(covmat)),function(x){dmvnorm(x=b.mle, mu=rep(0,R),Sigma=covmat[[x]] + V.gp.hat,tol=1)})
}


#'@title compute.hm.train.with.tol
#'@export
compute.hm.train.with.tol=function(train.b,se.train,covmat,A){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func.with.tol(b.mle=train.b[x,],R,V.gp.hat=diag(se.train[x,])^2,covmat)}))
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  train=lik.mat
  pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}


compute.lik.test.with.tol=function(b.test,J,se.test,covmat,A,pis){
  
  J=J
  R=ncol(b.test)
  
  if(file.exists(paste0("liketest",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func.with.tol(b.mle=b.test[x,],R,V.gp.hat=diag(se.test[x,])^2,covmat)}))
    saveRDS(lik.mat,paste0("liketest",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketest",A,".rds")))
  
  
  
  test=lik.mat
  write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
  post.weights=as.matrix(post.weight.func(pis,lik.mat))
  saveRDS(post.weights,paste0("post.weight.",A,".rds"))
  rm(post.weights) ## to conserve memory
  rm(lik.mat) ## to conserve memory
  
  
} 


col.func=function(lfsr,posterior.means,j){
  R=ncol(posterior.means)
  
  lfsr.mat=as.matrix(lfsr)
  col.mat=NULL
  
  for(r in 1:R){
    
    if (lfsr.mat[j,r]<=0.10) {
      col.mat[r]="green"
    } else if (lfsr.mat[j,r]<0.5) {
      col.mat[r]="orange"
    } else if (lfsr.mat[j,r]>=0.50) {
      col.mat[r]="red"
    } 
  }
  return(col.mat)
}

#' @title fixpoint.play
#' @details returns pis for squareEM recognizing that if the rowSums are 0 for certain likmat entries, these need to be removed
#' @return pi

fixpoint.play=function(pi, matrix_lik, prior){  
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  badguys=which(m.rowsum==0)##if the sum of the pi*lik is zero, the denomenator for class prob will be NaN
  if(length(badguys)>0){classprob = m[-badguys,]/m.rowsum[-badguys]}
  else{classprob = m/m.rowsum}#an n by k matrix
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
  #return(classprob)
}

fixtest=function(pi, matrix_lik, prior){  
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
  
}


#'@title compute.hm.train.bma.only
#'@export
compute.hm.train.bma.only=function(train.b,se.train,covmat,A){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=train.b[x,],V.gp.hat=diag(se.train[x,])^2,covmat)}))
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  train=lik.mat
  pis=mixEMbmaonly(matrix_lik=train,prior=rep(1,ncol(train)))
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}






#' @title get.prior.covar.with.heterogeneity
#' @param R number of tissues
#' @param omega.table L vector grid weights
#' @return a L x K list of covariance matrices
#' @export
get.prior.covar.with.heterogeneity <- function(R,omega.table)  {
  test=list()
  for(l in 1:nrow(omega.table)){
    test[[l]]=list()
    omega=omega.table[l,]
    test[[l]][[1]]=diag(R)*omega
    configs=matrix(0,nrow=R,ncol=R)
    for(r in 1:R){
      configs[r,r]=1}
    configs=rbind(configs,rep(1,R),rep(1,R),rep(1,R))##the last three with model varying degrees of heterogeneity
    for(c in 1:nrow(configs)) {
      mat=(configs[c,]%*%t(configs[c,]))
      test[[l]][[c+1]]=omega*mat}
    test[[l]][[c-1]]=0.5*test[[l]][[c-1]]+diag(0.5*omega,R)
    test[[l]][[c]]=0.75*test[[l]][[c]]+diag(0.25*omega,R)}
  
  return(U.0kl=test)
}


#' @title compute.covmat.with.heterogeneity.no.shrink
#' @param b.gp.hat a JxR matrix of betahats
#' @param sebetahat a JxR matrix of their standard errors
#' @return A list of covariance matrices
#' @export


compute.covmat.with.heterogeneity.no.shrink = function(b.gp.hat,sebetahat,A,zero=FALSE){
  
  #omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
  omega=10*c(0.01,0.04,0.16,0.64,2.56)
  omega.table=data.frame(omega)
  
  
  U.0=get.prior.covar.with.heterogeneity(R = ncol(b.gp.hat),omega.table = omega.table)
  R = ncol(b.gp.hat)
  
  covmat=unlist(U.0,recursive=F)
  if(zero==TRUE){covmat=c(covmat,list(matrix(rep(0,R*R),ncol=R,nrow=R)))}
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}


compute.covmat.with.heterogeneity = function(b.gp.hat,sebetahat,A,zero=FALSE){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
  
  omega.table=data.frame(omega)
  
  
  U.0=get.prior.covar.with.heterogeneity(R = ncol(b.gp.hat),omega.table = omega.table)
  R = ncol(b.gp.hat)
  
  covmat=unlist(U.0,recursive=F)
  if(zero==TRUE){covmat=c(covmat,list(matrix(rep(0,R*R),ncol=R,nrow=R)))}
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}



