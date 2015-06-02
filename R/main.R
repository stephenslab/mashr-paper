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
#' @param L integer number of gridweights
#' @return U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export
get.prior.covar.Ukl.with.rho <- function(P, lambda.mat, Q, factor.mat,omega.table,rhos)  {
    test=list()
    id=diag(1,R)
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
#' @export
get.prior.covar.Ukl <- function(P, lambda.mat, Q, factor.mat,omega.table,bma=TRUE)  {
  test=list()
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
      fact=factor.mat
      load=as.matrix(lambda[,q])
      fact=t((as.matrix(fact[q,])))
      rank.prox=load%*%fact
      a=(1/M*(t(rank.prox)%*% rank.prox))
      a[is.nan(a)] = 0
      a.norm=a/max(diag(a))
      test[[l]][[q+3]]=omega*a.norm
    }
    full.rank=as.matrix(lambda)%*%as.matrix(factor.mat)
    b=(1/M*(t(full.rank)%*%full.rank))
    b[is.nan(b)]=0
    b.norm=b/max(diag(b))
    test[[l]][[Q+4]]=omega*b.norm}
    
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


#' @title autoselect.mix.sd 
#' @param beta hat is an J x R matrix of betahats (or t statistics)
#' @param sebetahat is a J x R matrix of their standard errors (or 1s, if T staitsitcs)
#' @return omega, a list of stretch factors by which to scale each covariance matrix U_k
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
#' @export
mult.tissue.grid = function(mult,betahat,sebetahat){ R=ncol(betahat);mix.weights=unlist(sapply(seq(1:ncol(betahat)),function(r){autoselect.mix.sd(betahat = betahat[,r],sebetahat = sebetahat[,r],mult=2)}))
    mult=sqrt(2);sigmaamin=min(mix.weights);sigmaamax=max(mix.weights);npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult));
    omega=mult^((-npoint):0) * sigmaamax;return(omega)}






#' @title lik.func 
#' @details computes likelihood for each betahat
#' @param b.mle Rx1 vector of mles
#' @param V.gp.hat RxR matrix of standard errors
#' @param U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export

lik.func=function(b.mle,V.gp.hat,covmat)
{ sapply(seq(1:length(covmat)),function(x){dmvnorm(x=b.mle, sigma=covmat[[x]] + V.gp.hat)})
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
#' @describe generate a K x R matrix of post.weighted quantieis for each gene snp pair and sum them to get total weighted
#' @param post.means J x K x R arrays of posterior means for each snp in each component in each tissue
#' @param post.cov J x K x R arrays of posterior variance for each snp in each component in each tissue
#' @param post.ups J x K x R arrays of posterior upper tail probabilities for each snp in each component in each tissue
#' @param post.downs J x K x R arrays of posterior lower tail probabilities for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#' @return get vector of total weighted quanties for each gene SNP pair
#' @export

total.mean=function(j,post.means,post.weights){weightedmeans=(as.matrix(post.means[j,,1:R]*post.weights[j,]))
    colSums(weightedmeans)}

#' @title total.up
#' @export
total.up=function(j,post.ups,post.weights){colSums(as.matrix(post.ups[j,,1:R]*post.weights[j,]))}

#' @title total.down
#' @export
total.down=function(j,post.downs,post.weights){colSums(as.matrix(post.downs[j,,1:R]*post.weights[j,]))}

#' @title total.null
#' @export
total.null=function(j,post.nulls,post.weights){colSums(as.matrix(post.nulls[j,,1:R]*post.weights[j,]))}

#' @title total.covs.partone
#' @export
total.covs.partone=function(j,post.means,post.covs,post.weights){colSums(as.matrix((post.covs[j,,1:R]+post.means[j,,1:R]^2)*post.weights[j,]))}




#' @title compare.func
#' @details compares upper and lower tail probability at each tissue to determine larger
#' @param all.upper and all.lower are JxR matrices of upper and lower tail probabilities for all gene pairs
#' @return 1 x R vector of lfsr at each tissue

compare.func=function(j,all.upper,all.lower){
    as.matrix(apply(rbind(all.upper[j,],all.lower[j,]),2,function(j){1-max(j)}))}







#' @title compute.covmat.with.rho
#' @export

compute.covmat.with.rho = function(b.gp.hat,sebetahat,Q,X.c,lambda.mat,P,A,factor.mat){
  
    omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
    
    omega.table=data.frame(omega)
    rhos=seq(0.1,1,by=0.1)
    lambda=lambda
    A=A
    factor.mat=factor.mat
    X.c=X.c
    Q=Q
    
 
    U.0kl=get.prior.covar.Ukl.with.rho(P=2,lambda=lambda,Q=Q,factor.mat=factor.mat, omega.table=omega.table,rhos=rhos)
    
    covmat1=unlist(U.0kl,recursive=F)
    covmat=unlist(covmat1,recursive=F)
    saveRDS(covmat,paste0("covmat",A,".rds"))
    
    return(covmat)}


#' @title compute.covmat
#' @return A list of covariance matrices
#' @export

compute.covmat = function(b.gp.hat,sebetahat,Q,X.c,lambda.mat,P,A,factor.mat,bma=TRUE){
  
  omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
  
  omega.table=data.frame(omega)
 
  lambda=lambda
  A=A
  factor.mat=factor.mat
  X.c=X.c
  Q=Q
  
U.0kl=get.prior.covar.Ukl(P=2,lambda=lambda,Q=Q,factor.mat=factor.mat, omega.table=omega.table,bma)
    
   
    covmat=unlist(U.0kl,recursive=F)
  saveRDS(covmat,paste0("covmat",A,".rds"))
  
  return(covmat)}


#' @title compute.mixture.dist
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

#'@title
#'@export

post.array.per.snp=function(j,covmat,b.gp.hat,se.gp.hat){


  
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

total.mean.per.snp=function(j,post.weights,post.means){
  post.weights[j,]%*%post.means
}


total.null.per.snp=function(j,post.weights,post.nulls){
  post.weights[j,]%*%post.nulls
}


total.up.per.snp=function(j,post.weights,post.ups){
  post.weights[j,]%*%post.ups
}


total.down.per.snp=function(j,post.weights,post.downs){
  post.weights[j,]%*%post.downs
}


total.covs.partone.persnp=function(j,post.means,post.covs,post.weights){
  post.weights[j,]%*%(post.covs+post.means^2)
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
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=b.train[x,],V.gp.hat=diag(se.train[x,])^2,covmat)}))
    
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
#'@export
compute.lik.test=function(b.gp.hat,J,se.gp.hat,covmat,A,pis){
  
  J=J
  R=ncol(b.gp.hat)
  
  if(file.exists(paste0("liketest",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=b.gp.hat[x,],V.gp.hat=diag(se.gp.hat[x,])^2,covmat)}))
    
    saveRDS(lik.mat,paste0("liketest",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketest",A,".rds")))
  
  

  test=lik.mat
  write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
  post.weights=as.matrix(post.weight.func(pis,lik.mat))
  saveRDS(post.weights,paste0("post.weight.",A,".rds"))
  rm(post.weights) ## to conserve memory
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
#'@export

factor_sim=function(n,d=3,betasd,esd=0.3,K=10){
  
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
  sebetahat=abs(matrix(runif(n*d,esd-0.05,esd+0.05),ncol=d))
  betahat = beta + sebetahat
  tstat=betahat/sebetahat
  lambda=possible_loadings[z,]
  return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=sebetahat,factors=F,lambda=lambda,tstat=tstat))
}







