#' @title post.b.gpkl.cov
#' @param U.0k.l let U.0k.l represent a specific matrix in U.0kl (.e.g, U.0kl[[l]][[k]])
#' @return post.b.gpkl.cov ## returns an R*R posterior covariance matrix for jth gene snp pair
#' @return post.b.gpkl.mean return a 1 * R vector of posterior means for a given prior covariance matrix
#' @export
post.b.gpkl.cov <- function(V.gp.hat.inv, U.0k.l){
  U.gp1kl <- U.0k.l %*% solve(V.gp.hat.inv %*% U.0k.l + diag(nrow(U.0k.l)))
  return(U.gp1kl)
}


post.b.gpkl.mean <- function(b.mle, V.gp.hat.inv, U.gp1kl){
  mu.gp1kl <- U.gp1kl %*% V.gp.hat.inv %*% b.mle
  return(mu.gp1kl)
}

#' @title Get prior covariances
#' @description Compute the prior covariance matrix for each gene SNP pair in the rows and componenets in the columns
#' @param P integer number of PCs
#' @param L integer number of gridweights
#' 
#' @return U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair
#' @export
get.prior.covar.Ukl <- function(P, lambda.mat, Q, factor.mat,omega.table)  {
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
      fact=(as.matrix(fact[q,]))
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
    test[[l]][[Q+4]]=omega*b.norm}}
  return(U.0kl=test)
}



#' @title select.grid
#' @description mult.tissue.grid computes a set of omega 'stretch' factors in log 2 multiples of sd; the max and min for each tissue are chosen by autoselect.mix.sd
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


mult.tissue.grid = function(mult,betahat,sebetahat){ R=ncol(betahat);mix.weights=unlist(sapply(seq(1:ncol(b.gp.hat)),function(r){autoselect.mix.sd(betahat = betahat[,r],sebetahat = sebetahat[,r],mult=2)}))
mult=sqrt(2);sigmaamin=min(mix.weights);sigmaamax=max(mix.weights);npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult));
omega=mult^((-npoint):0) * sigmaamax;return(omega)}






##' @title lik.func computes likelihood for each betahat
##' @param b.mle Rx1 vector of mles
##' @param V.gp.hat RxR matrix of standard errors
##' @param U.0kl L dimensional list of K dimensional list with prior covairnace matrix for each grid weight, prior covariance pair

lik.func=function(b.mle,V.gp.hat,covmat)
{ sapply(seq(1:length(covmat)),function(x){dmvnorm(x=b.mle, sigma=covmat[[x]] + V.gp.hat)})
}

##' @title hm.weight.gen
##' @param lik.mat = JxK matrix of marginal likelihoods 
##' @param Jtrain = size of train set
##' @param Jtest=size of test set
##' @return compute HM weights on training and test set for compoarsion

hm.weight.gen=function(lik.mat,jtrain,jtest){

        train=lik.mat[1:jtrain,]
          test=lik.mat[jtrain:nrow(lik.mat),]
            pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
        return(pis$pihat)
}


##' @title total.lik.func
##' @param test = J x R matrix of likielihoods from test data set
##' @param pis HM weights
##' @return likelihood of whole data set
##' export
##' 
total.lik.func=function(test,pis){
  
sum(log(test%*%pis))}



##' post.weight.func converts the matrix of likelihood for each gene snp pairs to matrix of posterior weights ##` for each componenet
##' @param pis = object from EM output with prior weight P(Z=K) as computed from 
##' @param lik.mat = a JxK matrix of likelihoods (may be training set) for the P(D|Z=K)
##' @return a 1xK vector of posterior weight for each gene snp pait


post.weight.func=function(pis,lik.mat){d=t(apply(lik.mat,1,function(x){x*pis$pihat}))
                                       
                                       marg=rowSums(d)
                                       return(d/marg)}


##' @title creat an array of posterior quantities
##' @param b.mle =  Rx1 vector of beta.hats
##' @param V.gp.hat = Rx1 vector of standard errors
##' @param covmat = LxK dimenstional (unlisted list) of prior covariance matrices
##' @return post.means JxKxR array of posterior means, correspodning to the posterior mean ##' for the Jth individual in the Kth compoenent across all R tissues
##' @return post.covs JxKxR array of posterior vars, correspodning to the posterior vars ##' for the Jth individual in the Kth compoenent across all R tissues
##' @return post.ups JxKxR array of posterior tailprobs, corresponding to the marginal 
##' upper tail probability for the Jth individual in the Kth compoenent across all R 
##' @return post.ups JxKxR array of posterior tailprobs, corresponding to the marginal 
##' upper  tail probability for the Jth individual in the Kth component across all R 
##' @return post.nulls JxKxR array of posterior nullprobs, corresponding to the marginal 
##'  "null probability"" for the Jth individual in the Kth component across all R 


post.array.generator=function(b.gp.hat,se.gp.hat,covmat){
  
R=ncol(b.gp.hat)
J=nrow(b.gp.hat)
K=length(covmat)
post.means=array(NA,dim=c(J,K,R))
post.covs=array(NA,dim=c(J,K,R))
post.ups=array(NA,dim=c(J,K,R))
post.nulls=array(NA,dim=c(J,K,R))
post.downs=array(NA,dim=c(J,K,R))



for(j in 1:J){
  b.mle=as.vector(t(t.j[j,]))##turn i into a R x 1 vector
  V.gp.hat=diag(v.j[j,])^2
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




saveRDS(post.means,paste0("post.means",A,".rds"))
saveRDS(post.covs,paste0("post.covs",A,".rds"))
saveRDS(post.ups,paste0("post.ups",A,".rds"))
saveRDS(post.downs,paste0("post.downs",A,".rds"))
saveRDS(post.nulls,paste0("post.nulls",A,".rds"))
}


##ask how to write in a function that reads in arrays


#' Total mean
#'
#' Generate a K x R matrix for each gene snp pair of weighted quantities
#'
#' Generate a K x R matrix of post.weighted quantities for each gene snp pair and sum the to get total weighted posterior quantities
#' 
#' @param post.means J x K x R arrays of posterior means for each snp in each component in each tissue
#' @param post.cov J x K x R arrays of posterior variance for each snp in each component in each tissue
#' @param post.ups J x K x R arrays of posterior upper tail probabilities for each snp in each component in each tissue
#' @param post.downs J x K x R arrays of posterior lower tail probabilities for each snp in each component in each tissue
#' @param post.weights J x K matrix of posterior weights for each componenet and for each snp
#'
#' @return Get R-vector of total weighted quantities for each gene SNP pair
total.mean <- function(j, post.means, post.weights) {
    weightedmeans <- (as.matrix(post.means[j, , 1:R] * post.weights[j, ]))
    colSums(weightedmeans)
}

total.up <- function(j, post.ups, post.weights) {
    colSums(as.matrix(post.ups[j, , 1:R] * post.weights[j, ]))
}

total.down <- function(j, post.downs, post.weights) {
    colSums(as.matrix(post.downs[j, , 1:R] * post.weights[j, ]))
}

total.null <- function(j, post.nulls, post.weights) {
    colSums(as.matrix(post.nulls[j, , 1:R] * post.weights[j, ]))
} 

##' @title compares upper and lower tail probability at each tissue to determine larger
##' @param all.upper and all.lower are JxR matrices of upper and lower tail probabilities for all gene pairs 
##' @return 1 x R vector of lfsr at each tissue

compare.func=function(j){
  as.matrix(apply(rbind(all.upper[j,],all.lower[j,]),2,function(j){1-max(j)}))}

##' @title gen.XXX to compute total weight mean (1 x R vector) for each gene-snp pair
##' @return JxR matrix of grand means, grand nulls, grand lfsr 
gen.mean= function(post.means,post.weights,A){
 all.mus= t(sapply(seq(1:J),function(j){total.mean(j)}))
write.table(all.mus,paste0("post.mean.",A,".txt"))}


gen.null.probs=function(post.nulls,post.weights,A){
      all.nuller=t(sapply(seq(1:J),function(j){total.null(j)}))
      write.table(all.nuller,paste0("post.null.",A,".txt"))}

gen.tail.probs=function(post.ups,post.downs,post.weights,A){
      all.lower=t(sapply(seq(1:J),function(j){total.down(j)}))
  
      all.upper=t(sapply(seq(1:J),function(j){total.up(j)}))
  
      lfsr.mat=t(sapply(seq(1:J),function(j){compare.func(j)}))
  
  
      write.table(all.lower,paste0("post.low.",A,".txt"))
      write.table(all.upper,paste0("post.ups.",A,".txt"))
      write.table(lfsr.mat,paste0("lfsr.",A,".txt"))
}






compute.covmat = function(b.gp.hat,sebetahat,Q,X.c,lambda.mat,P,A,factor.mat){
    omega=mult.tissue.grid(mult=sqrt(2),b.gp.hat,sebetahat)
    omega.table=data.frame(omega)
    lambda=lambda
    A=A
    factor.mat=factor.mat
    X.c=X.c
    Q=Q
    
    U.0kl=get.prior.covar.Ukl(P=2,lambda=lambda,Q=Q,factor.mat=factor.mat, omega.table=omega.table)
    
    covmat=unlist(U.0kl,recursive=F)
    
    saveRDS(covmat,paste0("covmat",A,".rds"))
    
    return(covmat)}



compute.mixture.dist=function(b.gp.hat,se.gp.hat,covmat,A=t.stat.inference){
    
    
    
    
    lik.mat=t(sapply(seq(1:J),function(x){lik.func(b.mle=t.j[x,],V.gp.hat=diag(v.j[x,])^2,covmat)}))
    
    saveRDS(lik.mat.simulated.sfa,paste0("likelihood",A,".rds"))
    
    
    pis=hm.weight.gen(lik.mat,J/2,J.2)
    
    
    total.lik.func(lik.mat,pis)
    
    post.weights=as.matrix(post.weight.func(pis,lik.mat))
    
    saveRDS(post.weights,paste0("post.weight.",A,".rds"))
    
    
    post.array.generator(b.gp.hat,se.gp.hat,covmat)
    
    
    gen.mean(post.means,post.weights,A)
    gen.null.probs(post.nulls,post.weights,A)
    gen.tail.probs(post.means,post.weights,A)
}





