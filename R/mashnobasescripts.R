

#' @title post.cov.with.proj
#' @param L projection matrix
#' @param b.gp.hat projected L
#' @param tinv inverse of marginal covariance of b.mle
#' @param U.k prior covariance of true effect
#' @export

post.cov.with.proj=function(tinv,U.k,L){
  B.jk=U.k-U.k%*%t(L)%*%tinv%*%L%*%U.k
  return(B.jk)}


#' @title post.mean.with.proj
#' @param L projection matrix
#' @param b.mle projected L
#' @param tinv inverse of marginal covariance of b.mle
#' @param U.k prior covariance of true effect
#' @export

post.mean.with.proj=function(b.mle, tinv,U.k,L){
  b.jk=U.k%*%t(L)%*%tinv%*%b.mle
  return(b.jk)}

#' @title post.array.per.snp.no.baseline
#' @param L projection matrix
#' @param b.gp.hat projected L
#' @param se.gp.hat standard errors of unprojected estimates
#'@export

post.array.per.snp.no.baseline=function(j,covmat,b.gp.hat,se.gp.hat,L){
  
  
  R=ncol(L)
  b.mle=as.vector(t(b.gp.hat[j,]))
  K=length(covmat)
  post.means=array(NA,dim=c(K,R))
  post.covs=array(NA,dim=c(K,R))
  post.ups=array(NA,dim=c(K,R))
  post.nulls=array(NA,dim=c(K,R))
  post.downs=array(NA,dim=c(K,R))
  V.gp.hat=diag(se.gp.hat[j,])^2
  LVL=L%*%V.gp.hat%*%t(L)
  tinvlist=lapply(covmat,function(x){solve(LVL+L%*%x%*%t(L))})
  for(k in 1:K){
    
    ###generate the invertedt
    U.gp1kl <- post.cov.with.proj(tinv = tinvlist[[k]],U.k = covmat[[k]],L=L)
    mu.gp1kl <- as.array(post.mean.with.proj(b.mle = b.mle,tinv = tinvlist[[k]],U.k = covmat[[k]],L = L))
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
  return(list(post.means=post.means,post.ups=post.ups,post.downs=post.downs,post.covs=post.covs,post.nulls=post.nulls,tinvlist=tinvlist))
}


#'@title total.quant.per.snp.nobaseline
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles (let's use these as the transformed, consistnent with w = LChat, so centered)
#'@param se.gp.hat JxR matrix of standard errors of uncentered estimates
#'@param pis Kx1 matrix of prior weights
#'@param L projection matrix: default IDentity; can be the centering matrix that was used to center uncentered averages
#'@param column dimension of b.gp.hat (in the untransformed case, this is the same as the dimension of the true effects, but not in transformed case)
#'@return writes the posterior weighted quantities to a file

#'@export

total.quant.per.snp.no.baseline=function(j,covmat,b.gp.hat,se.gp.hat,pis,A,checkpoint=FALSE,L=diag(1,d)){
  gene.snp.name=rownames(b.gp.hat)[j]
  V.gp.hat=diag(se.gp.hat[j,])^2
  b.mle=b.gp.hat[j,]
  LVL=L%*%V.gp.hat%*%t(L)###redfine V.jhat as the marginal varianc of Lchat which is LVL'
  LSigL=lapply(covmat,function(x){L%*%x%*%t(L)})##redefine Sigma_p as a list of the variance of Lv (i.e., LUkL')
  log.lik.snp=log.lik.func(b.mle,LVL,LSigL)
  log.lik.minus.max=log.lik.snp-max(log.lik.snp)
  
  exp.vec=exp(log.lik.minus.max)
  post.weights=t(exp.vec*pis/sum(exp.vec*pis))
  
  all.arrays=post.array.per.snp.no.baseline(j=j,covmat=covmat,b.gp.hat=b.gp.hat,se.gp.hat=se.gp.hat,L=L)
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
  else{return(list(posterior.means=all.mus,posterior.nulls=all.nulls,posterior.downs=all.downs,posterior.ups=all.ups,lfsr=lfsr,marginal.var=marginal.var,post.weights=post.weights,LSigL=LSigL,LVL=LVL,tinvlist=all.arrays$tinvlist))}
}





#' @title compute.hm.train.log.lik.pen.with.L
#' @description = takes a matrix of training sumamry statistics and their standard errors and computes the likelihood matrix according to a list of covariance matrices, using the exponent of the log likelihood - max (llj)
#' @param train.b =  JxR matrix of training beta hats (here, already transformed)
#' @param se.train = JxR matrix of training standard errors
#' @param covmat = LxK dimenstional (unlisted list) of prior covariance matrices
#' @param  A  output file name
#' @param  Pen likelihood penalty, default1
#' @param  L = projection matrix
#' @return An object containing pis and model fit from the EM, and a pdf barplot
#' @export




compute.hm.train.log.lik.pen.with.L=function(train.b,se.train,covmat,A,pen,L){
  
  J=nrow(train.b)
  R=ncol(train.b)
  ###redfine V.jhat as the marginal varianc of Lchat which is LVL'
  LSigL=lapply(covmat,function(x){L%*%x%*%t(L)})##redefine Sigma_p as a list of the variance of Lv (i.e., LUkL')
  
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){
      V.gp.hat=diag(se.train[x,])^2;
      LVL=L%*%V.gp.hat%*%t(L);
      log.lik.func(b.mle=train.b[x,],V.gp.hat=LVL,covmat=LSigL)}))##be sure to use log lik function
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  #train=lik.mat###but this matrix is the loglik matrix
  train=t(apply(lik.mat,1,function(x){
    e=exp(x-max(x))
    return(e)}
  ))
 
  pis=mixEM(matrix_lik=train,prior=c(rep(1,ncol(train)-1),pen))
 
  
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}



#'@title compute.lik.test.withL
#'@param b.test JxR transformed estimates
#'@param J number of gene snp pairs to consider
#'@param se.test JxR matrix of standard errors of the test set
#'@param covmat K list of covariance matrices
#'@param L projection matrix
#'@param pis K matrix of HM weights form compute.hm.train
#'@export

compute.lik.test.withL=function(b.test,J,se.test,covmat,A,pis,L){
  
  J=J
  R=ncol(b.test)
  
  
  
  LSigL=lapply(covmat,function(x){L%*%x%*%t(L)})##redefine Sigma_p as a list of the variance of Lv (i.e., LUkL')
  
  
  if(file.exists(paste0("liketest",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){
      V.gp.hat=diag(se.test[x,])^2;
      LVL=L%*%V.gp.hat%*%t(L);
      log.lik.func(b.mle=b.test[x,],V.gp.hat=LVL,covmat=LSigL)}))##be sure to use log lik function
    saveRDS(lik.mat,paste0("liketest",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketest",A,".rds")))
  
  
  
  test=exp(lik.mat)
  write.table(total.lik.func(test,pis),paste0("total.lik.",A,".txt"))
  #post.weights=as.matrix(post.weight.func(pis,lik.mat))
  #saveRDS(post.weights,paste0("post.weight.",A,".rds"))
  #rm(post.weights) ## to conserve memory
  rm(lik.mat) ## to conserve memory
  
  
}


#' @title init.covmat.with.projection
#' @param t.stat matrix of transformed strong statistics (MxR) from which to derive covariance matrices (note this P will need to be RxR)
#' @param factor.mat KxR matrix of factors that have computed on strong projected ts from SFA
#' @param lambda.mat KxR matric of loadings that have been computed strong projected ts , also from SFA
#' @param K number of components to fit (i.e., the first dimension of the array), default is 3+Q
#' @param P rank of SVD approximation
#' @return KxRxR matric of prior covariance matrices to initialize the EM
#' @export


init.covmat.single.with.projection=function(t.stat=t.stat,factor.mat=factors,lambda.mat=lambda,P,Q){
  K=3+Q
  R=ncol(t.stat)#number of tissues
  true.covs=array(NA,dim=c(K,R,R))
  
  X.t=as.matrix(t.stat)
  X.real=X.t
  X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
  
  M=nrow(X.c)
  data.prox=((t(X.c)%*%X.c)/(M-1))
  true.covs[1,,]=data.prox
  full.rank=lambda.mat%*%factor.mat
  if(K>1){
    sfa.prox=(t(full.rank)%*%full.rank)/(M-1)
    true.covs[2,,]=sfa.prox
    
    svd.X=svd(X.c)
    v=svd.X$v;u=svd.X$u;d=svd.X$d
    cov.pc=1/M*v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P])%*%t(v[,1:P]%*%diag(d[1:P])%*%t(u[,1:P]))
    true.covs[3,,]=cov.pc
    
    
    
    if(Q!=0){for(q in 1:Q){
      
      load=as.matrix(lambda.mat[,q])
      fact=as.matrix(factor.mat[q,])
      rank.prox=load%*%t(fact)
      a=(1/(M-1)*(t(rank.prox)%*% rank.prox))
      
      true.covs[3+q,,]=a
    }}}
  return(true.covs)
}

#' @title deconvolution.em.bovy.with.L
#' @details wrapper to compute denoised estimates of the fuller rank covariance matrices
#' @param t.stat the R centered matrix to initalize covariance matrices (should be w plus the additional column)
#' @param w, the Mx(R-1) matrix of projected strong t statistics to fit Bovy
#' @param vj the array of LVL' (JxR-1xR-1)
#' @details note that w and vj, L should all by R-1, while t.stat and factor mat should be R
#' @param P rank of PC approxiatmion
#' @param Q rank of SFA approximation
#' @param L projection matrix, this can be Rx(R-1)
#' @return a 2 element list of K pis and the KxRxR true.covariance arrays
#' @export

deconvolution.em.with.bovy.with.L=function(t.stat,factor.mat,v.j,lambda.mat,P,L,Q,w){
  K=3+Q
  
  init.cov=init.covmat.single.with.projection(t.stat=t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,P=P,Q=Q)
  init.cov.list=list()
  for(i in 1:K){init.cov.list[[i]]=init.cov[i,,]}
  #head(init.cov.list)
  mean.mat=matrix(rep(0,ncol(t.stat)*nrow(t.stat)),ncol=ncol(t.stat),nrow=nrow(t.stat))
  
  ydata=  w
  xamp= rep(1/K,K)
  xcovar= init.cov.list
  fixmean= TRUE     
  ycovar=  v.j     
  xmean=   mean.mat   
  projection=list();for(l in 1:nrow(t.stat)){projection[[l]]=L}
  
  e=extreme_deconvolution(ydata=ydata,ycovar=ycovar,xamp=xamp,xmean=xmean,xcovar=init.cov.list,fixmean=T,projection=projection)
  
  true.covs=array(dim=c(K,R,R))
  for(i in 1:K){true.covs[i,,]=e$xcovar[[i]]}
  pi=e$xamp
  max.step=list(true.covs=true.covs,pi=pi)
  return(max.step)}




#' @title convertstandarderrors
#' @param s.j matrix of unscaled standard errors
#' @param L projection matrix, this can be Rx(R-1)
#' @return the standard errors of wj, diag(LVL)
#' @export

convertstandarderrors=function(s.j,L){
  t(apply(s.j,1,function(x){
  v.j=diag(x)^(2)
  lvl=L%*%v.j%*%t(L)
  return(sqrt(diag(lvl)))}))}

#' @title genlvlarray
#' @param s.j matrix of unscaled standard errors
#' @param L projection matrix, this can be Rx(R-1)
#' @return the JxRxR (or R-1) array of marginal variances
#' @export

genlvlarray=function(s.j,L){
  LVLarray=array(dim=c(nrow(s.j),nrow(L),nrow(L)))
  for(i in 1:nrow(s.j)){
    v=diag(s.j[i,])^2
    LVLarray[i,,]=L%*%v%*%t(L)
  }
  return(LVLarray)}






