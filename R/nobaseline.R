library("mvtnorm")
###Recall that T is the variance of the marginal distribution of $hat{b}$ (i.e., integrating over the uncertainty in B)
#tinv=solve(U.k+V.j.hat)





#'@title total.quant.per.snp.nobaseline
#'@param covmat a K list of covariance matrices
#'@param b.gp.hat JxR matrix of mles
#'@param se.gp.hat JxR matrix of standard errs
#'@param pis Kx1 matrix of prior weights
#'@param projection matrix
#'@return writes the posterior weighted quantities to a file
#'@export

total.quant.per.snp.nobaseline=function(j,covmat,b.gp.hat,se.gp.hat,pis,A,projmat=TRUE,checkpoint=FALSE){
  gene.snp.name=rownames(b.gp.hat)[j]
  V.gp.hat=diag(se.gp.hat[j,])^2
  b.mle=b.gp.hat[j,]
  R=ncol(b.mle)
  #   lik.snp=lik.func(b.mle,V.gp.hat,covmat)
  #   if(sum(lik.snp*pis)==0){post.weights=rep(1/length(pis),length(pis))}
  #   else{post.weights=t(lik.snp*pis/sum(lik.snp*pis))}
  if(projmat==TRUE){
  L=diag(R)-1/R*as.vector(rep(1,R))%*%t(as.vector(rep(1,R)))
  }
  if(projmat==FALSE){
    L=diag(R)
  }
  
  b.mle=L%*%(matrix(as.numeric(b.mle)))##center the betahats
  V.gp.hat=L%*%V.gp.hat%*%t(L)###redfine V.jhat
  covmat=lapply(covmat,function(x){L%*%x%*%t(L)})##redefine Sigma_p
  log.lik.snp=log.lik.func(t(b.mle),V.gp.hat,covmat)
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

compute.hm.train.log.lik.nobaseline=function(train.b,se.train,covmat,A,pen=FALSE,projmat=TRUE){
  
  J=nrow(train.b)
  R=ncol(train.b)
  
  if(projmat==TRUE){
    L=diag(R)-1/R*as.vector(rep(1,R))%*%t(as.vector(rep(1,R)))
  }
  if(projmat==FALSE){
    L=diag(R)
  }
  
  train.b=t(L%*%t(train.b))
  covmat=lapply(covmat,function(x){L%*%x%*%t(L)})##redefine Sigma_p
  
  
  if(file.exists(paste0("liketrain",A,".rds"))==FALSE){
    lik.mat=t(sapply(seq(1:J),function(x){log.lik.func(b.mle=train.b[x,],V.gp.hat=L%*%diag(se.train[x,])^2%*%t(L),covmat)}))##be sure to use log lik function
    
    saveRDS(lik.mat,paste0("liketrain",A,".rds"))}
  
  else(lik.mat=readRDS(paste0("liketrain",A,".rds")))
  
  #train=lik.mat###but this matrix is the loglik matrix
  train=t(apply(lik.mat,1,function(x){
    e=exp(x-max(x))
    return(e)}
  ))
  #pis=mixEM.normlik(matrix_lik=train,prior=rep(1,ncol(train)))##here the matrix_lik is log normalized
  if(pen==TRUE){
    pis=mixEM(matrix_lik=train,prior=c(rep(1,ncol(train)-1),10))
  }
  else{
    pis=mixEM(matrix_lik=train,prior=rep(1,ncol(train)))
  }
  
  saveRDS(pis,paste0("pis",A,".rds"))
  
  pdf(paste0("pis",A,".pdf"))
  barplot(t(as.matrix(pis$pihat)))
  dev.off()
}
