
#' @title sim.array.generator
#' @details simulates from multivariate normal according to componenet responsibility; counts number of times there exist discordant signs among the simulations for each gene snp pair
#' @param b.j.hat = 1xR vector of MLEs
#' @param se.j.hat=1xR vector of their standard errors
#' @param covmat K list of prior covariance matrices
#' @param pi K vector of prior weights
#' @param sim number of simulations per gene snp pair
#' @return Proportion of Time signs disagree
#' @export

sim.array.generation = function(j,b.j.hat,se.j.hat,covmat,pis,sim){
K=length(covmat)
b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
V.j.hat=diag(se.j.hat[j,]^2)
lik=lik.func(b.mle,V.j.hat,covmat)
post.weights=lik*pis/sum(lik*pis)
component=apply(rmultinom(100,1,prob = post.weights),2,function(x){which(x==1)})##choose a componenet according to responsibility
tinv=lapply(seq(1:sim),function(sim){k=component[sim];solve(covmat[[k]]+V.j.hat)})##generate a list of inverted Ks for all the simulations
b.j.=lapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.mean(b.mle,tinv=tinv[[sim]],covmat[[k]])})##for each component, compute posterior mean}
B.j.=lapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.cov(tinv=tinv[[sim]],covmat[[k]])})## a list of posterior covariances
 #and ask if they are either all pos or all neg
simulations=sapply(seq(1:sim),function(x){ 
  dat=mvrnorm(1,mu = b.j.[[x]],Sigma = B.j.[[x]],tol = 10)##for each simulation, generate a multivariate normal
              #mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  pos=sum(dat>0);neg=sum(dat<0);pos*neg!=0})##for each simulation, ask if they are all one direction; if not, assign1 for heterogeneity
 
return(sum(simulations)/length(simulations)
)}

# tim=proc.time()
# heterogeneity=sapply(seq(1:nrow(t.stat)),function(j){sim.array.generation(j,b.j.hat=maxz,v.j,covmat,pis,sim=100)})
# proc.time()-tim

#write.table(heterogeneity,"het.table.txt")
#' @title heter.mixture
#' @param lfsr = Jx(R+1)
#' @param se.j.hat=1xR vector of their standard errors
#' @param max.step=list containing two arrays: true.covs= KxRxR arrays of prior covariance matrices
#' @param pi K vector of prior weights
#' @return list of 3 arrays: JxKxR conditional posterior means, JxKxRxR posterior covariance matries and JxK normalized likelihoods (P(K|Data))
#' @export

heter.mixture=function(lfsr.mat,thresh,posterior.means){
lfsr=lfsr.mat[,-1]
genesnpnames=posterior.means[,1]
pm=posterior.means[,-1]
twos=apply(lfsr,1,function(x){sum(x<=thresh)>1})##return posterior means where the qtl was 'active in at least two tissues'
inds=which(twos==TRUE)##save indices
ptrue=pm[twos,]
rownames(ptrue)=inds 
##grab only those posterior means considered signficant at a threshold and count how many discordancies
lfs=lfsr[twos,];t=sapply(seq(1:nrow(ptrue)),function(x){l=lfs[x,];p=ptrue[x,];plow=p[which(l<thresh)];pos=sum(plow>0);neg=sum(plow<0);pos*neg!=0})
t=as.matrix(t)
##return indices
rownames(t)=inds
##heterogeneous
goodguys=inds[t==TRUE]
return(list(goodguys,genesnpnames[goodguys]))}



#' @title het.norm.array.generator
#' @details simulates from multivariate normal according to componenet responsibility; counts number of times there exist discordant signs among the simulations for each gene snp pair
#' @param b.j.hat = 1xR vector of MLEs
#' @param se.j.hat=1xR vector of their standard errors
#' @param se.beta.hat=1xR vector of standard errors of betas for conversion
#' @param covmat K list of prior covariance matrices
#' @param pi K vector of prior weights
#' @param sim number of simulations per gene snp pair
#' @return Proportion of Time signs disagree
#' @export

het.norm.array.generator = function(j,b.j.hat,se.j.hat,se.beta.hat,covmat,pis,sim,thresh){
  K=length(covmat)
  b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
  V.j.hat=diag(se.j.hat[j,]^2)
  log.lik.snp=log.lik.func(b.mle,V.j.hat,covmat)
  log.lik.minus.max=log.lik.snp-max(log.lik.snp)
  #log.pi=log(pis)
  #s=log.lik.minus.max+log.pi
  exp.vec=exp(log.lik.minus.max)
  post.weights=t(exp.vec*pis/sum(exp.vec*pis))
  standard.error=as.matrix(se.beta.hat[j,])
  component=apply(rmultinom(sim,1,prob = post.weights),2,function(x){which(x==1)})##choose a componenet according to responsibility
  tinv=lapply(seq(1:sim),function(sim){k=component[sim];solve(covmat[[k]]+V.j.hat)})##generate a sim length list of inverted Ks for all the simulations
  b.j.=t(sapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.mean(b.mle,tinv=tinv[[sim]],covmat[[k]])}))##for each component, compute posterior mean}
  b.j.list=lapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.mean(b.mle,tinv=tinv[[sim]],covmat[[k]])})##for each component, compute posterior mean}
  B.j.=lapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.cov(tinv=tinv[[sim]],covmat[[k]])})## a list of posterior covariances
  post.b.dat=t(sapply(seq(1:sim),function(x){mvrnorm(1,mu = b.j.[x,],Sigma = B.j.[[x]],tol = 10)*standard.error}))###generate a multivariate normal vector and multiply by standard error of betahat to concert to b
  consistency=t(sapply(seq(1:sim),function(x){het.func(het.norm(effectsize=t(post.b.dat[x,])),threshold = thresh)}))##convert to a 1 row matrix and count number of effects greater than 0.5
  
  return(list(consistency=consistency,component=component,b.j.=b.j.,B.j.=B.j.,post.b.dat=post.b.dat))}

# tim=proc.time()
# consistency.mat=matrix(NA,nrow=nrow(maxz),ncol=sim)
# component.mat=matrix(NA,nrow=nrow(maxz),ncol=sim)
# post.mean.j=array(dim = c(nrow(maxz),sim,ncol(maxz)))
# for(j in 1:nrow(maxz)){
#   
#   a=het.norm.array.generator(j = j,b.j.hat = maxz,se.j.hat = matrix(rep(1,ncol(maxz)*nrow(maxz)),
#                                                                     ncol = ncol(maxz)),se.beta.hat=se.beta.hat,pis = pis,sim = sim,
#                              thresh = thresh,covmat = covmat)
#   consistency.mat[j,]=a$consistency
#   component.mat[j,]=a$component
#   post.mean.j[j,,]=a$b.j.
#   post.b.data.j[j,,]=a$post.b.dat
# }


