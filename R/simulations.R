
#' @title sim.array.generator
#' @param b.j.hat = 1xR vector of MLEs
#' @param se.j.hat=1xR vector of their standard errors
#' @param max.step=list containing two arrays: true.covs= KxRxR arrays of prior covariance matrices
#' @param pi K vector of prior weights
#' @return list of 3 arrays: JxKxR conditional posterior means, JxKxRxR posterior covariance matries and JxK normalized likelihoods (P(K|Data))
#' @export

sim.array.generation = function(j,covmat,pis,sim){
K=length(covmat)
b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
V.j.hat=diag(se.j.hat[j,]^2)
lik=lik.func(b.mle,V.j.hat,covmat)
post.weights=lik*pis/sum(lik*pis)
component=apply(rmultinom(100,1,prob = post.weights),2,function(x){which(x==1)})##choose a componenet according to responsibility
tinv=lapply(seq(1:sim),function(sim){k=component[sim];solve(covmat[[k]]+V.j.hat)})##generate a list of inverted Ks for all the simulations
b.j.=lapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.mean(b.mle,tinv=tinv[[sim]],covmat[[k]])})##for each component, compute posterior mean}
B.j.=lapply(seq(1:sim),function(sim){ k=component[sim];post.b.jk.ed.cov(tinv=tinv[[sim]],covmat[[k]])})## a list of posterior covariances


#example: cov[43,44]
##for each simulation, generate a multivariate normal and ask if they are either all pos or all neg
simulations=sapply(seq(1:sim),function(x){ 
  dat=mvrnorm(1,mu = b.j.[[x]],Sigma = B.j.[[x]],tol = 10)
              #mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  pos=sum(dat>0);neg=sum(dat<0);pos*neg!=0})
##count the number of times they are both negative and positive 
return(sum(simulations)/length(simulations))}


twos=apply(lfsr,1,function(x){sum(x<=0.05)>1})##return posterior means where the qtl was 'active in at least two tissues'

inds=which(twos==TRUE)##save indices

ptrue=pm[twos,]
rownames(ptrue)=inds
lfs=lfsr[twos,];t=sapply(seq(1:nrow(ptrue)),function(x){l=lfs[x,];p=ptrue[x,];plow=p[which(l<0.05)];pos=sum(plow>0);neg=sum(plow<0);pos*neg!=0})


t=as.matrix(t)

##return indices
rownames(t)=inds
##heterogeneous
inds[t==TRUE]
