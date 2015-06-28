

##################
###################
#Testing Section###

##################
###################
J=100;R=6;K=3

source('~/matrix_ash/R/sarah\'s_mixem.R')
b.j.hat=matrix(rnorm(J*R),ncol=R)
se.j.hat=matrix(rep(1,J*R),ncol=R)
factor.mat=matrix(rnorm(K*R),ncol=R)
lambda.mat=matrix(rnorm(K*J),ncol=K)

par.init=list(true.covs=init.covmat(t.stat = b.j.hat,factor.mat = factor.mat,lambda.mat = lambda.mat,K = 3,P=2),pi=rep(1/K,K))###output a list of K covariance matrices and initial pis 
par.init.unlist=unlist(par.init)
K=dim(par.init$true.covs)[1]
prior=rep(1,K)



###To test, set 
max.step=par.init
(dim(max.step$true.covs))##check to make sure [K,R,R]
(length(max.step$pi))
niter=100
neglik=rep(0,niter)
max.step.unlist=par.init.unlist
##and then run the fixpoint function for the first iteration##
for(i in 1:niter){
  a=fixpoint.cov(max.step.unlist,b.j.hat,se.j.hat)
  neglik[i]=negpenlogliksarah(max.step.unlist=a,b.j.hat = b.j.hat,se.j.hat = se.j.hat)
  max.step.unlist=a
}

##the neg lik should be deceasing
plot(neglik)
##test to make sure that the results of computations are stored properly for the jth individual and kth componenet
test.funct=function(j,max.step.unlist,k,R){
  
  L=length(max.step.unlist)
  K=L/(R^2+1)
  dim.true.covs=c(K,R,R)
  pi.length=K
  max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])
  true.covs=max.step$true.covs
  pi=max.step$pi
  
  b.mle=as.vector(t(b.j.hat[j,]))##turn i into a R x 1 vector
  V.j.hat=diag(se.j.hat[j,]^2)
  lik=sapply(seq(1:K),function(k){dmvnorm(x=b.mle, sigma=true.covs[k,,] + V.j.hat)})##compute K element likeilihood for each idndiviual
  B.j.=(lapply(seq(1:K),function(k){
    tinv=solve(true.covs[k,,]+V.j.hat)##covariance matrix of the marginal distribution
    post.b.jk.ed.cov(tinv=tinv,true.covs[k,,])
  }
  ))##create a K dimensional list of covariance matrices 
  
  ##compute a K dimensional list of posterior means for each J
  b.j.=(lapply(seq(1:K),function(k){
    tinv=solve(true.covs[k,,]+V.j.hat)##covariance matrix of the marginal distribution
    post.b.jk.ed.mean(b.mle,tinv=tinv,true.covs[k,,])##for each component, compute posterior mean
  }
  ))
  
  
  
  test=em.array.generator(max.step = max.step,b.j.hat = b.j.hat,se.j.hat = se.j.hat)
  pm=test$post.means;pc=test$post.covs;q.mat=test$q.mat
  
  par(mfrow=c(1,2))
  plot(pm[j,k,],b.j.[[k]])##test to make sure posterior mean is stored properly
  plot(diag(pc[j,k,,]),diag(B.j.[[k]]))##test to make sure posterior covariance is stored properly
}

test.funct(j = 1,max.step.unlist = par.init.unlist,k=3,R=R)
#s=squarem(par=par.init.unlist,b.j.hat=b.j.hat,se.j.hat=se.j.hat,fixptfn=fixpoint.cov, objfn=negpenlogliksarah)
max.step.unlist=s$par
dim.true.covs=c(K,R,R)
max.step = list(true.covs = array(max.step.unlist[1:prod(dim.true.covs)], dim = dim.true.covs), pi = max.step.unlist[(prod(dim.true.covs)+1):(prod(dim.true.covs)+pi.length)])


dim.true.cov.fun=function(b.j.hat,se.j,hat,max.step.unlist){
  L=length(max.step.unlist)
  K=L/(R^2+1)
  dim.true.covs=c(K,R,R)
  pi.length=K
  return(list(dim.true.covs=dim.true.covs,pi.length=pi.length))
  }
