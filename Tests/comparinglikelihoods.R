##running EM

rm(list=ls())
setwd("~/Dropbox/EMSTuff/")
source("~/matrix_ash/R/main.R")
source("~/matrix_ash/R/mixEm.R")
source("~/matrix_ash/R/sarah's_mixem.R")
library("bigmemory")
library("SQUAREM")
library("mvtnorm")


##First, we load the strong t statistics from which to derive the covariance matrices#

t.stat=na.omit(read.table("~/Dropbox/AllGeneSNPStuff/max_tscore_eQTL.table.txt"))[,-c(1,2,47)]
v.j=matrix(rep(1,nrow(t.stat)*ncol(t.stat)),nrow=nrow(t.stat),ncol(t.stat))
lambda.mat=as.matrix(read.table("~/Dropbox/AllGeneSNPStuff//tri_gtex_allstrongt_lambda.out"))
factor.mat=as.matrix(read.table("~/Dropbox/AllGeneSNPStuff/tri_gtex_allstrongt_F.out"))


##Now, we perform the matric deconvolution in using the EM algorithm to estimate a 'denoised' empirical covariance matrix using the strongest 1000 gene snp pairs
max.step=deconvolution.em(t.stat = t.stat,factor.mat = factor.mat,lambda.mat = lambda.mat,K = 1,P=2,permsnp = 1000)


##we now feed this object into our compute covariance matrices structure to produce the full set of covariance matrices
covmat=compute.hm.covmat(t.stat,v.j,Q,lambda.mat,P,A,factor.mat,max.step=max.step)
    
u###We now proceeed as before, in estimating the weights on a set of largely null training data, 
start="~/Dropbox/cyclingstatistician/beta_gp_continuous/matched/firstbatch"

t.stat=na.omit(read.table(paste0(start,"t.stat.txt"),header=F,skip=1)[,-c(1,2)])
v.j=matrix(rep(1,nrow(t.stat)*ncol(t.stat)),nrow=nrow(t.stat),ncol(t.stat))

b.train=t.stat[1:1000,]
se.train=v.j[1:1000,]
compute.hm.train(train.b = b.train,se.train = se.train,covmat = covmat,A="1000trained") ##compute the HM weights on training data



###We can then compute the likelihood on the test data###
b.test=t.stat[1001:nrow(t.stat),]
se.test=v.j[1001:nrow(v.j),]

A="1000trained"
pis=readRDS(paste0("pis",A,".rds"))$pihat


compute.lik.test(b.test = b.test,J = nrow(b.test),se.test = se.test,covmat,A,pis)


###Now compute posterior quantities###
A="1000trained"
pis=readRDS(paste0("pis",A,".rds"))$pihat

rownames(b.test)=genesnpnames
weightedquants=lapply(seq(1:100),function(j){total.quant.per.snp(j,covmat,b.gp.hat=b.test,se.gp.hat = se.test,pis,A="Simulations",checkpoint = FALSE)})


