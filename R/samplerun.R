rm(list=ls())
source("~/matrix_ash/R/main.R")
source("/mnt/lustre/home/surbut/prior_mix_multi/matched/mixEm.R")
start="firstbatch"

b.gp.hat=na.omit(read.table(paste0("/mnt/lustre/home/surbut/prior_mix_multi/matched/",start,"beta.hat.unstd.txt"),header=F,skip=1)[,-c(1,2)])
se.gp.hat=na.omit(read.table(paste0("/mnt/lustre/home/surbut/prior_mix_multi/matched/",start,"se.beta.hat.txt"),header=F,skip=1)[,-c(1,2)])
t.j=na.omit(read.table(paste0("/mnt/lustre/home/surbut/prior_mix_multi/matched/",start,"t.stat.txt"),header=F,skip=1)[,-c(1,2)])
v.j=matrix(rep(1,nrow(t.j)*ncol(t.j)),nrow=nrow(t.j),ncol(t.j))
t.stat=na.omit(read.table("/mnt/lustre/home/surbut/prior_mix_multi/matched/filtered_SNP_batcht.stat.txt",header=F,skip=1)[,-c(1,2,3,48)])

library("mvtnorm")


R=ncol(b.gp.hat)#number of tissues
X.t=as.matrix(t.stat)
#X.real=X.t[which(truth$config!=0),]
X.real=X.t
X.c=apply(X.real,2,function(x) x-mean(x)) ##Column centered matrix of t statistics
R=ncol(X.t)
M=nrow(X.c)


lambda=as.matrix(read.table("/mnt/lustre/home/surbut/prior_mix_multi/matched/tri_gtex_strongmay10_lambda.out"))
factor=as.matrix(read.table("/mnt/lustre/home/surbut/prior_mix_multi/matched/tri_gtex_strongmay10_F.out"))
covmat=compute.covmat(t.j,v.j,Q=5,X.c,lambda,P=2,A="fiveeighteen",factor)
compute.mixture.dist(t.j,J=nrow(t.j),v.j,covmat,A="fiveeighteen")
compute.total.quant(A="fiveeighteen",J=nrow(t.j))
checkfunc(1,t.j,v.j,A="fiveeighteen",100)

##And if you wanted to plot##

genesnpnames=na.omit(read.table("~/Dropbox/AllGeneSNPStuff/max_beta_eQTL.table.txt"))[,2]
tissuenames=read.table("~/matrix_ash/tissuenames.txt")[,1]
marginal.var=read.table(paste0("marginal.var."A,".txt"))
lfsr.mat=read.table(paste0("lfsr.",A,".txt"))
posterior.means=read.table(paste0("post.mean.",A,".txt"))
plotting.func(186,posterior.means,lfsr.mat,marginal.var,genesnpnames,tissue.names)
