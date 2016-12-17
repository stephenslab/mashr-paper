#' @title compute.mag.by.sharing
#' @details computes the matrix of pairwise sharing by magnitude
#' @export

##thresh=0.05
##pm.mash.beta=pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
###lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]

compute.mag.by.sharing=function(lfsr.mash,thresh,pm.mash.beta){
  
shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(pm.mash.beta)
for(i in 1:ncol(lfsr.mash)){
  for(j in 1:ncol(lfsr.mash)){
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    #a=(intersect(sig.row,sig.col))
    #quotient=abs(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])##divide effect sizes
    ##divide effect sizes
    shared.fold.size[i,j]=mean(quotient>0.5&quotient<2)
    
  }}
return(shared.fold.size)}


#' @title compute.sharing.by.sig
#' @details computes the matrix of pairwise sharing by significiance
#' @export
##thresh=0.05
##lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]
compute.sharing.by.sig=function(lfsr.mash,thresh){
shared=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared)=rownames(shared.fold.size)=colnames(pm.mash.beta)
for(i in 1:ncol(lfsr.mash)){
  for(j in 1:ncol(lfsr.mash)){
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    intersect.indices=intersect(sig.row,sig.col)###gives you the shared genes
    sign.int=sign(post.means[intersect.indices,i]*post.means[intersect.indices,j])
    shared[i,j]=length(which(sign.int>0))/length(union(sig.row,sig.col))
    #shared[i,j]=length(sign.int)/length(union(sig.row,sig.col))
  }}
return(shared)
}


#' @title compute.sharing.by.sign
#' @details computes the matrix of pairwise sharing by sign
#' @export
#pm.mash.beta=pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
#lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]

compute.sharing.by.sign=function(lfsr.mash,thresh,pm.mash.beta){
shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(pm.mash.beta)

for(i in 1:ncol(lfsr.mash)){
  for(j in 1:ncol(lfsr.mash)){
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    #a=(intersect(sig.row,sig.col))
    #quotient=abs(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])##divide effect sizes
    ##divide effect sizes
    shared.fold.size[i,j]=mean(quotient>0)
    
  }}
return(shared.fold.size)
}