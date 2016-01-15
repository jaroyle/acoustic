# A general question: Do you deal with DA at the level of the sample (how many samples were missed?)
# or do you deal with it at the level of the cluster.... how many sources were missed?
# this is not clear .....
#
#


mcmc.fn<- function(traps, obs.dB, xlim, ylim,nsim,nburn,cluster=FALSE, diag.plot=FALSE,clust.prior=TRUE){

# cluster = FALSE just estimates the number of sounds that were made
# if cluster=TRUE then we will try to group the sounds into unique locations


#lik( data | origin.loc , ID)
#Pr(0) = Pr(S< -3)
#Pr(non-zero) = Pr(S)  # normal likelihood
 traps<- as.matrix(traps)
##
##
## This algorithm estimates the number of unique CALL LOCATIONS
##
##
ntraps<-nrow(traps)
cutpoint<- -3
alpha<- -1
beta<- -1
sigma.s<- 0.5
psi<- 0.5
buffer<- 2  # should be input
g0<- 1
psi.clust<- 0.8

nind<- nrow(obs.dB)
M<- nrow(obs.dB) + 200    # max population size of 'calls'

if(cluster){
   Nclust<- 400              # Max populatoin size of clusters
}else{
   Nclust<- M
}

U<- matrix(NA,nrow=Nclust,ncol=2)
for(i in 1:Nclust){
  U[i,]<- c(runif(1,xlim[1],xlim[2]),runif(1,ylim[1],ylim[2]))
}

ID<- rep(NA, M)
# associate each sample with closest
for(i in 1:nind){    # nind should be nsamp
loc.tmp<- apply(matrix(traps[obs.dB[i,]!=0,],ncol=2,byrow=FALSE),2,mean)
print(loc.tmp)
dvec<-   sqrt(  (loc.tmp[1] - U[,1] )^2 + (loc.tmp[2] - U[,2])^2)
ID[i]<-  (1:Nclust)[dvec==min(dvec)][1]
}
for(i in (nind+1):M){
  ID[i]<- sample(1:Nclust, 1 )
}
# Note: Not every ID might be represented. There are some clusters that have "0 membership"
# Cluster locations
if(!cluster) ID<- 1:M   # everyone is their own cluster

# Samples == observed calls. But there were some calls not detected.
#  so the primary estimation problem is estimating the number of CALLS
#  there is uncertain ID about the calls.  So this is what we do the DA on. 
#
obs.dB<- rbind(obs.dB, matrix(0,nrow=M-nrow(obs.dB), ncol=ntraps) )
z<- c(rep(1,nind),rep(0, M-nind))

zero.guys<- c(rep(0, nind) , rep(1, M-nind)) 
Ulong<- U[ID,]  # string out the S so that there is one for each sample

 # In this situation the ID[i] is the "cluster" of sample i and each cluster has a common S[i,]
# string D out so that there is one row for each observation
D<- e2dist(U, traps) 

out<- matrix(NA,nrow=nsim,ncol=7)
colnames(out)<-c("alpha","beta","sigma.s","g0","psi.clust","psi","N")
IDout<- matrix(NA, nrow=nsim, ncol=M)
Uout<- array(NA, dim=c(nsim,Nclust,2))
zout<- matrix(NA,nrow=nsim,ncol=M)

for(sim in 1:nsim){

cat("iteration: ", sim,fill=TRUE)

# Basic parameter updates do not change with the introduction of cluster structure
mu<- alpha+ beta*D[ID,]
part2 <- dnorm(obs.dB,  mu , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]   # This is all Eq. 2 from Dawson and Efford 2009
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.mat<- loglik
loglik.curr<- sum(loglik[z==1,])

alpha.c<- rnorm(1,alpha, .3)
mu.c<- alpha.c + beta*D[ID,]
part2 <- dnorm(obs.dB,  mu.c , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu.c, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.cand<- sum(loglik[z==1,])
if(runif(1)<exp(loglik.cand-loglik.curr)){
 alpha<- alpha.c
 loglik.mat<- loglik
 loglik.curr<- loglik.cand
 mu<- mu.c
}

beta.c<- rnorm(1,beta, .3)
mu.c<- alpha+ beta.c*D[ID,]
part2 <- dnorm(obs.dB,  mu.c , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu.c, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.cand<- sum(loglik[z==1,])
if(runif(1)<exp(loglik.cand-loglik.curr)){
 loglik.mat<- loglik
 beta<- beta.c
 loglik.curr<- loglik.cand
 mu<- mu.c
}
 
sigma.c<- rnorm(1,sigma.s, .1)
part2 <- dnorm(obs.dB,  mu , sigma.c,log=TRUE)
part1<-pnorm( cutpoint, mu, sigma.c,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.cand<- sum(loglik[z==1,])
if(runif(1)<exp(loglik.cand-loglik.curr)){
 loglik.mat<- loglik
 sigma.s<- sigma.c
 loglik.curr<- loglik.cand
}

##
## update cluster membership
#  each sample needs to have it's cluster membership updated   SHOULD THIS BE OVER M ????
#  NOTE: the zero samples need to be placed too!  

if(cluster==TRUE  ){
# Ulong = reordered AND stretched-out. U = Nclust x 2 ,  Ulong = M x 2

Ulong <- U[ID,] 

for(s in 1:M){ 
if(z[s]==0) next

ID.cand<- ID
    
# pick a nearby cluster as a candidate:
dvec<- sqrt(  (Ulong[s,1] - U[,1])^2 + (Ulong[s,2] - U[,2])^2)  # length Nclust
prob.vec<- exp( -(dvec^2)/(2*1.0) ) 
prob.vec<- (prob.vec)/sum(prob.vec)  # 
# consider placing the sample with a new cluster
ID.cand[s]<- sample(1:Nclust, 1, prob=prob.vec)

J.to<- prob.vec[ID.cand[s]]
J.from<- sqrt(  (U[ID.cand[s],1] - U[,1])^2 + (U[ID.cand[s],2] - U[,2])^2 )
J.from<- exp( -(J.from^2)/(2*1.0) )
J.from<- J.from/sum(J.from)
J.from<-J.from[ID[s]]
adjust<- J.from/J.to

if(ID.cand[s]==ID[s]) next
 
#guys.in.current.ID <- ID == ID[s]   # Note this includes "s" obviously
#guys.in.candidate.ID <- ID == ID.cand[s]   # This is the current state of the candidate cluster, does not include "s"
#guys.index<- z == 1 & (guys.in.current.ID | guys.in.candidate.ID)    # these are all the guys in the from and to clusters right now
#obs.dB.current<- obs.dB[guys.index,]  # all the data for samples assigned to cluster ID[s]

obs.dB.current<- obs.dB[s,]
mu.current<-  (alpha+ beta*D[ID[s],])  #####[guys.index, ]
part2 <- dnorm(obs.dB.current,  mu.current , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu.current, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=1,ncol=ntraps)
ll.part1<-  ( part1[obs.dB.current==0] )
loglik[obs.dB.current==0]<-  ll.part1
ll.part2<-  ( part2[obs.dB.current!=0] )
loglik[obs.dB.current!=0]<- ll.part2
loglik<-rowSums(loglik)
loglik.current<- sum(loglik)
n.current <- c( sum(ID[z==1]==ID[s]) , sum(ID[z==1]==ID.cand[s]) )

if(diag.plot){
cat("----------------------------------------------------------",fill=TRUE)
cat("ID[s]: ", ID[s], fill=TRUE)
cat("n.current: ", n.current, fill=TRUE)
cat("loglik: ", loglik.current,fill=TRUE)
plot(gr,xlim=xlim,ylim=ylim)
text(S[c(ID[s],ID.cand[s]),],as.character(c(ID[s],ID.cand[s]))) 

for(xi in 1:n.current[1]){
  xx<- matrix(obs.dB[ID==ID[s],],ncol=ntraps,byrow=FALSE)[xi,]
  c.locs<- gr[xx<0,]
  points(c.locs,pch=20,col=c("black","blue","green","cyan","red","orange")[xi])
  text(gr[obs.dB[s,] <0,],"X")  # locations of sample being considered for swapping
 for(pt in 1:nrow(c.locs)){
 lines(  rbind(c.locs[pt,],
  S[ID[s],]) )
 }

} 
 }

## compute the loglike for the 2 clusters being affected, AFTER the swap
#guys.in.current.after <- ID.cand == ID[s]         
#guys.in.candidate.after <- ID.cand == ID.cand[s]  
#guys.index2<- z == 1 & (guys.in.current.after | guys.in.candidate.after)    # these are all the guys in the from and to clusters right now
#obs.dB.current2<- obs.dB[guys.index2,]  # all the data for samples assigned to cluster ID[s]

obs.dB.current2<- obs.dB[s,]
mu.current2<- (alpha+ beta*D[ID.cand[s],])  #####[guys.index2, ]
part22 <- dnorm(obs.dB.current2,  mu.current2 , sigma.s,log=TRUE)
part12<-pnorm( cutpoint, mu.current2, sigma.s,log=TRUE)
loglik2<- matrix(0,nrow=1,ncol=ntraps)
ll.part12<-  (part12[obs.dB.current2==0] )
loglik2[obs.dB.current2==0]<-  ll.part12
ll.part22<-  ( part22[obs.dB.current2!=0] )
loglik2[obs.dB.current2!=0]<-   ll.part22
loglik2<- rowSums(loglik2)
loglik.cand <- sum(loglik2)
n.prop<-     n.current + c( -1, +1) 

if(diag.plot){
  cat("ID.cand[s]: ",ID.cand[s], fill=TRUE)
  cat("n.prop: ", n.prop, fill=TRUE)
  cat("loglik: ",loglik.cand,fill=TRUE)
}

if(ID[s] == ID.cand[s]){
 if( loglik.cand != loglik.current) browser()
}

# Zero-inflated Poisson cluster size model . This has to be a model for the AUGMENTED population of clusters
if(clust.prior){
prior.curr<- sum( log( psi.clust*dpois(n.current, lambda=g0) + as.numeric(n.current==0)*(1-psi.clust) )  )
prior.cand<- sum( log( psi.clust*dpois(n.prop, lambda=g0) + as.numeric(n.prop==0)*(1-psi.clust) )  )
}else{
 prior.curr<- prior.cand<- 0
}

if(runif(1)<exp( (loglik.cand + prior.cand)-(loglik.current+prior.curr) )*adjust ){
  ## Not the right stuff to update ... need to fix this
   ID[s]<- ID.cand[s]
}

}   # end loop over samples
}  # end if cluster == TRUE

 n.c<- table(ID[z==1])
 n.c.big<- rep(0,Nclust)
 n.c.big[as.numeric(names(n.c))]<- n.c

#lik.curr<- sum( log( psi.clust*dpois(n.c.big, lambda=g0) + as.numeric(n.c.big==0)*(1-psi.clust) )  )
#g0.cand<- rnorm(1, g0, .1)
#if(g0>0){
#lik.cand<- sum( log( psi.clust*dpois(n.c.big, lambda=g0.cand) + as.numeric(n.c.big==0)*(1-psi.clust) )  )
# if(runif(1)<exp(lik.cand-lik.curr)){
#  g0<- g0.cand 
#  }
# }

g0<- sum(z)/(psi.clust*Nclust)
lik.curr<- sum( log( psi.clust*dpois(n.c.big, lambda=g0) + as.numeric(n.c.big==0)*(1-psi.clust) )  )
psi.clust.cand<- rnorm(1, psi.clust, .05)
if(psi.clust.cand>0 & psi.clust.cand <1){
g0.cand<- sum(z)/(psi.clust.cand*Nclust)
lik.cand<- sum( log( psi.clust.cand*dpois(n.c.big, lambda=g0.cand) + as.numeric(n.c.big==0)*(1-psi.clust.cand) )  )
 if(runif(1)<exp(lik.cand-lik.curr)){
  psi.clust<-psi.clust.cand
  g0<- g0.cand
  }
 }


mu<- alpha+ beta*D[ID,]
part2 <- dnorm(obs.dB,  mu , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]   # This is all Eq. 2 from Dawson and Efford 2009
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.mat<- loglik
 



U.cand<- cbind(rnorm(Nclust,U[,1],.4), rnorm(Nclust,U[,2],.4))
inbox<- U.cand[,1]< xlim[2] & U.cand[,1]> xlim[1] & U.cand[,2] < ylim[2] & U.cand[,2] > ylim[1]
U.cand[!inbox,]<- U[!inbox,]   
D.cand<- e2dist(U.cand, traps)
mu.c<- alpha+ beta*D.cand[ID,]
part2 <- dnorm(obs.dB,  mu.c , sigma.s,log=TRUE)
part1 <- pnorm( cutpoint, mu.c, sigma.s,log=TRUE)
loglik.cand<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik.cand[obs.dB==0]<- part1[obs.dB==0]
loglik.cand[obs.dB!=0]<- part2[obs.dB!=0]
#loglik.cand<- rowSums(loglik)
#loglik.cand<- rowSums(loglik)
if( sum(dim(loglik.cand)-dim(loglik.mat))!=0) browser()
loglik.diff<- loglik.cand - loglik.mat
loglik.diff<- rowSums(loglik.diff)
loglik.diff[z==0]<- 0     # I think this is right, sets rat = 1 so always accept?

#loglik.cand <- aggregate(loglik.cand, list(ID), sum)
#rat<- exp(loglik.cand[,2]- aggregate( rowSums(loglik.mat), list(ID), sum)[,2] )

loglik.diff<- aggregate(loglik.diff, list(ID), sum) 
loglik.diff2<- rep(0,Nclust)
loglik.diff2[loglik.diff[,1]]<- loglik.diff[,2]
 
rat<- exp( loglik.diff2 )
 
 
#rat[z==0]<- 1    # not sure about this. if cluster = false then this z=0 condition has to be set
                  # in general z=0 need to be removed from the aggregate above.. SEE ABOVE
swap<- runif(Nclust)< rat 
U[swap,]<- U.cand[swap,] 
D[swap,]<- D.cand[swap,]

# Recompute likelihood matrix
mu<- alpha+ beta*D[ID,]
part2 <- dnorm(obs.dB,  mu , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]   # This is all Eq. 2 from Dawson and Efford 2009
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.mat<- loglik

 

# does not need to change? This is vocalization level pr(detection)
prob0<- exp(rowSums(loglik.mat))
fc<- prob0*psi/(prob0*psi + (1-psi))
z[zero.guys==1]<- rbinom(sum(zero.guys), 1, fc[zero.guys==1])
# psi update does not need to change
psi<- rbeta(1, 1+ sum(z), 1+M-sum(z)) 

IDout[sim,]<- ID
Uout[sim,1:Nclust,1:2]<- U
out[sim,]<- c(alpha, beta, sigma.s, g0, psi.clust, psi, sum(z))
zout[sim,]<- z
}

return( list( parms = out[(nburn+1):nsim,],  
ID = IDout[(nburn+1):nsim,] , Uout=Uout[(nburn+1):nsim,,],
 zout=zout[(nburn+1):nsim,]   ) )
}


















 

