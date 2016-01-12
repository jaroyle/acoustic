# A general question: Do you deal with DA at the level of the sample (how many samples were missed?)
# or do you deal with it at the level of the cluster.... how many sources were missed?
# this is not clear .....
#
#


mcmc.fn<- function(traps, obs.dB, xlim, ylim,nsim,nburn,cluster=FALSE){

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

nind<- nrow(obs.dB)
M<- nrow(obs.dB) + 200    # max population size of 'calls'

if(cluster){
   Nclust<- 300              # Max populatoin size of clusters
}else{
   Nclust<- M
}

u<- matrix(NA,nrow=Nclust,ncol=2)
for(i in 1:Nclust){
  u[i,]<- c(runif(1,xlim[1],xlim[2]),runif(1,ylim[1],ylim[2]))
}
#for(i in 1:nrow(obs.dB)){
# u[i,]<-apply(traps[obs.dB[i,]!=0,],2,mean)
#}

ID<- rep(NA, M)
# associate each sample with closest
for(i in 1:nind){    # nind should be nsamp
loc.tmp<- apply(matrix(traps[obs.dB[i,]!=0,],ncol=2,byrow=FALSE),2,mean)
print(loc.tmp)
dvec<-   sqrt(  (loc.tmp[1] - u[,1] )^2 + (loc.tmp[2] - u[,2])^2)
ID[i]<-  (1:M)[dvec==min(dvec)][1]
}
for(i in (nind+1):M){
ID[i]<- sample(1:Nclust, 1 )
}
# Note: Not every ID might be represented. There are some clusters that have "0 membership"
#
# Cluster locations
S<- u
if(!cluster) ID<- 1:M   # everyone is their own cluster

# Samples == observed calls. But there were some calls not detected.
#  so the primary estimation problem is estimating the number of CALLS
#  there is uncertain ID about the calls.  So this is what we do the DA on. 
#
obs.dB<- rbind(obs.dB, matrix(0,nrow=M-nrow(obs.dB), ncol=ntraps) )
z<- c(rep(1,nind),rep(0, M-nind))

zero.guys<- c(rep(0, nind) , rep(1, M-nind)) 
u<- S[ID,]  # string out the S so that there is one for each sample

 # In this situation the ID[i] is the "cluster" of sample i and each cluster has a common S[i,]
# string D out so that there is one row for each observation
D<- e2dist(S, traps)[ID,]

out<- matrix(NA,nrow=nsim,ncol=5)
colnames(out)<-c("alpha","beta","sigma.s","psi","N")


for(sim in 1:nsim){


# Basic parameter updates do not change with the introduction of cluster structure
mu<- alpha+ beta*D
part2 <- dnorm(obs.dB,  mu , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]   # This is all Eq. 2 from Dawson and Efford 2009
loglik[obs.dB!=0]<- part2[obs.dB!=0]
loglik.curr<- sum(loglik[z==1,])

alpha.c<- rnorm(1,alpha, .3)
mu.c<- alpha.c + beta*D
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
mu.c<- alpha+ beta.c*D
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
##    only have to for each observed individual.... right?
##
if(cluster==TRUE){

for(s in 1:nind){  # each sample needs to have it's cluster membership updated
   # NOTE: the zero samples need to be placed too!  

ID.cand<- ID
    # u = location of cluster that sample s is currently in 
    #  S = cluster locations
# pick a nearby cluster as a candidate:
dvec<- sqrt(  (u[s,1] - S[,1])^2 + (u[s,2] - S[,2])^2)  # length Nclust
prob.vec<- exp( -(dvec^2)/(2*1.0) ) 
prob.vec<- (prob.vec)/sum(prob.vec)  # 
# consider placing the sample with a new cluster
ID.cand[s]<- sample(1:Nclust, 1, prob=prob.vec)

from.clust<- ID[s]   # cluster currently assigned to
prop.clust<- ID.cand[s]  # candidate cluster   NOTE: proposed cluster could be an empty cluster Thus no "splitting" is needed. 
 

guys.in.current.ID <- ID == ID[s]   # Note this includes "s" obviously
guys.in.candidate.ID <- ID == ID.cand[s]   # This is the current state of the candidate cluster, does not include "s"

guys.index<- z == 1 & (guys.in.current.ID | guys.in.candidate.ID)    # these are all the guys in the from and to clusters right now

obs.dB.current<- obs.dB[guys.index,]  # all the data for samples assigned to cluster ID[s]
mu.current<- alpha+ beta*D[guys.index, ]
part2 <- dnorm(obs.dB.current,  mu.current , sigma.s,log=TRUE)
part1<-pnorm( cutpoint, mu.current, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=sum(guys.index),ncol=ntraps)
loglik[obs.dB.current==0]<- part1[obs.dB.current==0]
loglik[obs.dB.current!=0]<- part2[obs.dB.current!=0]
loglik.current<- sum(loglik)

## Now I  need to compute the loglike for the 2 clusters being affected, AFTER the swap


if(runif(1)<exp(loglik.cand-loglik.curr)){
 loglik.mat<- loglik
 beta<- beta.c
 loglik.curr<- loglik.cand
 mu<- mu.c
}
 


#Need to compute the likelihood of the current clusters and the likelihood of the proposed clusters and
#then think about making the swap using a MH rule



}   # end loop over samples
}  # end if cluster == TRUE

#
#  the loglike matrix needs to be updated here....
# 
## Just updating logliklihood object
## This is a waste!
#part2 <- dnorm(obs.dB,  mu , sigma.s,log=TRUE)
#part1<-pnorm( cutpoint, mu, sigma.s,log=TRUE)
#loglik.mat<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
##
## 



S.cand<- cbind(rnorm(Nclust,S[,1],.4), rnorm(Nclust,S[,2],.4))
inbox<- S.cand[,1]< xlim[2] & S.cand[,1]> xlim[1] & S.cand[,2] < ylim[2] & S.cand[,2] > ylim[1]
S.cand[!inbox,]<- S[!inbox,]   
D.cand<- e2dist(S.cand[ID,], traps)
mu.c<- alpha+ beta*D.cand

part2 <- dnorm(obs.dB,  mu.c , sigma.s,log=TRUE)
part1 <- pnorm( cutpoint, mu.c, sigma.s,log=TRUE)
loglik<- matrix(0,nrow=nrow(obs.dB),ncol=ntraps)
loglik[obs.dB==0]<- part1[obs.dB==0]
loglik[obs.dB!=0]<- part2[obs.dB!=0]
#loglik.cand<- rowSums(loglik)
#loglik.cand<- rowSums(loglik)

if( sum(dim(loglik)-dim(loglik.mat))!=0) browser()
loglik.diff<- loglik - loglik.mat
loglik.diff<- rowSums(loglik.diff)
loglik.diff[z==0]<- 0     # I think this is right, sets rat = 1 so always accept?

#loglik.cand <- aggregate(loglik.cand, list(ID), sum)
#rat<- exp(loglik.cand[,2]- aggregate( rowSums(loglik.mat), list(ID), sum)[,2] )

loglik.diff<- aggregate(loglik.diff, list(ID), sum) 
rat<- exp( loglik.diff[,2] )
idx<- rep(1,Nclust)
idx[loglik.diff[,1]]<- rat
rat<-idx
#rat[z==0]<- 1    # not sure about this. if cluster = false then this z=0 condition has to be set
                  # in general z=0 need to be removed from the aggregate above.. SEE ABOVE
swap<- runif(Nclust)< rat 
S[swap,]<- S.cand[swap,] 
#######
# These are of differing dimensions which needs resolved if cluster=TRUE (swap = Nclust x 1)
########
loglik.mat[swap,]<- loglik[swap,]
D[swap,]<- D.cand[swap,]
mu<- mu.c[swap,]
 

# does not need to change? This is vocalization level pr(detection)
prob0<- exp(rowSums(loglik.mat))
fc<- prob0*psi/(prob0*psi + (1-psi))
z[zero.guys==1]<- rbinom(sum(zero.guys), 1, fc[zero.guys==1])
# psi update does not need to change
psi<- rbeta(1, 1+ sum(z), 1+M-sum(z)) 


out[sim,]<- c(alpha, beta, sigma.s, psi, sum(z))

}

return(out[(nburn+1):nsim,] )
}
 

