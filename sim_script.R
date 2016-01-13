library(scrbook)
simout<- matrix(NA,nrow=20,ncol=6)

for(iter in 1:20){
# trap locations (acoustic detectors)
gr<- expand.grid(2:11,2:11)

xlim<- c(0,13)
ylim<- c(0,13)

# Individuals
N<- 40
sigma.move<- 0.4   # Movement about home range center. The standard SCR scale parameter

# activity centers
Sx<- runif(N, xlim[1], xlim[2])
Sy<- runif(N, ylim[1], ylim[2])

# number of singing locations of each individual. Can be arbitrary, we 
# build the model conditional on the number of detections.  Note: I'm thinking a parametric model
# might be important even necessary for estimation.... This is an issue to keep in mind.
#nmoves<- rpois(N, lambda=3)
# Here I allow for at most 8 locations of singing
ncall.locs<- sample(1:8, N,replace=TRUE, prob=c(.3, rep(.1, 7)))

# This is the total number of locations from where calls were made.
N.sources<- sum(ncall.locs)


#  Now we loop through N and for each location of each individual we simulate up to 4 calls
#   made from that location. This is to mimic sampling over a time interval.
#   As noted above, a parametric model might be important/necessary here.
locs.x<- NULL
locs.y<- NULL
ind.ID<- source.ID<- NULL
for(i in 1:N){
  if(ncall.locs[i]==0) next   # this guy not detected at all
  ncalls.x<- rnorm(ncall.locs[i], Sx[i], sigma.move)
  ncalls.y<- rnorm(ncall.locs[i], Sy[i], sigma.move)
    # next line is the number of calls from each location
  nbr.calls<- sample(1:4, ncall.locs[i], replace=TRUE,prob=c(.25, .5, .2, .05) )
  locs.x<- c(locs.x,   rep(ncalls.x, nbr.calls) )  # nbr.calls at each location
  locs.y<- c(locs.y,   rep(ncalls.y, nbr.calls) ) 
  ind.ID<- c(ind.ID, rep(i, sum(nbr.calls) ) )
  source.ID<- c(source.ID,rep(paste(i,1:ncall.locs[i],sep="."), nbr.calls))
}
locs<- cbind(locs.x, locs.y)

# This is the total number of acoustic signals made ("tweets")
N.signals<- length(source.ID)

# Now simulate some acoustic sampling data

ID.true<- as.numeric(factor(source.ID))

# dist from location of each signal to recorders
D<- e2dist(locs, gr)

# decibels picked up in each recorder. There are 3 parameters to be estimated here
dB <-   -1 + -2*D + rnorm(prod(dim(D)), 0, .5)

# was messing with graphical summaries here, got bored.
if(1==2){
for(s in 1:nrow(locs)){
par(mar=c(3,3,3,6))
cc<-  dB[s,]> -3   # Cut off to hear a sound (due to background noise)
z<- dB[s,]
z[!cc]<- 0
if(sum(z) ==0) next
spatial.plot(gr,z,cx=6,col=terrain.colors(10))
##browser()
}
}

# A key feature of acoustic sampling is that there is a truncation at very low volumes due 
# to background ambient noise. Here I truncate at -3 dB

obs.dB<- dB
# truncation due to ambient noise.  A value of 0 = no discernable signal
obs.dB[obs.dB < -3]<- 0

# How many tweets did we pick up?
n.obs<- apply(obs.dB<0,1,sum)

# Matrix of observations
obs.dB<-obs.dB[n.obs>0,]
 
traps<- gr
# Run my mcmc function

 
out<- mcmc.fn(traps,obs.dB,xlim,ylim,900,100,cluster= TRUE)
simout[iter,]<- c(apply(out$parms,2,mean),N.signals)



}


s<- out$Sout
id<- out$ID
z<- out$zout

niter<- nrow(z)

nclust<- meanclust<- rep(NA, niter)
for(i in 1:niter){
 nclust[i]<- length(unique(id[i,][z[i,]==1]))
 meanclust[i]<-  mean(table(id[i,][z[i,]==1]))
}

stmp<- matrix(NA,nrow=niter,ncol=2)
for(i in 1:niter){
 stmp[i,1:2]<-  s[i,id[i,1],1:2]
}
plot(gr)
points(gr[obs.dB[1,]<0,],pch=20)
points(stmp,pch=20,col="red")
points(gr[obs.dB[1,]<0,],pch=20)


M<- ncol(id)
c1<- NULL
for(i in 1:niter){
 c1<-c(c1, (1:M)[id[i,]==1])
 
}












