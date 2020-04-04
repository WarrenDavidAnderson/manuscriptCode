
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure S2
############################################################################


dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/figS2"
setwd(dir)

#################################################
# create simulation
#################################################

library(pracma)

# function to compute inverse of cdf
inv.norm = function(rank=NULL,n=NULL){
  qnorm( rank/(n+1) )
}

# function to get observed and expected fold changes
get.obs.exp = function(n=NULL){
  # ndat = rnorm(n,0,1)
  # ndat = ndat[order(ndat)] 
  nm = c(1:n)
  muf = rep(0,length(nm))
  mum = rep(0,length(nm))
  x = rep(0,length(nm))
  y = rep(0,length(nm))
  etd = rep(0,length(nm))
  for(ii in nm){
    ndat = rnorm(n,0,1)
    ndat = ndat[order(ndat)]
    mum[ii] = mean(ndat[1:ii])
    muf[ii] = mean(ndat[(ii+1):n])
    arg = -sqrt(2)*erfinv(1-2*ii/(1+n))
    x[ii] = (1/sqrt(2*pi))*exp((-1/2)*(arg)^2)
    y[ii] = ii/(n+1)
    etd[ii] = x[ii] / (y[ii]*(1-y[ii]))
  }
  udat = data.frame(fm=nm/n, mum=mum, muf=muf, 
                    obs=muf-mum, exp=etd, 
                    numerator=x, 
                    denominator=y*(1-y))
  udat$obs[n] = udat$obs[n-1]
  return(udat)
} # get.obs.exp

# plot function
plt.fn = function(dat=NULL, col="red", fname=NULL){
  plot(dat$fm, dat$obs, pch=19,xlab="fraction males",ylab="max log fold change",ylim=c(1,5))
  lines(dat$fm, dat$exp, type="l", col=col, lwd=2)
  abline(h=4/sqrt(2*pi), lty=2, col="red", lwd=2)
}

# loop through various sample sizes
nsamp = seq(100,100000,by=100)
reso = rep(0,length(nsamp))
rese = rep(0,length(nsamp))
for(ii in 1:length(nsamp)){
  n = nsamp[ii]
  nm = 1
  ndat = rnorm(n,0,1)
  ndat = ndat[order(ndat)]
  mum = mean(ndat[1:nm])
  muf = mean(ndat[(nm+1):n])
  arg = -sqrt(2)*erfinv(1-2*nm/(1+n))
  x = (1/sqrt(2*pi))*exp((-1/2)*(arg)^2)
  y = nm/(n+1)
  etd = x / (y*(1-y))  
  reso[ii] = muf-mum
  rese[ii] = etd
}
nobs = data.frame(nsamp=nsamp, maxfcO=reso, maxfcE=rese)

# plot results
dat1 = get.obs.exp(n=100)
dat2 = get.obs.exp(n=1000)
dat3 = get.obs.exp(n=10000)

pdf("figS2.pdf")
par(mfrow=c(3,3))
plt.fn(dat1)
plt.fn(dat2)
plt.fn(dat3)
plot(nobs$nsamp,nobs$maxfcO,type="l",lwd=2,xlab="sample size",ylab="max log fold change",ylim=c(1,5))
lines(nobs$nsamp,nobs$maxfcE,type="l",col="red",lwd=2)
abline(h=4/sqrt(2*pi), lty=2, col="red", lwd=2)
dev.off()



