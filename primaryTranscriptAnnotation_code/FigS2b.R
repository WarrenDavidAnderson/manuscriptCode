

# exponential model for selecting fraction of gene end for peak search
add = 100000
tau = 50000
tmin = 0.3
tmax = 1
x = c(-add:0)
y0 = exp(-x/tau)
y = ( (y0-min(y0))/(max(y0)-min(y0)) ) * (tmax-tmin) + tmin

# selection at x = tau
xt = tau - add
yt = ( (exp(-xt/tau)-min(y0))/(max(y0)-min(y0)) ) * (tmax-tmin) + tmin


# examples - set clip
clip = add - tau

# code model
arg = -clip
y.model = ( (exp(-arg/tau)-min(y0))/(max(y0)-min(y0)) ) * (tmax-tmin) + tmin

# code format
ind = add - clip + 1
y.code = y[ ind ]

# simplified model1
arg = clip
y.simp1 = ( (exp(arg/tau)-min(y0))/(max(y0)-min(y0)) ) * (tmax-tmin) + tmin

# simplified model2
cp = c(0:add)
arg = clip
yr0 = exp(cp/tau)
y.simp2 = ( (exp(arg/tau)-min(yr0))/(max(yr0)-min(yr0)) ) * (tmax-tmin) + tmin

# show information
y.model
y.code
y.simp1
y.simp2

# Fig S2b
cp = c(0:add)
yr0 = exp(cp/tau)
yr = ( (yr0-min(yr0))/(max(yr0)-min(yr0)) ) * (tmax-tmin) + tmin
plot(cp,yr,type="l",xlab="bp clipped",ylab="fraction gene end")
abline(v=clip, lty=2); abline(h=y.simp2, lty=2)






