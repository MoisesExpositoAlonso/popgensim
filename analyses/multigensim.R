sourceCpp('multitools.cpp')

############ funcitons
source('r/multigensim.R')

############ starting parameters

# Example
p=1000
n=100
No=5000
Nmax=1e6
MU=5
maf=mafsim(p,'exponential',1)
X=Xsim(n,p,maf)
X[1:5,1:5]
s=ssim(p,svar = 0.1)
w=wsim(X,s)
d=0.5
rate=1
svar=0.1
tmax=5

############ simulations

pop<-multigenpopsimC(w,No = 100000,p = 0.3)
popline<-apply(pop,2,sum)
tmax=5



hist(wsim(X,s,mode=1,mu=1))
hist(wsim(X,s,mode=2,mu=1))
hist(wsim(X,s,mode=3,mu=1))

multigenpopsimC(w)


multievo(p,n,rate,svar,mode,No,Nmax,d,tmax)
