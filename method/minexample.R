


Eys= Ey_go(Go[,SNPs],s = PropoS(length(SNPs),0.1),mode=1)
hist(Eys)

test_Ey_go()

test_Likelihood(
          A=Go@address,
          y=Y$Fitness,
          h=Y$row,
          s=PropoS(length(SNPs),0.2),
          n=1:nrow(Go)-1,
          b=1,
          a=0.2,
          p=0.8,
          # mu=mean(Y$Fitness[Y$Fitness!=0]),
          m=topcols-1,
          Fitnessmode=2,
          verbose=TRUE)

LLGaussMix(
          y=130 /97,
          e = .1,
          v = .1 * 1,
          p = .8
  )

# ################################################################################
# ## MCMC doublecheck
# ################################################################################

sourceCpp('MCMC2.cpp')
sourceCpp('multitools.cpp')

r<-napMCMCR()

# ################################################################################
# ## MCMC with a minimum exmaple
# ################################################################################
n=500
m=1000
maf=mafsim(m)
X <- Xsim(n,m,maf)
# X=Go[1:n,1:m]
# apply(G,2,mean)
# apply(Go[,1:1000],2,mean)
# cor(sort(apply(G,2,mean)),sort(apply(Go[,1:1000],2,mean)))
# plot(sort(apply(G,2,mean)),sort(apply(Go[,1:1000],2,mean)))

svar=0.1
a=0
b=0
p=0
replicates=5

s= ssim(p,svar)

Ey=wsim(X,s,mode=2)
Fitness=sampleEys(Ey,a,b,p,rep = replicates)
# Fitness[is.na(Fitness)]<-0
h=sort(rep.int(1:n,replicates))
length(h)
length(Fitness)
plot(Fitness,Ey[h])


r0<-napMCMC(
          y=Fitness,
          h=h,
          A=X,
          s=rnorm(m,0,0.5),
          m=1:m,
          n=1:n,
          # mu=100,
          iterations = 5e4,
          TEST = TRUE, # DO THE TEST
          # verbose=FALSE,
          # debug=TRUE,
          Fitnessmode=2,
          Priormode=3,
          Proposalmode=3
    )
parchain<-as.mcmc(r0$parchain)
colnames(parchain)<- r0$parnames
plot(parchain)
plot(parchain[,7])


tail(r0$parchain)

schain=r0$chain
sinf<-MCMCglmm::posterior.mode(as.mcmc(schain))
sinf_range<-HPDinterval(as.mcmc(schain))
plot(sinf,s)




# schain<-as.mcmc(r0$chain)
# head(schain)
# tail(schain)
# plot(schain)
# plot(log(schain+1))
#

