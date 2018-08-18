
## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(moiR)
load_all('.')

###############################################################################
### With some real data

## Best positions
gwares<-readRDS('dataint/rFitness_mlp.rda')
gwares$order<-1:nrow(gwares)
maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))

subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[1000]),]
head(subg)
subg$rs<-paste(subg$chr,subg$pos,sep='_')

###############################################################################
## Subset of matrix
#data("genomes")

genomes<-readRDS('databig/genomes.rda')
G <- attachgenomes(genomes)
G[1:5,1:5]
Map <- genomes$map
head(Map)
Map$rs<-paste(Map$chr,Map$physical.pos,sep='_')
Map$order<-1:nrow(Map)

matched<-merge(Map, subg, by='rs')
dim(matched)
head(matched)

X=change01(inputena.mat(G[,matched$order],value = 0))
X[1:5,1:5]
# X<-inputena.mat(G[,sample(1:1000,size = 5)]) # this for random

################################################################################
## Get phenotypes
Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Fitness_mlp")], by.y="id",all.x=T)[,7]
Yabs<-Y
Y=relative(Y)
# Y=meanvarcent(Y)
Y[is.na(Y)]<-1

################################################################################
## Calculate starting values of S
s<-ridgegwa(X,Y,type = 'penalized')
hist(s)

hist(Y)

################################################################################
## Calculate observed LD change -> fast
# sourceCpp('src/ld.cpp')
# LDobs<-ldCnext(X,Y) ## much slower implementation
LDobs<-ldCnextmat(X,Y)
LDobs<-cleanld(LDobs)

summary(LDobs)

# LDobs[1:5,1:5]
# LDobs2[1:5,1:5]
# LDobs2[LDobs2>10]<-10
# LDobs[995:999,995:999]
# LDobs[1:5,1:5] * LDobs[1:5,1:5]
# hist(LDobs2,xlim = c(0,10),breaks=50)
# summary(fn(LDobs))
# LDobs[1,]
# LDobs[,1]
# View(LDobs)


################################################################################
## Calculate example of expected LD change that will be used in the object fun
sourceCpp('src/ld.cpp')
LDexpect<-ldCexpect(s,e=1)
hist(LDexpect,xlim=c(0,10),breaks=50)
LDexpect[LDexpect>10]<-10
LDexpect[LDexpect< (-10)]<- (-10)


################################################################################
### Object function
# Example of the starting likelihood
theta=c(1,s)
y=LDobs
ld.lik(theta, LDobs)

################################################################################
# With BFGS algorithm
res<-optim(par = c(1,s),
      fn = ld.lik,
      y=LDobs,
      method="BFGS")

# With
library(maxLik)
theta.start=c(1,s)
names(theta.start) = c("e",paste0("s",1:length(s)) )

theta.mle = maxLik(ld.lik, start=theta.start, x=LDobs )

save.image()

res$par[1]
hist(res$par[-1])

plot(res$par[-1],s)
cor.test(res$par[-1],s)


p1<-qplot(res$par[-1],s,xlab='Likelihood inference of s',ylab='Starting values from cGWA')
p1
save_plot('lik_s_inference.pdf',plot = p1)



# ldcreation(rnorm(10,1,1),rnorm(10,1,1),ep=rnorm(10,1,0.5))
#
# n=1e6
# hist(ldcreation(rnorm(n,1,1),rnorm(n,1,1),ep=rnorm(n,1,0.1)))
#
#
# ep=5
#
# s1.<-seq(-0.8,+0.8,by=0.05)
# s2.<-seq(-0.8,+0.8,by=0.05)
# ss<-expand.grid(s1.,s2.)
# s1=ss[,1]
# s2=ss[,2]
#
# ep=seq(0,3,by=0.1)
# abline(v = 1)
#
# rexp(100,0.1) %>% hist
#
#
# ldcreation(s1,s2,ep)
#
#
# ldcreation(0.2,0.5,ep=3)
# ldcreation(0.2,0.5,ep=1)
#
# s1=0.2
# s2=0.5
# ((1+s1)*(1+s2))^ep * ((1-s1)*(1-s2))^ep
#
# plot(ldcreation(0.2,0.5,ep=seq(0,+5,by=1)),
#     x= seq(0,+5,by=1))



