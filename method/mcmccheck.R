## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(coda)
library(Rcpp)

library(moiR)
load_all('.')


###############################################################################
## Genomematrix
################################################################################

genomes<-readRDS('../gws/databig/genome.rda')
# G <- attachgenomes(genomes)
Go<-bigmemory::attach.big.matrix('../gws/databig/gc.desc')


Go[1:10,1:10]

m=10

topcols<-SNPs<- SNPsub<-sample(1:ncol(Go),m)

newtopcols=topcols
################################################################################
## Fitness mlp
################################################################################
## In file data-cleaning/Yfitness.R
Y<-readRDS('data/Y.rds')


################################################################################
## Test functions
################################################################################
sourceCpp('MCMC.cpp')
# sourceCpp('MCMCmean.cpp')
# sourceCpp('MCMCmean_c.cpp')
# sourceCpp('MCMC2.cpp')


test_wfn()
test_Ey_go()
test_Prior(mode=3)
test_Prior(mode=3,variance = 0.5)

test_GProposal()

## Test selection proposals
test_Proposals()
test_Proposals(m=10,
               bw=5,
               min=0,
               max=+2.5,
               mode = 2,
               verbose=TRUE
                )

# R2=LDrelative(Gdouble@address,topcols-1,debug=TRUE)
# test_ProposalsLD(R2,verbose = TRUE,bw = 1)

LLGaussMix(
          y=0,
          e = .1,
          v = .1 * 1,
          p = .8
  )

LLGaussMix(
          y=0.01,
          e = .1,
          v = .05,
          p = .1
  )

test_Likelihood(
          A=Go@address,
          y=Y$oFitness,
          h=Y$row,
          s=runif(m,-1,+1),
          b=1.5,
          a=0.5,
          p=0.8,
          # mu=0.8,
          n=1:nrow(Go)-1,
          m=topcols-1,
          Fitnessmode=1,
          verbose=TRUE)


test_Prior(m,0,1,0,0.5,mode=3)


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
# ## MCMC check that Posterior = Prior when likelihood is 1
# ################################################################################

# # sourceCpp('MCMC.cpp')
# sourceCpp('MCMCmean.cpp')

r0<-napMCMC(
          y=Y$Fitness,
          h=Y$row,
          A=Go,
          s=PropoS(length(SNPsub),0.2),
          m=SNPsub-1,
          n=1:nrow(Go),
          mu=100,
          iterations = 1e3,
          TEST = TRUE, # DO THE TEST
          verbose=FALSE,
          debug=TRUE,
          Fitnessmode=3,
          Priormode=1, # need for test also
          # Priormode=3, # need for test also
          Proposalmode=3
    )
tail(r0$parchain)

parchain<-as.mcmc(r0$parchain)
colnames(parchain)<- r0$parnames
plot(parchain)
plot(parchain[,7])

schain<-as.mcmc(r0$chain)
head(schain)
tail(schain)
plot(schain)
plot(log(schain+1))

# pdf('figs/MCMC-TEST_lik=1.pdf')
# plot(r0)
# dev.off()





################################################################################
## MCMC minimal example
################################################################################

# sourceCpp('MCMC.cpp')
# sourceCpp('MCMCmean.cpp')
set.seed(1)
# pdf('figs/MCMC-TEST_simpleexample_many.pdf',height = 4,width = 4,useDingbats = FALSE)

## model of fitness
modelfitness=2
modelfitness=3
modelfitness=1
## Parameters
p=0.2
b=0.2
a=0.1
zeroprop=0.5
mu=1

## Select SNPs
# m=3
# m=10
m=100
# m=1000
# m=100000

newtopcols<-topcolsall[1:m]
s = rexp(m,0.1); s=s/max(s);
s= s* rbinom(n = length(s),size = 1,p=1-zeroprop)
if(modelfitness != 3) s= s* sample(c(-1,1),size = m,replace = TRUE)
qplot(s)
s
Go[,newtopcols]
Go[,newtopcols] %>% cor
Go[,newtopcols] %>% cor %>% upperTmat() %>% max
## Create the R2 matrix for LD informed proposals

# R2=LDrelative(Gdouble@address,newtopcols-1)

## Create expectations
Eys= Ey_go(Go[,newtopcols],s = s,mode=modelfitness)
qplot(Eys,xlab="True genotype fitness")

## with other models
qplot( Ey_go(Go[,newtopcols],s = s,mode=1),main=paste('mode = ',1),
       xlab="True genotype fitness")
qplot( Ey_go(Go[,newtopcols],s = s,mode=2),main=paste('mode = ',2),
       xlab="True genotype fitness")
qplot( Ey_go(Go[,newtopcols],s = s,mode=3),main=paste('mode = ',3),
       xlab="True genotype fitness")

## Sample 5 repliates from those expectations
Yobs<- sapply(1:nrow(Go), function(i) rnorm(5,Eys[i],Eys[i]*b))
Yobs<- sapply(1:nrow(Go), function(i) rnorm(5,Eys[i],a+Eys[i]*b))
Yobs[Yobs<0] <-0

Yobs<-c(Yobs)
Yobs[sample(1:length(Yobs),size = floor(p*length(Yobs)) )]<-0

# The ids of genotypes
Hids<-sort(rep(1:nrow(Go), 5))

## Buid data
newY<-data.frame(rFitness=Yobs,row=Hids)
qplot(newY$rFitness,xlab="Observed fitness (all replicates)")

################################################################################
## Run MCMC

sstart<-runif(m,0,1)
sstart<-runif(m,-1,1)
# sourceCpp('MCMC.cpp')
n=length(unique(newY$row))

test_Likelihood(y=newY$rFitness,
                h=newY$row,
                A=Go@address,
                s=s,
                m = newtopcols-1,
                n= unique(newY$row)-1,
                b = .5,
                a=0.1,
                p=.1,
                Fitnessmode = 1,
                TEST=FALSE,
                verbose = FALSE)
an<-gwsMCMC(
          y=newY$rFitness,
          h=newY$row,
          A=Go@address,
          s=rep(0,m),
          m=newtopcols-1, # important -1
          n= unique(newY$row)-1,
          bw=.1,
          b = .5, bmin=0,bmax=1,
          a=0.1, amin=0,amax=0.5,
          p = 0.1,
          mu=1,mumin=1,mumax=1,
          epi=1,epimin=1,epimax=1,
          nupdates=1,

          min= ifelse(modelfitness==3, 0, -1), # in inverse multiplicative is 0 minimum
          max=ifelse(modelfitness==3, 10, 1),

          # iterations = 10000,
          iterations = 10000,
          TEST = FALSE,
          verbose=FALSE,
          debug=FALSE,
          Fitnessmode=modelfitness,
          Priormode=1,
          Proposalmode=1 ## not LD informed
        )


# Global parameters
# qplot(sample(an$p,1000),geom="density",xlab="p of Bernoulli mortality", xlim=c(0,1),fill="grey")+
#   geom_vline(xintercept = PI, lty="dashed")
# qplot(sample(an$b,1000),geom="density",xlab="mean to sd ratio", xlim=c(0,1),fill="grey")+
#   geom_vline(xintercept = BE, lty="dashed")
# qplot(
#   # sample(an$a,1000),
#   tail(an$a,1000),
#   geom="density",xlab="mean to sd ratio", xlim=c(0,1),fill="grey")+
#   geom_vline(xintercept = AI, lty="dashed")

parchain<-as.mcmc(an$parchain)
colnames(parchain)<- an$parnames
plot(parchain)

schain<-as.mcmc(an$chain)
plot(schain[,1:10])


# Estimate of selection coefficients
schain=as.mcmc(an$chain)
sinf<-apply(schain,2,median)
qplot(sinf)
sinf_range<-HPDinterval(schain)
sinf

psel<-qplot(x=s,y=sinf,
      # xlim=c(0,1),
      # ylim=c(0,1)
      xlim=c(-1,1),
      ylim=c(-1,1)
      ) +
  geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )+
  geom_abline(intercept = 0,slope = 1,lty="dotted")+
  ylab("Inferred selection coefficients")+
  xlab("True selection coefficients")
psel
cor.test(s,sinf)

toplot<-data.frame(y= Ey_go(X=Go[,newtopcols],s = sinf,mode=modelfitness,mu = tail(an$parchain[,4],1) ),
      x=Eys
      )
toplot$xnozero<-toplot$x
toplot$xnozero[toplot$xnozero==0]<-NA
pind<-ggplot(data=toplot,aes(y=y,x=x)) +
geom_point(col="grey") +
stat_smooth(aes(y=y , x=xnozero), se=F,
            method="glm",lty="dashed",col="black")+
      ylim(c(0,max(c(toplot$y,toplot$x)) ) ) +
      xlim(c(0,max(c(toplot$y,toplot$x))))+
      ylab("Inferred individual fitness")+
      xlab("True individual fitness")+
  geom_abline(intercept = 0,slope = 1,lty="dotted")
pind

infhist<-qplot(Ey_go(X=Go[,newtopcols],s = sinf,mode=modelfitness,mu = tail(an$parchain[,4],1) ),xlab="")+
  xlim(c(0,max(c(toplot$y,toplot$x))))+coord_flip()
truehist<-qplot(Eys)+xlim(c(0,max(c(toplot$y,toplot$x))))

panel<-plot_grid(
          truehist,
          ggplot(),
          pind,
          infhist ,
          ncol=2,nrow=2,
          rel_widths =c(0.5,0.1),rel_heights = c(0.1,0.5)
        )


panel2<-plot_grid(psel,panel)
panel2
# dev.off()

pdf(paste0('trialMCMC_',runname,".pdf"), width = 20,height = 10)
  print(panel2)
dev.off()


################################################################################
#### TRIAL GWA ####
sourceCpp("gwa.cpp")
Gdouble<-bigmemory::attach.big.matrix('databig/godou.desc')
Gdouble[1:10,1:10]
Ymean<-Vy(newY$rFitness,newY$row)


betam= - BMmgwa(Gdouble@address,Ymean-1,newtopcols-1) # Not fair due to distribution
betac= - BMcgwa(Gdouble@address,Ymean-1,newtopcols-1) # Not fair due to distribution
# betam_log= BMmgwa(Gdouble@address,log(2-Ymean),newtopcols-1)
# betam=exp(betam_log) -.001

pdf('figs/MCMC-TEST_compareGWAc.pdf',height = 4,width = 4,useDingbats = FALSE)

qplot(x=s,y=betac,xlim=c(0,1),ylim=c(0,1)) + geom_abline(intercept = 0,slope = 1,lty="dotted")+
  ylab("Inferred selection coefficients")+
  xlab("True selection coefficients")+
  ggtitle("GWA")

X=Go[,newtopcols]

toplot<-data.frame(y= X %*% betac,
                    x=newY$rFitness)
toplot$xnozero<-toplot$x
toplot$xnozero[toplot$xnozero==0]<-NA
ggplot(data=toplot,aes(y=y,x=x)) +
ggtitle("GWA")+
geom_point(col="grey") +
stat_smooth(aes(y=y , x=xnozero), se=F,lty="dashed",col="black")+
      # ylim(c(0,1.2)) +
      xlim(c(0,1.2))+
      ylab("Inferred individual fitness")+
      xlab("True individual fitness")
dev.off()


########################################################
########################################################
########################################################
########################################################
########################################################

# burnin = floor(nrow(an$chain)*0.1)
# rmoc<-an$chain[-c(1:burnin),]
# rmoc<-an$chain
#
# mcmcres<-plot_grid(
#   qplot(fn(rmoc[,1]), xlim=c(-.1,1.1), xlab='s',main='SNP1') +geom_vline(xintercept = s[1], col='red',lty='dashed'),
#   qplot(fn(rmoc[,2]), xlim=c(-.1,1.1), xlab='s',main='SNP2') +geom_vline(xintercept = s[2], col='red',lty='dashed'),
#
#   qplot(fn(rmoc[,1]),x=1:nrow(rmoc), ylim=c(-.1,1.1),geom='line',xlab='iterations',ylab='s') +
#     # geom_hline(yintercept = sstart[1], col='grey',lty='dashed') +
#     geom_hline(yintercept = s[1], col='red',lty='dashed') ,
#   qplot(fn(rmoc[,2]),x=1:nrow(rmoc), ylim=c(-.1,1.1),geom='line',xlab='iterations',ylab='s') +
#     # geom_hline(yintercept = sstart[2], col='grey',lty='dashed')+
#     geom_hline(yintercept = s[2], col='red',lty='dashed')
# )
# mcmcres
# save_plot(filename = 'figs/MCMC-TEST_simpleexample.pdf',plot = mcmcres,base_height = 7,base_width = 7)

# pdf('figs/MCMC-TEST_simpleexample_many.pdf',height = 7,width = 7,useDingbats = FALSE)
# for(i in 1:m){
#   tmpp<-plot_grid(
#   qplot(fn(rmoc[,i]), xlim=c(-.1,1.1), xlab='s',main=paste('SNP',i)) +
#   geom_vline(xintercept = s[i], col='red',lty='dashed'),
#     qplot(fn(rmoc[,i]),x=1:nrow(rmoc), ylim=c(-.1,1.1),geom='line',xlab='iterations',ylab='s') +
#       geom_hline(yintercept = s[i], col='red',lty='dashed')
#                   )
#   print(tmpp)
# }
# dev.off()
# dev.off()


# test_Likelihood(y=newY$oFitness,
#                 h=newY$row,
#                 A=Go@address,
#                 s=s,
#                 m = topcols-1,
#                 # be = .1,
#                 be = 1,
#                 p=.01,
#                 Fitnessmode = modelfintes,
#                 TEST=FALSE,
#                 verbose = TRUE)

## Check, is the likelihood better or worse?

# Eys_inferred= Ey_go(Go[,topcols],s = sinf,mode=modelfintes)
# Eys_start= Ey_go(Go[,topcols],s = sstart,mode=modelfintes)
#
#
# hist(Eys,xlim=c(0,1))
# # hist(Eys_start,xlim=c(0,1))
# hist(Eys_inferred,xlim=c(0,1))
# hist(newY$oFitness,xlim=c(0,1.1))
#
# plot(newY$oFitness~Eys[newY$row])
# plot(newY$oFitness~Eys_inferred[newY$row])
# cor.test(newY$oFitness,Eys[newY$row])
# cor.test(newY$oFitness,Eys_inferred[newY$row])

# plot(newY$oFitness~Eys_start[newY$row])

#
# test_Likelihood(y=newY$oFitness,
#                 h=newY$row,
#                 A=Go@address,
#                 s=s,
#                 m = topcols-1,
#                 be = .1,
#                 pi=.01,
#                 Fitnessmode = modelfintes,
#                 TEST=FALSE,
#                 verbose = TRUE)
# test_Likelihood(y=newY$oFitness,
#                 h=newY$row,
#                 A=Go@address,
#                 s=sinf,
#                 m = topcols-1,
#                 be = .1,
#                 pi=.01,
#                 Fitnessmode = modelfintes,
#                 TEST=FALSE,
#                 verbose = TRUE)
#
#
#
# sum(sapply( 1:nrow(newY) ,
#                   function(i) LLGaussMix(newY$oFitness[i],Eys[newY$row][i],Eys[newY$row][i]*.1,pi=0.9 ) ) ,na.rm = TRUE)
# sum(sapply( 1:nrow(newY) ,
#                   function(i) LLGaussMix(newY$oFitness[i],Eys_inferred[newY$row][i],Eys_infered[newY$row][i]*.1,pi=0.9 )) ,na.rm = TRUE)
# sum(sapply( 1:nrow(newY) ,
#                   function(i) LLGaussMix(newY$oFitness[i],Eys_start[newY$row][i],Eys_start[newY$row][i]*.1,pi=0.9 )) ,na.rm = TRUE)
#
#
# 0.1 * 1.5*1.9

# newY$oFitness[i]
# i=1
# i=100
# i=sample(1:nrow(newY),1)
# newY$oFitness[i]
# alpha(Eys_start[newY$row][i],M2V(poly,Eys_start[newY$row][i]))
# beta(Eys_start[newY$row][i],M2V(poly,Eys_start[newY$row][i]))
#
# i=sample(1:nrow(newY),1)
# LLGammaR(
#   newY$oFitness[i],
#   Eys_infered[newY$row][i],
#   M2V(poly,Eys_infered[newY$row][i])
#          )
#
#
# i=sample(1:nrow(newY),1)
# LLGammaR(
#   newY$oFitness[i],
#   Eys[newY$row][i],
#   M2V(poly,Eys[newY$row][i])
#          )
# dgamma(
#   x= newY$oFitness[i],
#   shape=alpha(Eys[newY$row][i],M2V(poly,Eys[newY$row][i])),
#   rate=beta(Eys[newY$row][i],M2V(poly,Eys[newY$row][i])),
#   log=FALSE
# )
#
# dgamma(
#   x=0,
#   shape=alpha(Eys[newY$row][i],M2V(poly,Eys[newY$row][i])),
#   rate=beta(Eys[newY$row][i],M2V(poly,Eys[newY$row][i])),
#   log=FALSE
# )
#
# dnorm(1,1.0000002,0.1)
# rnorm(1000000,1.0000002,0.1) %>% hist(freq=FALSE)
# rnorm(1000000,2,0.2) %>% hist(freq=FALSE)
#
# dnorm(1,1.0000002,0.1,log=TRUE)
#
# dnorm(0,0,1)

###############

# Real<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys[newY$row][i],M2V(poly,Eys[newY$row][i]) ) ) ,na.rm = TRUE)
# Real
# Infe<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_infered[newY$row][i],M2V(poly,Eys_infered[newY$row][i])) ) ,na.rm = TRUE)
# Infe
# Start<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_start[newY$row][i],M2V(poly,Eys_start[newY$row][i])) ) ,na.rm = TRUE)
# Start
#
# Real<-prod(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys[newY$row][i],0.1 ,FALSE) ) ,na.rm = TRUE)
# Real
# Infe<-prod(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_infered[newY$row][i],0.1, FALSE ) ) ,na.rm = TRUE)
# Infe
# Start<-prod(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_start[newY$row][i],0.1 , FALSE ) ),na.rm = TRUE)
# Start


# Real<-prod(sapply( 1:nrow(newY) , function(i) dgamma(newY$oFitness[i],Eys[newY$row][i],0.1 ,FALSE) ) ,na.rm = TRUE)
# Real
# Infe<-prod(sapply( 1:nrow(newY) , function(i) dgamma(newY$oFitness[i],Eys_infered[newY$row][i],0.1, FALSE ) ) ,na.rm = TRUE)
# Infe
# Start<-prod(sapply( 1:nrow(newY) , function(i) dgamma(newY$oFitness[i],Eys_start[newY$row][i],0.1 ), FALSE ) ,na.rm = TRUE)
# Start


#
#
# Eys_start= Ey_go(Go[,topcols],s = runif(2),mode=modelfintes)
# sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_start[Y$row][i],M2V(poly,Eys_start[Y$row][i])) ) ,na.rm = TRUE)
#
#
# sourceCpp('MCMC.cpp')
# Real<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys[Y$row][i],.01) ) ,na.rm = TRUE)
# Real2<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys[Y$row][i],Vys[Y$row][i] ) ) ,na.rm = TRUE)
# Infe<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_infered[Y$row][i],.01 )) ,na.rm = TRUE)
# Start<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_start[Y$row][i],.01 )) ,na.rm = TRUE)
# Real
# Real2
# Infe
# Start
#
#
#
# Real<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys[Y$row][i],M2V(poly,Eys[Y$row][i]) ) ) ,na.rm = TRUE)
# Infe<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_infered[Y$row][i],M2V(poly,Eys_infered[Y$row][i])) ) ,na.rm = TRUE)
# Start<-sum(sapply( 1:nrow(newY) , function(i) LLGammaR(newY$oFitness[i],Eys_start[Y$row][i],M2V(poly,Eys_start[Y$row][i])) ) ,na.rm = TRUE)
# Real
# Infe
# Start

# for(v in seq(0.01,1,length.out = 10)){
#   Infe<-sum(sapply( 1:nrow(newY) , function(i) LLGammap(newY$oFitness[i],Eys_infered[Y$row][i],v) ) ,na.rm = TRUE)
#   Start<-sum(sapply( 1:nrow(newY) , function(i) LLGammap(newY$oFitness[i],Eys_start[Y$row][i],v) ) ,na.rm = TRUE)
#   Real<-sum(sapply( 1:nrow(newY) , function(i) LLGammap(newY$oFitness[i],Eys[Y$row][i],v) ) ,na.rm = TRUE)
#   # print(Infe<Start)
#   message(Infe," ",Start, " ", Real)
# }
#
# LLGamma(.9,1,0.3,calllgamma(1,0.3)) %>% exp
#
# LLGammap(.91,1,0.3) %>% exp
# rgamma(100,1,0.3) %>% hist
# pgamma(1,1,0.3,lower.tail = FALSE)
# pgamma(1,1,0.3,lower.tail =TRUE)
# pgamma(15,1,0.3,lower.tail = FALSE)
# # pgamma(15,1,0.3,lower.tail =TRUE)
#
# dgamma(1,1,0.3)
#
# LLGammap(.5,1,0.3) %>% exp
#
# xLLGammap(1,1,0.2,calllgamma(1,0.2))
#
# LLGammap(1,1,0.2,calllgamma(1,0.2))
# LLGammap(1,1,0.2) %>% exp
#
# exp( LLGammap(1,1,0.2,calllgamma(1,0.2)) )
#
# LLGamma(yh=newY$oFitness)
# LLik

################################################################################
## Various fitness models
################################################################################
# Rcpp::sourceCpp('matrixaccess.cpp')
# Rcpp::sourceCpp('MCMC.cpp')

# Checking genotypes all below 1 and
# apply(X,1,sum)
# apply(X,1,sum) %>% hist

# s= runif(10,0,1)
# Ey(X,s) %>% hist
# Ey_additive(X,s) %>% hist
# Ey_inversemult(X,s) %>% hist
# Ey_inversemult_center(X,s) %>% hist
# Ey_inversemult_center(Xc,s) %>% hist

# hist(momentrgamma(50))



################################################################################
## Dealing with fixed variances
################################################################################

# sourceCpp('MCMC.cpp')
# pmeanvar<-ggdotscolor(Vy(Y$oFitness,Y$row),My(Y$oFitness,Y$row),
#             xlab='Fitness variances',
#             ylab='Fitness means') %>%
#   addggregression(formula=y~poly(x,2),se=FALSE,lty='dashed',col='darkgrey') +
#   xlim(c(0,0.3))+
#   ylim(c(0,0.5))
# save_plot(file='figs/fitness_mean_variance_pergeno.pdf',pmeanvar,
#   base_height=6,base_width=6)

# datamod<-data.frame(variance=Vy(Y$oFitness,Y$row),mean=My(Y$oFitness,Y$row))
# lmod<-lm(data=dplyr::filter(datamod,variance>0), variance~poly(mean,2))
# summary(lmod)
# qplot(y=predict(lmod) ,x= dplyr::filter(datamod,mean>0)$mean)
# abs(predict(lmod,newdata = dplyr::filter(datamod,variance<0)))

# sourceCpp('MCMC.cpp')
# goodvars<-Vyimpute(Y$oFitness, Y$row)
# goodvars

# qplot(y=goodvars ,x= My(Y$oFitness,Y$row))

#>>>>>> There are still variances that are lower than 0, I just make them positive

################################################################################
## Microbenchmark of class fitness function implementation
################################################################################
# sourceCpp('MCMC.cpp')
# library(microbenchmark)

# m=1000
# n=10
# s=runif(m,0,1)
# microbenchmark(Ey_additive(Go[1:n,1:m],s),
#                Ey_general(Go[1:n,1:m],s,1),
#                Ey_go(Go[1:n,1:m],s,1),
#                times=100
#                )

# m=ncol(Go)
# m=10
# n=nrow(Go)
# s=runif(m,0,1)


#### The different implementations
# >> Expectation of fitness for the whole genome takes about 1 minute
# microbenchmark(Ey_additive(Go[1:n,1:m],s),
#                Ey_general(Go[1:n,1:m],s,1),
#                Ey_go(Go[1:n,1:m],s,1),
#                times=10
#                )

#### Matrix vs pointer implementation
# sink('tables/expectation_benchmark_accessor.txt')
# timestamp()
# microbenchmark(times=1,
#   Ey_go(Go[],runif(ncol(Go),0,1),mode=1),
#   Ey_go_ma(Go@address,runif(ncol(Go),0,1),mode=1)
#   )
# microbenchmark(times=10,
#     Ey_go_ma(Go@address,runif(ncol(Go),0,1),mode=1)
# )
# sink()

##>>>>>>>>>> Discovered that classes are more efficient even having switch statemetns

################################################################################
## Microbenchmark of likelihood calculation parts
################################################################################

# sourceCpp('MCMC.cpp')
# m=100000
# m=ncol(Go)
# n=nrow(Go)
# s=runif(m,0,1)

#### Increase of time with number of SNPs
# n=nrow(Go)
# sink('tables/likelihood_time_benchmark.txt')
# timestamp()
# microbenchmark(times=10,
#   snps10=LLik(Y$oFitness,Y$row,Go[1:n,1:10],runif(10,0,1),Myvar,debug = FALSE,mode= 1),
#   snps50=LLik(Y$oFitness,Y$row,Go[1:n,1:50],runif(50,0,1),Myvar,debug = FALSE,mode= 1),
#   snps100=LLik(Y$oFitness,Y$row,Go[1:n,1:100],runif(100,0,1),Myvar,debug = FALSE,mode= 1),
#   snps1000=LLik(Y$oFitness,Y$row,Go[1:n,1:1000],runif(1000,0,1),Myvar,debug = FALSE,mode= 1),
#   snps10000=LLik(Y$oFitness,Y$row,Go[1:n,1:10000],runif(10000,0,1),Myvar,debug = FALSE,mode= 1),
#   snpsall=LLik(Y$oFitness,Y$row,Go[1:n,1:ncol(Go)],runif(ncol(Go),0,1),Myvar,debug = FALSE,mode= 1)
#   )
# sink()

#### Evaluation of likelihood genome-wide
# sink('tables/whole-genome_expectation_likelihood_benchmark.txt')
# cat(timestamp())
# microbenchmark(times=10,
#   Ey_go_ma(Go@address,runif(ncol(Go),0,1),mode= 1),
#   goLLik(Y$oFitness,Y$row,Go@address,runif(ncol(Go),0,1),Myvar,debug = FALSE,mode= 1)
# )
# sink()

#>>>>>>>>> Discovered that passing the memory address is much more efficient.

