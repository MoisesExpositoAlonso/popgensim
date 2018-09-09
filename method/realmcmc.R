## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(coda)

library(moiR)
load_all('.')
library(Rcpp)


###############################################################################
## Genomematrix
################################################################################
library(bigmemory)

genomes<-readRDS('../gws/databig/genome.rda')
# G <- attachgenomes(genomes)
#
# # Gi<-bigmemory::attach.big.matrix('databig/gi.desc')
# # Gm<-bigmemory::attach.big.matrix('databig/gm.desc')
# Go<-bigmemory::attach.big.matrix('databig/go.desc')
# Godouble<-bigmemory::attach.big.matrix('databig/godou.desc')
# Gi<-bigmemory::attach.big.matrix('databig/gi.desc')
# Gc<-bigmemory::attach.big.matrix('databig/gc.desc')

# G <- attachgenomes(genomes)

Go<-bigmemory::attach.big.matrix('../gws/databig/gc.desc')
load('../gws/data/map.rda')
head(map)

SNPsnames<-readRDS('../gws/dataint/topgamma_mlp.rda')
SNPs<-which(map$SNP %in% SNPsnames)
head(SNPs)

################################################################################
## Fitness models
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
## Fitness mlp
################################################################################
## In file data-cleaning/Yfitness.R
Y<-readRDS('data/Y.rds')

hist(Y$rFitness, breaks=1000)
hist(Y$rFitness, breaks=100)

################################################################################
## MCMC test
################################################################################

r0<-napMCMC(
            y=Y$Fitness,
            h=Y$row,
            A=Go,
            s=PropoS(2,0.5),
            b=0.3,bmin=0,bmax=1,
            a=0.1,amin =0, amax=5,
            mu=1, mumin=0.01,mumax=2,
            epi=1.01,epimin = 0, epimax=2,
            p=0.8,
            m=1:10,
            n=(1:nrow(Go)),
            iterations = 1e5,
            TEST = TRUE,
            Fitnessmode=1,
            Priormode=1,
            Proposalmode=1
            )
testrun<-as.mcmc(r0$parchain)
colnames(testrun)<-r0$parnames
plot(testrun)

pp<-parameterplot_(r0$parchain, r0$parnames)
pp
# save_plot(pp,filename = "figs/test_mcmc.pdf", ncol=2, nrow=7,base_width = 3.5,base_height = 1.8)



################################################################################
## MCMC real data
################################################################################
sourceCpp('MCMC.cpp')
# sourceCpp('MCMCmean_c.cpp') ## Trying out the scaling of phenotypes directly
sourceCpp('gwa.cpp')


SNPsub<- sort(sample(SNPs,100))
SNPsub<- SNPs


Fmode=2

r1<-napMCMC(
          y=Y$Fitness,
          h=Y$row,
          A=Go,
          s=PropoS( length(SNPsub), 0.1),
          b=1,bmin=0,bmax=1.5,
          a=.1,amin =0, amax=0.5,
          # mu=100,
          # mu=1,mumin=1,mumax=1,
          # mu=0.5,mumin=0.5,mumax=0.5,
          # mu=0.3,mumin=0.1,mumax=20,
          mu=97,mumin=1,mumax=120,
          epi=1, epimin=1,epimax=1,
          # epi=1.01, epimin=0.1,epimax=3,
          svar=0.2,svarmax = 0.25,
          ss=0.2,ssmin=0,ssmax=1,
          p=0.8,

          # m=SNPsub-1,
          # n=(unique(Y$row)) -1,
          m=SNPsub,
          n=(unique(Y$row)),

          iterations = 1e4,
          Fitnessmode=Fmode,
          Priormode=3,
          Proposalmode=3
)

### Parameters
# parameterplot_(r1$parchain, r1$parnames)
parinf<-parametersummary(r1$parchain, r1$parnames)
print(parinf)

### SNPs
sinf<-parametersummary(r1$chain , SNPsub)
# smode<-MCMCglmm::posterior.mode(r1$chain)
hist(sinf$statistics[,1])

### Individuals
# Ypred<-Ey_go(Go[,SNPsub], sinf$statistics[,1], Fmode, epi= parinf$statistics['epi',1] ) * parinf$statistics["mu",1]
# Ypred<-Ey_go(Go[,SNPsub], sinf$statistics[,1], Fmode)
# hist(Ypred, breaks=100)
# hist(Y$rFitness, breaks=100)
#
#
# toplot<-data.frame(x=Y$Fitness/ parinf$statistics["mu",1], y=Ypred[Y$row])
# toplot$xnozero <- toplot$x
# toplot$xnozero[toplot$xnozero==0] <-NA
# pnap<-ggplot(toplot, xlab='True individual fitness', ylab='Inferred fitness') + geom_point(aes(y=y,x=x), color='darkgrey') + stat_smooth(aes(y=y,x=xnozero), method='glm',se=F)
# print(pnap)


### Individuals
# Ypred<-Ey_go(Go[,SNPsub], sinf$statistics[,1], Fmode, epi= parinf$statistics['epi',1] ) * parinf$statistics["mu",1]
Ypred<-Ey_go(Go[,SNPsub], sinf$statistics[,1], Fmode) * parinf$statistics["mu",1]
hist(Ypred, breaks=100)
hist(Y$rFitness, breaks=100)

# toplot<-data.frame(x=Y$Fitness , y=Ypred[Y$row])
# toplot$xnozero <- toplot$x
# toplot$xnozero[toplot$xnozero==0] <-NA
# # pnap3<-ggplot(toplot, xlab='True individual fitness', ylab='Inferred fitness') + geom_point(aes(y=y,x=xnozero), color='darkgrey') %>% addggregression()# + stat_smooth(aes(y=y,x=xnozero), method='glm',se=F)
# pnap3<-ggdotscolor(x=toplot$xnozero, y=toplot$y, xlab='True individual fitness', ylab='Inferred fitness') %>% addggregression()# + stat_smooth(aes(y=y,x=xnozero), method='glm',se=F)
pnap<-iplotreal(Y$Fitness, Ypred[Y$row])
print(pnap)


################################################################################

Fmode=3

r1<-napMCMC(
          y=Y$Fitness,
          h=Y$row,
          A=Go,
          s=PropoS( length(SNPsub), 0.2),
          b=1,bmin=1,bmax=1,
          a=0.1,amin =0.1, amax=0.1,
          # mu=100,
          # mu=1,mumin=1,mumax=1,
          # mu=0.5,mumin=0.5,mumax=0.5,
          # mu=0.3,mumin=0.1,mumax=20,
          mu=100,mumin=10,mumax=150,
          epi=1, epimin=1,epimax=1,
          # epi=1.01, epimin=0.1,epimax=3,
          svar=0.2,svarmin=0.2,svarmax = 0.2,
          ss=0.01,ssmin=0.01,ssmax=0.01,
          p=0.85,

          # m=SNPsub-1,
          # n=(unique(Y$row)) -1,
          m=SNPsub,
          n=(unique(Y$row)),

          iterations = 1e4,
          Fitnessmode=Fmode,
          Priormode=3,
          Proposalmode=3
)

### Parameters
# parameterplot_(r1$parchain, r1$parnames)
parinf<-parametersummary(r1$parchain, r1$parnames)
print(parinf)
pinf<-parinf$statistics[,1]

### SNPs
sinf<-parametersummary(r1$chain , SNPsub)
s<-sinf$statistics[,1]
# smode<-MCMCglmm::posterior.mode(r1$chain)
hist(sinf$statistics[,1])

### Individuals
# Ypred<-Ey_go(Go[,SNPsub], sinf$statistics[,1], Fmode, epi= parinf$statistics['epi',1] ) * parinf$statistics["mu",1]
Ypred<-Ey_go(Go[,SNPsub], sinf$statistics[,1], Fmode)
Ypred <- Ypred* parinf$statistics["mu",1]
# hist(Ypred, breaks=100)
# hist(Y$rFitness, breaks=100)

pnap3<-iplotreal(Y$Fitness, Ypred[Y$row])
print(pnap3)

test_Likelihood(
          A=Go@address,
          y=Y$Fitness,
          h=Y$row,
          s=s,
          n=1:nrow(Go)-1,
          m=topcols-1,
          b=pinf['b'],
          a=pinf['a'],
          p=pinf['p'],
          mu=pinf['mu'],
          # mu=mean(Y$Fitness[Y$Fitness!=0]),
          Fitnessmode=Fmode,
          verbose=TRUE)

# toplot<-data.frame(x=Y$Fitness , y=Ypred[Y$row])
# toplot$xnozero <- toplot$x
# toplot$xnozero[toplot$xnozero==0] <-NA
# # pnap3<-ggplot(toplot, xlab='True individual fitness', ylab='Inferred fitness') + geom_point(aes(y=y,x=xnozero), color='darkgrey') %>% addggregression()# + stat_smooth(aes(y=y,x=xnozero), method='glm',se=F)
# pnap3<-ggdotscolor(x=toplot$xnozero, y=toplot$y, xlab='True individual fitness', ylab='Inferred fitness') %>% addggregression()# + stat_smooth(aes(y=y,x=xnozero), method='glm',se=F)



#####################################################################################


## GWA
starting<-BMgwa(Go[Y$row,SNPs], Y$Fitness, type=1)

gwapred<-(Go[Y$row,SNPs] %*% starting$coefficients) + starting$meanasintercept

toplot<-data.frame(x=Y$Fitness, y=gwapred)
toplot$xnozero <- toplot$x
toplot$xnozero[toplot$xnozero==0] <-NA
pgwa<-ggdotscolor(x=toplot$xnonzero, y=toplot$y, xlab='True individual fitness', ylab='Inferred fitness')  %>% addggregression()# + stat_smooth(aes(y=y,x=xnozero), method='glm',se=F)
print(pgwa)

iplotreal(Y$Fitness, gwapred)


# pgwa<-iplotreal(Y$rFitness, Y$row, BMpred(Go@address , starting$coefficients, 1:nrow(Go)-1, SNPs-1, starting$meanasintercept)  )

# lmod<-lm(My(fn(Y$rFitness),fn(Y$row)) ~ Go[,SNPs] )


# save_plot(plot_grid(pnap$pind,pgwa$pind),file="outplot.pdf", ncol=2,nrow=1, base_width = 5, base_height = 5)
save_plot(plot_grid(pnap$pind,pnap3$pind,pgwa$pind),file="outplot.pdf", ncol=3,nrow=1, base_width = 5, base_height = 5)


#
# cor.test(smode, starting$coefficients)
# plot(sinf$statistics[,1] , starting$coefficients)
# dev.off()
#
# plot()
# # plot(r0$chain[,1])

# ## Global parameters
# r0$parchain[,1] %>% plot
# r0$parchain[,2] %>% plot
# r0$parchain[,3] %>% plot
# r0$parchain[,4] %>% plot
# r0$parchain[,5] %>% plot


# # as.mcmc(r0$p) %>% plot(main='pi constant (frac of zeros)')
# # as.mcmc(r0$b) %>% plot(main='b constant = variance/mean')
# # as.mcmc(r0$a) %>% plot(main='a constant = variance intercapt')

# as.mcmc(r0$posterior) %>% plot(main='Posterior')
# as.mcmc(r0$chain) %>% plot(main='selection')

# ## Extract selection coefficients
# rcle<-removeburnin(r0,0.5)
# # rcle<-r0
# sinf<-MCMCglmm::posterior.mode(as.mcmc(rcle$chain))
# sinf<-apply(rcle$chain,2,mean)
# sinf_range<-HPDinterval(as.mcmc(rcle$chain))


# ## Plots

# qplot(x=topcols,y=sinf,ylim=c(0,1)) +
#   geom_segment(aes(x=topcols,xend=topcols,y=sinf_range[,1],yend=sinf_range[,2]) )+
#   # geom_abline(intercept = 0,slope = 1,lty="dotted")+
#   ylab("Inferred selection coefficients")+
#   xlab("SNP Position")

# toplot<-data.frame(y= Ey_go(Go[,topcols],s = sinf,mode=1)[Y$row],x=Y$oFitness)
# toplot$xnozero<-toplot$x
# toplot$xnozero[toplot$xnozero==0]<-NA
# ggplot(data=toplot,aes(y=y,x=x)) +
#   geom_point(col="grey") +
#   stat_smooth(aes(y=y , x=xnozero), se=F,lty="dashed",col="black")+
#         ylim(c(0,1.2)) + xlim(c(0,1.2))+
#         ylab("Inferred individual fitness")+
#         xlab("True individual fitness")+
#     geom_abline(intercept = 0,slope = 1,lty="dotted")


################################################################################
## MCMC with cross-validation
################################################################################
















################################################################################
## MCMC with mlp
################################################################################
#
#
# sourceCpp('MCMC.cpp')
# m=10
# topcols<-top_cols(genomes,p = m)
# length(topcols)
# n=nrow(Go)
# s=rep(0,m)
# s=runif(m)
#
#
# ################################################################################
# sourceCpp('MCMC.cpp')
# # Multplicative oriented
# it=1e4
# system.time(
#   r1<-as.mcmc(gwsMCMC(Y$oFitness, Y$row, Go@address,s,m=topcols, iterations = it,verbose=TRUE,TEST = FALSE, mode=1))
# )
# r1<-as.mcmc(r1)
# pdf(paste0('figs/MCMCmulti',it,'pdf')); plot(r1) ; dev.off();
#
# r1b<-as.mcmc(gwsMCMC(Y$oFitness, Y$row, Go@address,runif(m),m=topcols, iterations = it,verbose=FALSE,TEST = FALSE, mode=1))
#
# pdf('figs/MCMC_multi_1e4_converge.pdf'); gelman.plot(mcmc.list(r1,r1b)) ; dev.off();
#
# ################################################################################
# # Additive
# r2<-as.mcmc(gwsMCMC(Y$oFitness, Y$row, Go@address,s,m=topcols, iterations = it,verbose=FALSE,TEST = FALSE, mode=2))
# pdf(paste0('figs/MCMCadd',it,'pdf')); plot(r2) ; dev.off();
#
# ################################################################################
# # Inverse multiplicative
# r3<-as.mcmc(gwsMCMC(Y$oFitness, Y$row, Go@address,s,m=topcols, iterations = it,verbose=FALSE,TEST = FALSE, mode=5))
# pdf(paste0('figs/MCMCadd',it,'pdf')); plot(r3) ; dev.off()




# sourceCpp('MCMC.cpp')
# system.time(r2<-gwsMCMC(Y$oFitness, Y$row, Go@address,s,m=topcols, iterations = 1e5,verbose=FALSE, mode=1))
# r2<-as.mcmc(r2)
#
# # Additive oriented
# r2<-gwsMCMC(Y$oFitness, Y$row, X,s,
#             iterations = 1e5,verbose=FALSE, mode=2)  %>% as.mcmc
# # plot(r2)
#
# # Inverse Multiplicative oriented
# r3<-gwsMCMC(Y$oFitness, Y$row, X,s,
#             iterations = 1e5,verbose=FALSE, mode=2)  %>% as.mcmc





#lgamma((Ey(X[unique(Y$row),], s)^2+.Machine$double.xmin) / (Vy(Y$oFitness,Y$row)+.Machine$double.xmin))
#
#
#getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
#}
#sourceCpp('MCMC.cpp')
#plot(Vy(Y$oFitness,Y$row))
#mean(Vy(Y$oFitness,Y$row))
#getmode(Vy(Y$oFitness,Y$row))
#Vy(Y$oFitness,Y$row)
#Y$row %>% unique
#genomes$fam[398,]$sample.ID
#var(Y$oFitness)
#filter(Y,sample.ID==genomes$fam$sample.ID[unique(Y$row)[424]] ) %>% select(oFitness) %>% var
#filter()
#var(c(0.1))
#
#s=runif(10,0,1)
#Ey(X[unique(Y$row),], s)
##Esub(Ey(X[unique(Y$row),], s), Y$row-1, unique(Y$row))
#length(Ey(X[unique(Y$row),], s))
#
#
#Ey(X[unique(Y$row),], s)
#hsub(Y$row) %>% max
#unique(Y$row) %>% length
#
#
#length(Ey(X[unique(Y$row),], s))
#
#fit_lik(Y$oFitness, Y$row, X[unique(Y$row),],runif(10,0,1), Vy(Y$oFitness, Y$row))
#
## Multiplicative oriented
#r1<-gwsMCMC(Y$oFitness, Y$row, X[unique(Y$row),],runif(10,0,1),
#            iterations = 1e5,verbose=FALSE) %>% as.mcmc
#
#Ey(X[unique(Y$row),], runif(10,0,1))
#
#X
#
#dim(X[unique(Y$row),])
#plot(r1)


# # Additive oriented
# r2<-gwsMCMC(Y$oFitness, Y$row, X,runif(10,0,1),
#             iterations = 1e5,verbose=FALSE)  %>% as.mcmc
# plot(r2)
#
# # Inverse Multiplicative oriented
# r3<-gwsMCMC(Y$oFitness, Y$row, X,runif(10,0,1),
#             iterations = 1e5,verbose=FALSE)  %>% as.mcmc




# r1<-gwsMCMC(Y$rFitness, Y$row, X,rnorm(10,0,1),iterations = 1e6,verbose=FALSE) %>% as.mcmc()
# r2<-gwsMCMC(Y$rFitness, Y$row, X,rnorm(10,0,1),iterations = 1e6,verbose=FALSE) %>% as.mcmc()
# r3<-gwsMCMC(Y$rFitness, Y$row, X,rnorm(10,0,1),iterations = 1e6,verbose=FALSE) %>% as.mcmc()
#
# r1[,9] %>% hist(.,breaks=1000)
# r2[,9] %>% hist(.,breaks=1000)
# r3[,9] %>% hist(.,breaks=1000)
#
# r1[,1] %>% hist(.,breaks=4000)
# r1[,2] %>% hist
#
# rl<-mcmc.list(r1,r2,r3)
# gelman.diag(rl)
# gelman.plot(rl)
#
# r1[,7] %>% hist
# r1[,7] %>% plot
# r3[,7] %>% plot

# gammavalues(means^2/ Vy(Y$Fitness, Y$row))
# means=Ey( G[,1:10],rnorm(10,0,1))
# means[means<0]<-1e-38
# variances<-Vy(Y$Fitness, Y$row)
# variances[variances==0]<-1e+6
# gammavalues(means^2/ Vy(Y$Fitness, Y$row)+1e-38)

#
# r<-gwsMCMC(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1),iterations = 1e5,verbose=FALSE)
# r2<-gwsMCMC(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1),iterations = 1e5,verbose=FALSE)
# r3<-gwsMCMC(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1),iterations = 1e5,verbose=FALSE)
# r4<-gwsMCMC(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1),iterations = 1e6,verbose=FALSE)
#
# mt<-function(mcmc,last=10000){
#   apply(tail(mcmc,n=last),2,mean)
# }
# plot(mt(r4),mt(r))
# cor.test(mt(r4),mt(r))
#
# plot(mt(r4),mt(r2))
# plot(mt(r4),mt(r3))
#
# r10<-gwsMCMC(Y$Fitness, Y$row, G[,1:50],rnorm(50,0,1),
#              iterations = 1e5,verbose=FALSE)
# r11<-gwsMCMC(Y$Fitness, Y$row, G[,1:50],rnorm(50,0,1),
#              iterations = 1e5,verbose=FALSE)
#
# qplot(mt(r10,100),mt(r11,100))
#
# qplot(r10[10000,],r11[10000,])
# cor.test(r10[1e5,],r11[1e5,])
#
# cor.test(mt(r10),mt(r11))
# hist(r10[sample(1:1e6,10000),1])
# hist(r11[sample(1:1e6,10000),1])
#
# hist(r10[sample(1:1e6,10000),2])
# hist(r11[sample(1:1e6,10000),2])



#################################################################################
####### Simulated data
#################################################################################

# X<-G[,1:10]
# s<-rnorm(10,0,1.5)
# eh<-buildpheno(X,s)
# eh[eh<0]<- 1e-38
# eh%>% hist
# # yh1<-sapply(eh, function(r) {rexp(1,1/r)})
# # yh2<-sapply(eh, function(r) {rexp(1,1/r)})
# # yh3<-sapply(eh, function(r) {rexp(1,1/r)})
# yh1<-eh + rexp(length(eh),2)
# yh2<-eh + rexp(length(eh),2)
# yh3<-eh + rexp(length(eh),2)
# h<-rep(1:515,3)
# yh<-c(yh1,yh2,yh3)
# hist(yh)
#
# rx<-gwsMCMC(yh, h, G[,1:10],rnorm(10,0,1),iterations = 1e6,verbose=FALSE)
#
# library(coda)
# rxx<-coda::as.mcmc(rx)
# pdf('trialMCMC.pdf')
# plot(rxx)
# dev.off()
#
# plot(apply(rx,2,mean),s)
#
# s[10]

# #r<-gwsMCMC(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1),iterations = 1e6,verbose=FALSE)
# s_prior_unif(rnorm(10,0,1))
# fit_lik(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1),variances)
#
#
# r[1,] == r[10000,]
# plot(r[1,], r[10000,])
#
# plot(r[-c(1:1000),1])
# plot(r[-c(1:1000),2])
#
#
#
# plot(r[-c(1:10000),1])
# pred<-buildpheno(G[,1:10],r[1e6,])[Y$row]
# pred[pred<0]<-0
# plot(Y$Fitness,pred)
# hist(sapply(pred,function(r) rexp(1,r)))
#
#
# hist(r[-c(1:10000),1],breaks=5000)
# hist(r[-c(1:500),2])
#
# fit_lik(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1))
# fit_lik(Y$Fitness, Y$row, G[,1:10],rnorm(10,0,1))
#
# eh=-0.607184
# vh=2433.49
# yh=58
# Rcpp::cppFunction('
# double likfunction(double eh,double vh,double yh){
#         return ( pow(eh,2) - log(yh) - ((yh*eh)/vh) - log(tgamma(pow(eh,2)/vh)) - pow(eh,2)/vh * log(vh/eh));
# }')
#
# likfunctionR<-function(eh,vh,yh){
#   eh^2 - log(yh) - ((yh*eh)/vh) - lgamma(eh^2/vh) - eh^2/vh * log(vh/eh)
# }
# yh=1.17549e-38
# vh=4090
# eh=1e-38
# sourceCpp('MCMC.cpp')
# likfunction(eh,vh,yh)
# likfunctionR(eh,vh,yh)
# LnGamma(eh^2/vh + 1e-38 )
#
# vh/eh
#
# lgamma(eh^2/vh)
# LnGamma(eh^2/vh)
# log()
# gamma(3)
# lgamma(0.2)
#
# Vy(Y$Fitness,Y$row)
# means<-Ey(G[,1:10],rnorm(10,-.1,1))
# hist(means)
# means[means<0] <-0
# hist(means)
#
# unique(Y$row) %>% length()
# dim(G)

# imputeNAfit<-function(x){
#   x[is.na(x)]<-0
#   return(x)
# }
# subsetFitness=
#   fieldfilter(field.c,code='mlp') %>%
#   mutate(Fitness = imputeNAfit(Fitness)) %>%
#   mutate(Fitness = Fitness / mean(Fitness)) %>%
#   group_by(id) %>%
#   summarise(meanFitness=mean(Fitness,na.rm=TRUE),
#             varFitness=var(Fitness,na.rm=TRUE)) %>%
#   as.data.frame
#
# genomes<-curatedmerge(genomes, subsetFitness,c('Ymlp','Yvarmlp'))
# genomes$fam$Ymlp[is.na(fn(genomes$fam$Ymlp))]<-1
# genomes$fam$Yvarmlp[is.na(fn(genomes$fam$Yvarmlp))]<-0
#
#
# hist(dry$Fitness_mlp / max(dry$Fitness_mlp,na.rm=TRUE ))
# hist(dry$Fitness_mlp / (min(dry$Fitness_mlp,na.rm=TRUE )+.Machine$double.xmin) )
#
# hist(fn(genomes$fam$Ymlp))
# hist(fn(genomes$fam$Yvarmlp))
#
# X<-G[1:20,1:5]
# s<-rnorm(5,0,1)
#
# meancent.mat(X) * s
#
# (meancent.mat(X) * s) %>% apply(.,2,sum)
# ((X) * s) %>% apply(.,2,sum)
#
# # head(genomes$fam)
# # genomes$fam$Ymlp %>% is.na %>% table
# # genomes$fam$Yvarmlp %>% is.na %>% table
#
# library(Rcpp)
# sourceCpp('MCMC.cpp')
# findunique(c(1,1,2,3))
#
# r<-gwsMCMC(fn(genomes$fam$Ymlp),
#            fn(genomes$fam$Yvarmlp),
#            G[,1:10],
#            s_start(10,-5,+5),
#            iterations=1e6,
#            verbose=FALSE
# )
# r2<-gwsMCMC(fn(genomes$fam$Ymlp),
#            fn(genomes$fam$Yvarmlp),
#            G[,1:10],
#            s_start(10,-5,+5),
#            iterations=1e6,
#            verbose=FALSE
# )
# sourceCpp('MCMC.cpp')
# r3<-gwsMCMC(fn(genomes$fam$Ymlp),
#             fn(genomes$fam$Yvarmlp),
#             G[,1:50],
#             cgwa(G[,1:50],fn(genomes$fam$Ymlp)),
#             iterations=1e6,
#             verbose=FALSE
# )
#
#
# rs<-apply(r,2,mean)
# r2s<-apply(r2,2,mean)
# r3s<-apply(r3,2,mean)
#
# gwares<-gwa(G[,1:50],fn(genomes$fam$Ymlp))
# cgwares<-cgwa(G[,1:50],fn(genomes$fam$Ymlp))
#
# plot(gwares,cgwares)
# plot(r3s,cgwares)
# plot(r3s,gwares)
#
# hist(r3s)
#
# plot(rs,r2s)
# cor.test(rs,r2s)
#
# hist(r3[,1],breaks=50,col=transparent('black'))
# hist(r2[,1],breaks=50,col=transparent('blue'),add=TRUE)
#
# hist(r[,2],breaks=50,col=transparent('black'))
# hist(r2[,2],breaks=50,col=transparent('blue'),add=TRUE)
#
# plot(r[,2],type="lines")
#
# hist(r[,3],breaks=50,col=transparent('black'))
# hist(r2[,3],breaks=50,col=transparent('blue'),add=TRUE)
#
# hist(r[,4],breaks=50,col=transparent('black'))
# hist(r2[,4],breaks=50,col=transparent('blue'),add=TRUE)
#
#
# gwa(G[,1:10],fn(genomes$fam$Ymlp))
# plot(gwa(G[,1:10],fn(genomes$fam$Ymlp)),rs)
# cor.test(gwa(G[,1:10],fn(genomes$fam$Ymlp)),rs)
# cor.test(gwa(G[,1:10],fn(genomes$fam$Ymlp)),r2s)
#
# plot(
#   buildpheno(G[,1:50], r3s),
#   fn(genomes$fam$Ymlp)
# )
#
# buildpheno(G[,1:10], rs)
#
#
# plot(fn(genomes$fam$Ymlp),buildpheno(G[,1:10], r3s,type = 'multiplicative'))
# plot(fn(genomes$fam$Ymlp),buildpheno(G[,1:10], gwa(G[,1:10],fn(genomes$fam$Ymlp)),type = 'additive'))
# plot(fn(genomes$fam$Ymlp),buildpheno(G[,1:10], cgwa(G[,1:10],fn(genomes$fam$Ymlp)),type = 'additive'))
#
#
#
# r3<-gwsMCMC(fn(genomes$fam$Ymlp),
#             fn(genomes$fam$Yvarmlp),
#             G[,1:50],
#             cgwa(G[,1:50],fn(genomes$fam$Ymlp)),
#             iterations=1e6,
#             verbose=FALSE
# )


#### To delete
# likwrapper<-function(par,data){
#   y=mydata$y
#   yvar=mydata$yvar
#   X=mydata$X
#   s=par
#   return(fit_lik(y,yvar,X,s))
# }
# optim(par= gwa(G[,1:50],fn(genomes$fam$Ymlp)),
#       y=data,
#       fn = likwrapper,
#       method = 'BFGS'
#         )
# sourceCpp('MCMC.cpp')

# fit_lik(fn(genomes$fam$Ymlp),
#         fn(genomes$fam$Yvarmlp),
#         G[,1:50],
#         gwa(G[,1:50],fn(genomes$fam$Ymlp)))
#
#
# fit_lik(fn(genomes$fam$Ymlp),
#         fn(genomes$fam$Yvarmlp),
#         G[,1:10],
#         runif(10,-5,+5))
# Gamma_lik(0+minc , 4e-313+minc,0+minc)
# LnGamma(0+minc)
# rgamma(100,shape = 0 ,rate = 4.24909e-313)
# rgamma(100,shape = 0 ,rate = 4.24909e-313)
#
#
# RGamma_lik(1e-32,1e-32,1e-32)
# rgamma(100,scale = 1e-32,shape = 1e-32)
#
# exp(-1)
# dgamma(x = 1e-32,shape = 1e-32,rate = 1e-32,log = TRUE)
# Gamma_lik(1e-32,1e-32,1e-32)



#######
# Feedback from Peter
# Hamilton
# NUTS algorithm
#######
# sourceCpp('MCMC.cpp')
#
# s_prior_unif(c(.4,.3,.2))
# s_prior_unif(t(matrix(c(.4,.3,.2))))
#
# .t()
# Gamma_lik(alpha,beta,obs)
#
# sourceCpp('MCMC.cpp')
# r<-gwsMCMC(fn(genomes$fam$Ymlp),
#             fn(genomes$fam$Yvarmlp),
#             G[,1:10],
#             s_start(10,-5,+5),
#             iterations=1e5,
#             verbose=FALSE
#           )
# hist(r[,2])
#
# apply(r,2,hist)
# apply(r,2,summary)
# apply(r,2,mean)
# apply(r_2,2,mean)
#
# plot(apply(r,2,mean),apply(r_2,2,mean))
# # r_2<-r
# dev.off()
#
# par(mfrow=c(2,5))
# for(i in 1:ncol(r)){
#   hist(r[,i],breaks=100)
# }
# par(mfrow=c(2,5))
# for(i in 1:ncol(r)){
#   plot(r[,i],breaks=100,type='lines')
# }
# dev.off()
#
#
# apply(r,2,mean)
#
# View(r)
# head(r)
#
# hist(r[,10],breaks=1000)
#
# apply(r,2,mean)
# apply(r,2,function(i) t.test(i)$p.value)
#
#
# s_proposal_unif(c(.1,.2,.5))

# logacc=0
# for(i in 1:10){
#   logacc=sum(logacc + log(1/10))
# }
