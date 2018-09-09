################################################################################
## The aim is to simulate data under various models and re-estimate it
################################################################################

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


# sourceCpp('MCMC.cpp'); # old
sourceCpp('MCMC2.cpp');
sourceCpp('gwa.cpp');

####************************************************************************####
#### Genomematrix #####


genomes<-readRDS('../gws/databig/genome.rda')
# G <- attachgenomes(genomes)

Go<-bigmemory::attach.big.matrix('../gws/databig/gc.desc')



###************************************************************************####
### Single Run for profiling ####

set.seed(1)

## Preparate all parameters
p=0.1
b=0.1
a=0
ss=0.1
svar=0.1
m=10
epistasis=1
prophidden=0
training=0.9
training=1

sdat<-simulate_data(Go=Go,p=p,b=b,a=a,
										m=m,
										prophidden = prophidden,
										ss=ss,
										trainingp = training,
										svar = svar
										)

sdat$s

# res<-napMCMC(A=Go@address,
#              h=sdat$Y$row,
#              s = rep(0,m),
#              m=sdat$newtopcols-1,
#              n=sdat$training-1,
#              svar = 0.1,
#              ss=0.1,
#              y = sdat$Y$Y3,
#              Fitnessmode = 3,
#              # Priormode=3,
#              Priormode=3,
#              Proposalmode=3,
#              iterations = 1e3
#           )
res<-napMCMC(A=Go,
             h=sdat$Y$row,
             s = rep(0,m),
             m=sdat$newtopcols-1,
             n=sdat$training-1,
             svar = 0.1,
             ss=0.1,
             y = sdat$Y$Y3,
             Fitnessmode = 3,
             Priormode=3,
             Proposalmode=3,
             iterations = 1e3, TEST=T
          )

sinf<-apply(res$chain,2,median)
splot_(s = sdat$s,sinf=sinf)
parameterplot_(res$parchain,res$parnames)


#
# gres<-BMgwa2(X0= Go[sdat$training,sdat$newtopcols],
#        y = My(sdat$Y$Y3,sdat$Y$row),
#        type=3
#        )
#
# hist(sdat$s)
#
# splot_(s = sdat$s * sdat$ori,sinf=gres$coefficients)
#
# plot( sdat$Ey3, Go[, sdat$newtopcols] %*% gres$coefficients )
# plot( sdat$Ey3[sdat$testing], BMpred(Go@address,x = gres$coefficients, myrows = sdat$testing -1, mycols = sdat$newtopcols-1,intercept = 1) )

# gres1<-BMgwa2(X0= Go[sdat$training,sdat$newtopcols],
#        y = My(sdat$Y$Y3,sdat$Y$row),
#   type=1
#        )
#
# splot_(s = sdat$s,sinf=gres$coefficients)
# splot_(s = sdat$s,sinf=gres1$coefficients)
#
# cor.test(sdat$s,gres)
# lm(gres ~ sdat$s) %>% summary
#
# gwaind<-Go[sdat$training , sdat$newtopcols] %*% gres$coefficients + gres$meanasintercept
# plot(sdat$Ey3,gwaind);abline(0,1)
#
# # cor.test(sdat$s,sinf)
# # res$parchain[,6] %>% summary
# # res$parchain[,6] %>% qplot
#
#
# #### GRID Run ####

set.seed(1)

p=0.5
b=0.5
# a=0
a=.1
# ss=0.5
ss=0
svar=0.15
mu=1
m=100
m=1000
epistasis=1
prophidden=0
training=0.8
# training=1

runname<-paste0("grid",
                '_m=',m,
               'p=',p,
               '_b=',b,
               '_a=',a,
               '_ss=',ss,
               "_svar=",svar,
               "_hidden=",prophidden,
               '_mu=',mu,
               '_e=',epistasis,
               "_train=",training
								)
runname

sdat<-simulate_data(Go=Go,
                    p=p,
                    b=b,
                    a=a,
                    mu=mu,
                    epi = epistasis,
                    m=m,
                    prophidden = prophidden,
                    svar = svar,
                    ss=ss,
                    trainingp = training
                    )

sdat<-run_MCMC_GWA(sdat,
                    Go,
                    simumodes=c(1,2,3),
                    fitnessmodes=c(1,2,3),
                    priormodes=c(3),
                    proposalmodes=c(3),
                    gwatypes=c(1,3),
                    iterations=1e4,
                    burnin=0.1,
                    donap=TRUE,
                    dogwa=TRUE,
                    doplot=TRUE,
                    pdfname='out')
str(sdat)



# # qplot(sdat$s)
# # qplot(sdat$Ey1)
# # # qplot(sdat$Y$Y1)
# # qplot(sdat$Ey2)
# # # qplot(sdat$Y$Y2)
# # qplot(sdat$Ey3)
# # qplot(sdat$Y$Y3)#
# sdat<-run_MCMC_GWA(sdat,
#                     Go,
# 										simumodes=c(1,2,3),
# 										fitnessmodes=c(1,2,3),
# 										priormodes=c(4),
# 										proposalmodes=c(3),
# 										gwatypes=c(1,3),
#                     # iterations=5e3,
# 										iterations=1e4,
# 										# iterations=1e5,
#                     burnin=0.1,
#                     donap=TRUE,
# 									  dogwa=TRUE,
# 									  doplot=TRUE,
# 									  pdfname=runname)
#
# sdat<-run_MCMC_GWA(sdat,
#                     Go,
# 										simumodes=c(3),
# 										fitnessmodes=c(3),
# 										priormodes=c(4),
# 										proposalmodes=c(3),
# 										gwatypes=c(1,3),
#                     # iterations=5e3,
# 										iterations=1e4,
# 										# iterations=1e5,
#                     burnin=0.1,
#                     donap=TRUE,
# 									  dogwa=FALSE,
# 									  doplot=TRUE,
# 									  pdfname=runname)
#
# extract(sdat,3,3,4,3)$plot_s


# ####************************************************************************####
# #### All parameters in grid ####
#
# ### SIMULATION PARAMETERS
# Ps= c(0,.5,.85)
# Bs= c(0.1,0.5)
# As= c(0,0.15)
ps= c(0.5,0.0)
bs= c(0.25)
as= c(0.1)
sss=c(0,0.5,0.8)
svars=c(0.1,0.3)
mus=c(1)
prophid<-c(0.2,0.0)
ms=c(100,10,1000)
# ms=c(100,10)
episs=1


params<-expand.grid(ps,bs,as,sss,svars,mus,ms, prophid,episs);
dim(params)
colnames(params) <- c('p','b','a','ss',"svar",'mu',"m","hidden",'epistasis');


iterations=100
# iterations=5e4
iterations=1e4
i=1

for(i in 1:nrow(params)){
# ## Grid Parameters
p=params[i,'p']
b=params[i,'b']
a=params[i,'a']
ss=params[i,'ss']
svar=params[i,'svar']
prophidden=params[i,'hidden']
m=params[i,'m']
mu=params[i,'mu']
mu=1
epistasis=1
trainingp=0.8
# training=1


# Define name
runname<-paste0("GRID",
                '_m=',m,
               '_p=',p,
               '_b=',b,
               '_a=',a,
               '_ss=',ss,
               "_svar=",svar,
               "_hidden=",prophidden,
               '_mu=',mu,
               '_e=',epistasis,
               "_train=",trainingp
								)
runname

# simulate
sdat<-simulate_data(Go=Go,
										p=p,
										b=b,
										a=a,
                           mu=mu,
										epistasis = epistasis,
										m=m,
										prophidden = prophidden,
										svar = svar,
										ss=ss,
										trainingp = trainingp
										)

sdat<-run_MCMC_GWA(sdat,
                        Go,
                        simumodes=c(1,2,3),
                        fitnessmodes=c(1,2,3),
                        # simumodes=c(3),
                        # fitnessmodes=c(1,3),
                        priormodes=c(3),
                        proposalmodes=c(3),
                        gwatypes=c(1,3),
                        iterations=iterations,
                        # iterations=1e4,
                        burnin=0.1,
                        donap=TRUE,
                        dogwa=TRUE,
                        doplot=TRUE,
                        pdfname=paste0("figs/mcmc/",runname))


}

# grid1<-plotgrid(sdat,
# simumodes=c(3),
# fitnessmodes=c(1,3),
# priormodes=c(3),
# proposalmodes=c(3),
# gwatypes=c(3),
# plottype='s')
#
# grid2<-plotgrid(sdat,
# simumodes=c(3),
# fitnessmodes=c(1,3),
# priormodes=c(3),
# proposalmodes=c(3),
# gwatypes=c(3),
# plottype='i')
#
#
# allgrids<-plotgrid(sdat,
# simumodes=c(1,2,3),
# fitnessmodes=c(1,2,3),
# # simumodes=c(3),
# # fitnessmodes=c(1,3),
# priormodes=c(3),
# proposalmodes=c(3),
# gwatypes=c(1,3),
# plottype='all',
# pdf='out.pdf')


