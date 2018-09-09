## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(coda)
library(Rcpp)

library(moiR)
# install('~/moiR')
load_all('.')


####************************************************************************####
#### Genomematrix #####


genomes<-readRDS('databig/genome.rda')
# G <- attachgenomes(genomes)

Go<-bigmemory::attach.big.matrix('databig/gc.desc')


####************************************************************************####
#### Fitness mlp ####

## In file data-cleaning/Yfitness.R
# Y<-readRDS('data/Y.rds')


# sourceCpp('MCMC.cpp')
sourceCpp('MCMCmean.cpp');
sourceCpp('gwa.cpp');


####************************************************************************####
#### Trial random forest -- does not work #####
# library(randomForest)
# rf<-randomForest(Go[sdat$training,sdat$newtopcols],y = My(sdat$Y$Y3,sdat$Y$row))
# plot( predict(rf,newdata = Go[sdat$testing,sdat$newtopcols] )  ~ sdat$Ey3[sdat$testing])
# cor.test( predict(rf,newdata = Go[sdat$testing,sdat$newtopcols] ) , sdat$Ey3[sdat$testing] )

####************************************************************************####
#### Single Run for profiling ####
#
# set.seed(1)
#
# ## Preparate all parameters
PI=0.1
BE=0.1
AI=0
ss=0.1
svar=0.1
svar=0.5
# svar=0.3
m=10
epistasis=1
prophidden=0
training=0.9
training=1

# # simulate
sdat<-simulate_data(Go=Go,PI=PI,BE=BE,AI=AI,
										m=m,
										prophidden = prophidden,
										ss=ss,
										trainingp = training,
										svar = svar
										)

## Example one run

res<-napMCMC(A=Go@address,
             h=sdat$Y$row,
             s = rep(0,m),
             m=sdat$newtopcols-1,
             n=sdat$training-1,
             svar = 0.1,
             ss=0.1,
             y = sdat$Y$Y3,
             Fitnessmode = 3,
             Priormode=3,
             # Priormode=4,
             Proposalmode=3,
             iterations = 5e3
          )
sinf<-apply(res$chain,2,median)
splot_(s = sdat$s,sinf=sinf)
parameterplot_(res$parchain,res$parnames)

gres<-BMgwa2(X0= Go[sdat$training,sdat$newtopcols],
       y = My(sdat$Y$Y3,sdat$Y$row),
       type=3
       # type=1
       )
gres1<-BMgwa2(X0= Go[sdat$training,sdat$newtopcols],
       y = My(sdat$Y$Y3,sdat$Y$row),
  type=1
       )

splot_(s = sdat$s,sinf=gres$coefficients)
splot_(s = sdat$s,sinf=gres1$coefficients)

cor.test(sdat$s,gres)
lm(gres ~ sdat$s) %>% summary

gwaind<-Go[sdat$training , sdat$newtopcols] %*% gres$coefficients + gres$meanasintercept
plot(sdat$Ey3,gwaind);abline(0,1)

# cor.test(sdat$s,sinf)
# res$parchain[,6] %>% summary
# res$parchain[,6] %>% qplot


#### GRID Run ####

set.seed(1)

## Preparate all parameters
PI=0.5
BE=0.1
AI=0
ss=0.1
svar=0.1
mu=1
m=10
# m=100
epistasis=1
prophidden=0
training=0.8
# training=1

# Define name
runname<-paste0("grid",
                '_m=',m,
               'p=',PI,
               '_b=',BE,
               '_a=',AI,
               '_ss=',ss,
               "_svar=",svar,
               "_hidden=",prophidden,
               '_mu=',mu,
               '_e=',epistasis,
               "_train=",training
								)
runname

# simulate
sdat<-simulate_data(Go=Go,PI=PI,BE=BE,AI=AI,
                    MU=mu,
										epistasis = epistasis,
										m=m,
										prophidden = prophidden,
										svar = svar,
										ss=ss,
										trainingp = training
										)

# qplot(sdat$s)
# qplot(sdat$Ey1)
# # qplot(sdat$Y$Y1)
# qplot(sdat$Ey2)
# # qplot(sdat$Y$Y2)
# qplot(sdat$Ey3)
# qplot(sdat$Y$Y3)

sdat<-run_MCMC_GWA(sdat,
                    Go,
										simumodes=c(1,2,3),
										fitnessmodes=c(1,2,3),
										priormodes=c(3),
										proposalmodes=c(3),
										gwatypes=c(1,3),
                    # iterations=5e3,
										iterations=1e4,
										# iterations=1e5,
                    burnin=0.1,
                    donap=TRUE,
									  dogwa=TRUE,
									  doplot=TRUE,
									  pdfname=runname)
# plotgrid(sdat,
#           simumodes=c(3),
# 					fitnessmodes=c(3),
# 					priormodes=c(3),
# 					proposalmodes=c(3),
# 					gwatypes=c(1,3)
# 					)




# # saveRDS(file=paste0(runname,".rda"),sdat)
# Simumode=1
# Simumode=2
# Fitnessmode=1
# Priormode=3
# Proposalmode=3
# gwatype=3
#
#
# extract(sdat,1,Fitnessmode,3,Proposalmode)$plot_s
# extract(sdat,2,Fitnessmode,3,Proposalmode)$plot_s
# extract(sdat,Simumode,gwatype = 1)$plot_s
# extract(sdat,Simumode,gwatype = 1)$plot_s
#
# extract(sdat,1,1,1,1)$plot_i
# extract(sdat,Simumode,Fitnessmode,3,Proposalmode)$plot_i
# extract(sdat,Simumode,gwatype = 1)$plot_i
# extract(sdat,Simumode,gwatype = 3)$plot_i
#
# p1<-qplot(main="1")
# 		p2<-qplot(main="2")
# 		p3<-qplot(main="3")
# 		p4<-qplot(main="4")
# plot_grid(p1,p2,p3,p4)
# plot_grid(p1,p2,p3,p4,byrow=T)

# sdat$results[[Simumode]][[Fitnessmode]][[Priormode]][[Proposalmode]] [["plot_s"]]
# sdat<-splot(sdat,Simumode = 1,gwatype = 1)
# sdat$results[[Simumode]][["GWA1"]] [["plot_s"]]
# sdat$results[[Simumode]][["GWA1"]] [["R2_B_i"]]
#
# sdat<-iplot(sdat,Simumode = 1,Fitnessmode = 1,Priormode = 1,Proposalmode = 1)
# sdat$results[[Simumode]][[Fitnessmode]][[Priormode]][[Proposalmode]] [["plot_i"]]
# sdat<-iplot(sdat,Simumode = 1,gwatype = 1)
# sdat$results[[Simumode]][["GWA1"]] [["plot_i"]]
# sdat$results[[Simumode]][["GWA1"]] [["R2_B_i"]]


#
# plot(apply(sdat$chainF1_S1$chain,2,median),
#      x=sdat$s)
# plot(sdat$rGWA_S1, x=sdat$s)
#
# plot((Go[,sdat$newtopcols] %*% sdat$rGWA_S1)+1 ,
#   x=sdat$Ey1 )
# plot(Ey_go(Go[,sdat$newtopcols],apply(sdat$chainF1_S1$chain,2,median),mode=1),
#   x=sdat$Ey1 )

# mymcmc<-sim_model_MCMC(
#                       Go = Go, ## for fitness inverse mult
#                       PI=PI,BE=BE,AI=AI,
#                       m=m,ss=ss,
#                       prophidden =  0,
#                       fitnessmode = fitnessmode,
#                       fitnessmodesim = fitnessmodesim,
#                       proposalmode=proposalmode,
#                       # priormode=priormode,
#                       priormode=1,
#                       epistasis,
#                       iterations=iterations,
#                       file2sink=paste0('figs/mcmc/MCMC_',runname,".log"),
#                       wannasink=FALSE
#                       )

# # Plot results
# splot<-plot.s.MCMC(mcmcobj = mymcmc)
# splot$res
# splot$plot
# gplot<-plot.global.MCMC(mymcmc,PI,BE,AI)
# gplot
#
#
# ####************************************************************************####
# #### All parameters in grid ####
#
# ### SIMULATION PARAMETERS
# Ps= c(0,.5,.85)
# Bs= c(0.1,0.5)
# As= c(0,0.15)
Ps= c(0,0.5)
Bs= c(0.25)
As= c(0.0)
sss=c(0.1,0.8)
svars=c(0.1,0.3)
mus=c(1)
# ms=c(10,100)
ms=c(10)
episs=1


params<-expand.grid(Ps,Bs,As,sss,svars,mus,ms, prophid,modfits,propos,episs);
dim(params)
colnames(params) <- c('p','b','a','ss',"svar",'mu',"m","hidden",'epistasis');

# results<-data.frame(params)
# results$r_s<-NA
# results$b_s<-NA
# results$r_ind<-NA
# results$b_ind<-NA
# rescols<-c('r_s','b_s','r_ind','b_ind')



iterations=1e4

for(i in 1:nrow(params)){
# ## Grid Parameters
PI=params[i,'PI']
BE=params[i,'BE']
AI=params[i,'AI']
ss=params[i,'ss']
ss=params[i,'svar']
prophidden=params[i,'hidden']
m=params[i,'m']
mu=1
epistasis=1
training=0.8
# training=1

# Define name
runname<-paste0("grid",
                '_m=',m,
               'p=',PI,
               '_b=',BE,
               '_a=',AI,
               '_ss=',ss,
               "_svar=",svar,
               "_hidden=",prophidden,
               '_mu=',mu,
               '_e=',epistasis,
               "_train=",training
								)
runname

# simulate
sdat<-simulate_data(Go=Go,PI=PI,BE=BE,AI=AI,
                    MU=mu,
										epistasis = epistasis,
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
										priormodes=c(4),
										proposalmodes=c(3),
										gwatypes=c(1,3),
                    # iterations=5e3,
										iterations=1e4,
										# iterations=1e5,
                    burnin=0.1,
                    donap=TRUE,
									  dogwa=TRUE,
									  doplot=TRUE,
									  pdfname=runname)

}

# for(i in 1:nrow(params)){
#
# ## Grid Parameters
# PI=params[i,'PI']
# BE=params[i,'BE']
# AI=params[i,'AI']
# ss=params[i,'ss']
# m=params[i,'m']
# fitnessmode=params[i,'fitnessmode']
# proposalmode=params[i,'proposalmodel']
# epistasis=params[i,'epistasis']
#
#
# # Define name
# runname<-paste0('p=',PI,
#                '_b=',BE,
#                '_a=',AI,
#                '_ss=',ss,
#                '_m=',m,
#                "_fitsim=",fitnessmodesim,
#                "_fit=",fitnessmode,
#                '_prop=',proposalmode,
#                '_e=',epistasis,
#                '_iter=',iterations)
# runname
# message("running ",runname)
#
# # Run MCMC
#
# mcmc<-sim_model_MCMC(
#                       # Go = Go,
#                       Go = Go, ## for fitness inverse mult
#                       # Gdouble = Gdouble,
#                       PI=PI,BE=BE,AI=AI,
#                       m=m,ss=ss,
#                       prophidden =  0.0,
#                       fitnessmode = fitnessmode,
#                       fitnessmodesim = fitnessmodesim,
#                       proposalmode=proposalmode,
#                       epistasis,
#                       iterations=iterations,
#                       file2sink=paste0('figs/mcmc/MCMC_',runname,".log"))
#
# # Plot results
# splot<-plot.s.MCMC(mcmc)
# splot$res
# splot$plot
# gplot<-plot.global.MCMC(mcmc,PI,BE,AI)
# gplot
#
# dev.off()
#
# # mcmc$chain$parchain %>% head
# # mcmc$chain$parchain %>% tail
#
# pdf(paste0('figs/mcmc/MCMC_',runname,".pdf"), width = 20,height = 10)
#   print(splot$plot)
#   print(gplot)
# dev.off()
# results[i,rescols] <-splot$res
#
# }
#
# write.tsv(results,file = "tables/MCMC_simulateddata.tsv")


