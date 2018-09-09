#######################################################################################################################################
####  MOTIVATION ####
# The aim is to model environmentally the presence and absence of alleles, and use the fitted models to predict the future change of
# alleles.
# There are several implementations and options, but the most important is that several assumptions of migration can be made:
# no migration - it fits lat long in the model
# intermediate - controled by the values in 3 genetic PCA axis that have been previously modeled environmentally
# universal migration - no controls, thus alleles can be predicted anywhere as long as the environmental change is large enough.

#######################################################################################################################################

##PARSE COMMAND LINE
library('argparser')

cat("\n")
p <- arg_parser("Simulation or fitness and inference with a Non-Additive Polygenic MCMC:")
p <- add_argument(p, "--p", help="Proportion of zeros", default=0.5)
p <- add_argument(p, "--b", help="Mean-sd relationship", default=0.2)
p <- add_argument(p, "--a", help="Intercept sd", default=0.1)
p <- add_argument(p, "--ss", help="Sparsity of zeros in selection coefficients", default=0)
p <- add_argument(p, "--svar", help="Sd of the log(1+s) Normal distribution", default=0.1)
p <- add_argument(p, "--mu", help="Mean fitness of the population", default=1)
p <- add_argument(p, "--m", help="Number of SNPs", default=10)
p <- add_argument(p, "--epi", help="Epistatic parameter", default=1)
p <- add_argument(p, "--hidden", help="Proportion of hidden genotypes that are also causal", default=0)
p <- add_argument(p, "--training", help="Proportion of the data to train", default=0.8)
p <- add_argument(p, "--iterations", help="Iterations of the MCMC", default=1e4)
p <- add_argument(p, "--randomgenome", help="Want to generate a random genome without LD", default=FALSE)
argv<-parse_args(p)

# argv <-list(gen="gene",nameanalysis="run2_gene_pca",PCAcontrol=T,limitlayers=T,geocontrol=F,runtype="classification")
# cat(paste("\n ----------------------------------------------------- \n"))
# cat("\n Interpreting arguments:\n")
# genetics=argv$gen
#   cat(paste("\n  --> The SNPs subset modeled will be of type: ",genetics,"\n"))
# nameanalysis=argv$nameanalysis
#   cat(paste("\n  --> The run name is: ",nameanalysis,"\n"))
# makeparallel=argv$makeparallel
#   cat(paste("\n    Computation in parallel?",makeparallel,"\n"))
# PCAcontrol=argv$PCAcontrol
#   cat(paste("\n    Control predictions by PCA structure?",PCAcontrol,"\n"))
# geocontrol=argv$geocontrol
#   cat(paste("\n    Control predictions by geographic limits?",geocontrol,"\n"))
# runtype=argv$runtype
#   cat(paste("\n    The modeling will be of type: ",runtype,"\n"))
# limitlayers=argv$limitlayers
# cat(paste("\n    The prediction will be limited in only important layers: ",limitlayers,"\n"))
# cat(paste("\n ----------------------------------------------------- \n"))


p=argv$p
b=argv$b
a=argv$a
ss=argv$ss
svar=argv$svar
mu=argv$mu
m=argv$m
epi=argv$epi
prophidden=argv$hidden
training=argv$training
iterations=argv$iterations

#####**********************************************************************#####

## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(coda)
library(Rcpp)

library(moiR)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# install('.')
# library('gws')
load_all('.')

# sourceCpp('MCMC2.cpp');
sourceCpp('MCMC.cpp')
sourceCpp('gwa.cpp');

####************************************************************************####
#### Genomematrix #####

genomes<-readRDS('../gws/databig/genome.rda')
# G <- attachgenomes(genomes)
Go<-bigmemory::attach.big.matrix('../gws/databig/gc.desc')


####************************************************************************####
#### GRID Run ####

runname<-paste0("grid",
                '_m=',m,
               'p=',p,
               '_b=',b,
               '_a=',a,
               '_ss=',ss,
               "_svar=",svar,
               "_hidden=",prophidden,
               '_mu=',mu,
               '_e=',epi,
               "_train=",training
                )
runname
message(runname)


sdat<-simulate_data(Go=Go,
                    p=p,
                    b=b,
                    a=a,
                    mu=mu,
                    epi = epi,
                    m=m,
                    prophidden = prophidden,
                    svar = svar,
                    ss=ss,
                    trainingp = training
                    )

# sdat$Ey1 %>% hist
# sdat$Y[,1] %>% hist
#
# sdat<-run_MCMC_GWA(sdat,
#                         Go,
#                         simumodes=c(2),
#                         fitnessmodes=c(2),
#                         priormodes=c(3),
#                         proposalmodes=c(3),
#                         gwatypes=c(3),
#                         iterations=iterations,
#                         burnin=0.1,
#                         donap=TRUE,
#                         dogwa=F,
#                         doplot=T,
#                         pdfname=paste0("out"))
#
# plot_grid(
# plotgrid(sdat,
#     simumodes=c(2),
#     fitnessmodes=c(2),
#     priormodes=c(3),
#     proposalmodes=c(3),
#     gwatypes=c(3),
# plottype = 's'
# ) ,
# plotgrid(sdat,
#     simumodes=c(2),
#     fitnessmodes=c(2),
#     priormodes=c(3),
#     proposalmodes=c(3),
#     gwatypes=c(3),
# plottype = 'i'
# )
# )


sdat<-run_MCMC_GWA(sdat,
                        Go,
                        simumodes=c(1,2,3),
                        fitnessmodes=c(1,2,3),
                        priormodes=c(3),
                        proposalmodes=c(3),
                        gwatypes=c(1,3),
                        iterations=iterations,
                        burnin=0.1,
                        donap=TRUE,
                        dogwa=TRUE,
                        doplot=TRUE,
                        pdfname=paste0("figs/mcmc/",runname))
plotgrid(sdat,
         simumodes=c(1,2,3),
        fitnessmodes=c(1,2,3),
        priormodes=c(3),
        proposalmodes=c(3))
plotgrid(sdat,
         simumodes=c(1,2,3),
        fitnessmodes=c(1,2,3),
        priormodes=c(3),
        proposalmodes=c(3),
        plottype='i')

