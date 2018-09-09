## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(moiR)
load_all('.')


################################################################################
# Fitness mlp

Y=dry[,c("Fitness_mlp")]
Y=merge(data.frame(genomes$fam),dry[,c('id',"Fitness_mlp")],by.x='sample.ID',by.y='id',all.x=TRUE)
head(Y)

###############################################################################
## Genomematrix
genomes<-readRDS('databig/genome.rda')
G <- attachgenomes(genomes)
G[1:10,1:10]

# Gi<-deepcopy(G)
# library(Rcpp)
# sourceCpp('matrixaccess.cpp')
# imNA(Gi)
# Gi[1:5,1:5]



################################################################################
# Fitness mlp top hits

gwares<-readRDS('dataint/rFitness_mlp.rda')
maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))

subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[501]),]
head(subg)
subg$rs<-paste(subg$chr,subg$pos,sep='_')

bsr<-gws::.read_assoc('output/rFitness_mlp.param.txt')
bsr<-gws::.read_assoc('output/Fitness_mlp.param.txt')
head(bsr)
hist(bsr$effectsize)
gwar<-readRDS('dataint/rFitness_mlp.rda')
gwar<-readRDS('dataint/Fitness_mlp.rda')
head(gwar)

hist(gwar$beta)


################################################################################
# Predictions of fitness

betagwa<-merge(genome$map,gwar[,c('chr','pos','beta')],by.x=c('chromosome','physical.pos'),by.y=c('chr','pos'),all.x=TRUE)
betagwa$beta[is.na(betagwa$beta)]<-0 # some SNPs are not estiamted in GEMMA
head(betagwa)
hist(betagwa$beta)
dim(betagwa)
dim(G)
betabs<-merge(genome$map,bsr[,c('chr','pos','beta')],by.x=c('chromosome','physical.pos'),by.y=c('chr','pos'),all.x=TRUE)
betabs$beta[is.na(betabs$beta)]<-0 # some SNPs are not estiamted in GEMMA
hist(betabs$beta)

Ygwa<-spred(G@address,betagwa$beta)
Ybs<-spred(G@address,betabs$beta)

summary(betagwa$beta)
summary(betabs$beta)

hist( G[,1:100] %*% betagwa$beta[1:100] )
hist( G[,1:100] %*% betabs$beta[1:100] )

hist(Ygwa)
hist(Ybs)
hist(relative(Y$Fitness_mlp))

qplot(Ygwa,x=relative(Y$Fitness_mlp)) + stat_smooth(color='darkred')
qplot(Ybs,x=relative(Y$Fitness_mlp)) + stat_smooth(color='darkred')


################################################################################
# Conditional GWA
Y$Fitness_mlp

res<-cGWA(xpA = G@address,B = Y$Fitness_mlp)


# Map <- genomes$map
# head(Map)
# Map$rs<-paste(Map$chr,Map$physical.pos,sep='_')
# Map$order<-1:nrow(Map)
#
# # Merge with top hits
# matched<-merge(Map, subg, by='rs')
# dim(matched)
# head(matched)



