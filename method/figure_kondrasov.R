
devtools::load_all('.')

################################################################################
##### Hamming distances between acessions
################################################################################

MEANNOZERO<-function(x){
  xnona<-na.omit(x)
  xnona<-xnona[xnona!=0]
  if(!is.null(length(xnona))){
    xmean=mean(xnona)
  }else{
    xmean=NA
  }
  return(xmean)
}
Y<-readRDS('data/Y.rds')
Ymean<-tapply(Y$Fitness, Y$sample.ID, MEANNOZERO)
hist(Ymean)
table(is.na(Ymean))

# Ydist<-dist(Y$Fitness)
# dim(as.matrix(Ydist))
# Ydist[1:5,1:5]
data('dry')
data('fam')
head(fam)
Y<-merge(fam, dry, by='id',all.x=T)
dim(Y)
Ymean<-Y[,'Fitness_mli']

Ydist<-as.matrix(dist(Ymean))
Ydist[1:5,1:5]
Ylower<-Ydist[lower.tri(Ydist)]

#####**********************************************************************#####
####  genetics
# D<-make_hamming_distances(path='../plink',plinkfile = '515g')
D<-make_hamming_distances(path='plinks',plinkfile = '515g--mli-bslmm-BAY')
D[1:5,1:5]
Dlower<-D[lower.tri(D)]

length(Dlower)
length(Ylower)

#####**********************************************************************#####

####  the type of comparison
kgroupdist<-matrix(ncol=ncol(D), nrow=nrow(D))
kgroupdist[]<-0
kgroupdist[which(Y$kgroup ==5),]<-1
kgroupdist[,which(Y$kgroup ==5)]<-1
  allrel<-expand.grid(which(Y$kgroup ==5),which(Y$kgroup ==5))
kgroupdist[allrel[,1],allrel[,2]]<-2
table(kgroupdist)
klower<-kgroupdist[lower.tri(kgroupdist)]


#####**********************************************************************#####
####  toplot data
toplot=data.frame(Dlower,Ylower,type=klower)

# toplot$Ycut<- cut(toplot$Ylower ,10)
# toplot$Dcut <- cut(toplot$Dlower ,10)
#
#
# toplotcut<-lapply(1:10, function(i){
#   tmp<-dplyr::filter(toplot,Ycut)
# })

# toplot=dplyr::filter(data.frame(Dlower,Ylower), Ylower!=0)
# toplot$type<-toplot$Dlower > 4.7e5

ggplot(data=toplot) +
  geom_point(aes(y=Ylower, x=Dlower, group=type))+
  xlab('# mutations distance')+
  ylab('# offspring distance')+
  # xlim(c(0,4.5e5))+
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F) +
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F,  formula = 'y ~ 0+  I(x^2)', method='glm') +
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F,  formula = 'y ~ 0+  I(x^0.1)', method='glm') +
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F,  formula = 'y ~ 0+ I(x)', method='glm', lty='dashed')

glm(data=filter(toplot,type==0),formula='Ylower ~ 0+  I(Dlower^2)')$aic
glm(data=filter(toplot,type==0),formula='Ylower ~ 0+  I(Dlower^1.5)')$aic
glm(data=filter(toplot,type==0),formula='Ylower ~ 0+  I(Dlower^0.1)')$aic
glm(data=filter(toplot,type==0),formula='Ylower ~ 0+  I(Dlower)')$aic


ggplot(data=filter(toplot,type==0,Dlower<4.6e5)) +
  geom_point(aes(y=Ylower, x=Dlower, group=type))+
  xlab('# mutations distance')+
  ylab('# offspring distance')+
  # xlim(c(0,4.5e5))+
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F,  formula = 'y ~ 0+  I(x^20)', method='glm') +
  # stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
  #              se=F,  formula = 'y ~ 0+  I(1/-x)', method='glm') +
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F,  formula = 'y ~ 0+  I(1/exp(x))', method='glm') +
  stat_smooth( aes(y=Ylower, x=Dlower,group=type,color=as.factor(type)),
               se=F,  formula = 'y ~ 0+ I(x)', method='glm', lty='dashed')


glm(data=filter(toplot,type==1),formula='Ylower ~ 0+  I(Dlower^2)') %>% summary
glm(data=filter(toplot,type==1),formula='Ylower ~ 0+ I(Dlower)') %>% summary
glm(data=filter(toplot,type==2),formula='Ylower ~ 0+  I(Dlower^2)') %>% summary
glm(data=filter(toplot,type==2),formula='Ylower ~ 0+ I(Dlower)') %>% summary
glm(data=filter(toplot,type==0),formula='Ylower ~ 0+  I(Dlower^2)') %>% summary
glm(data=filter(toplot,type==0),formula='Ylower ~ 0+ I(Dlower)') %>% summary

glm(data=filter(toplot,type==0,Dlower<4.6e5),formula='Ylower ~ 0+  I(Dlower^3)') %>% summary
glm(data=filter(toplot,type==0,Dlower<4.6e5),formula='Ylower ~ 0+ I(Dlower)') %>% summary




#####**********************************************************************#####
####  below wrong
###############################################################################
### With some real data

## Best positions
# gwares<-readRDS('dataint/Fitness_mlp.rda')
# gwares$order<-1:nrow(gwares)
# maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
# maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))
#
# subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[1000]),]
# head(subg)
# subg$rs<-paste(subg$chr,subg$pos,sep='_')
#

################################################################################
##### Raw counts of all SNPs
################################################################################
#
# library(gws)
# library(bigmemory)
# library(Rcpp)
#
# # Get the big matrix
# # genome<-readRDS('databig/genome.rda')
# Go<-bigmemory::attach.big.matrix('databig/go.desc')
#
# # Do stuff implemented in matrixaccess
#
# sourceCpp("matrixaccess.cpp")
# # mutations<-BigRowSums(Go@address)
# # mutations<-data.frame(row=1:nrow(Go),mutations=mutations$mutations)
# # saveRDS(file='data/mutations.rds',mutations)
# mutations<-readRDS('data/mutations.rds')
# head(mutations)
#
# ################################################################################
# ## Get phenotypes
# Y<-readRDS('data/Y.rds')
#
# ## As per mean?
# # sourceCpp("MCMC.cpp")
# # Ymean<-My(Y$Fitness, Y$row)
# # Ymean[Ymean==-9]<-0
#
# ################################################################################
#
# toplot<-merge(Y,by.x='row',mutations,by.y='row',all.y=TRUE)
# toplot$Fitnessnonzero<-ifelse(toplot$Fitness == 0, NA, toplot$Fitness)
#
# head(toplot)
#
# toplot$mutbins<-toplot$mutations
# toplot$mutbins<-round(toplot$mutbins /10000 )*10000
#
# # filter(Y,Fitness==max(Fitness))
# # filter(Y,Fitness>max(Fitness)-10000)
# # filter(Y,sample.ID==6169)
# # filter(Y,sample.ID==6169)
#
# myplot<-
# ggplot(toplot,aes(y=Fitnessnonzero,x=mutbins)) +
#   xlab('Number of mutations')+
#   ylab('Absolute lifetime fitness (# offspring)') +
#   # ylab('Relative lifetime fitness (# offspring / max)') +
# #   geom_hex(alpha=0.8)+
# # 	scale_fill_gradientn(colours = moiR::jet.colors(),trans = "log")+
# 	geom_point(col=transparent('grey50'), aes(y=Fitness,x=mutations) , shape=16)+
#   stat_smooth(data=filter(toplot, mutations> 1 & oFitness !=0))
#   stat_summary(data=filter(toplot, mutations> 1 & oFitness !=0),geom = 'line',
#                fun.y = "max" )
# myplot


# stat_smooth(data=dplyr::filter(toplot, oFitness>0),
  #              aes(y=oFitness,x=mutations),
  #             formula = y ~ x,
  #             # formula = y ~ I(x^2)+x,
  #             se=TRUE,lty='dashed',col='black')
  # # stat_smooth(data=filter(toplot, oFitness>0),
  #             aes(y=oFitness,x=mutations),
  #             se=TRUE,lty='dashed',col='black')
# myplot
# save_plot(myplot,filename = 'figs/kondrasov_nmut_fitness.pdf',base_height = 6,base_width = 6)
# saveRDS(file='figsfiles/kondrasov.rda',myplot)
