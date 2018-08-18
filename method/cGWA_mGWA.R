## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(coda)

# library(moiR)
load_all("~/moiR")
load_all('.')
library(Rcpp)


###############################################################################
## Genomematrix
################################################################################

genomes<-readRDS('databig/genome.rda')
G <- attachgenomes(genomes)

# Gi<-bigmemory::attach.big.matrix('databig/gi.desc')
# Gm<-bigmemory::attach.big.matrix('databig/gm.desc')
# Go<-bigmemory::attach.big.matrix('databig/go.desc')
# Go<-bigmemory::attach.big.matrix('databig/toy.desc')
# Go<-bigmemory::attach.big.matrix('databig/toy.desc')
God<-bigmemory::attach.big.matrix('databig/godou.desc')
Go<-bigmemory::attach.big.matrix('databig/gc.desc')

m=10
# topcols=1:10
# topcols<-top_cols(genomes,p = m)
## topcols=1:ncol(Go) ##< careful, this is genome-wide

################################################################################
## Fitness mlp
################################################################################
Y<-readRDS('data/Y.rds')
head(Y)


###############################################################################
## GWA
################################################################################
# load('.RData')
sourceCpp("gwa.cpp")

topcols<-1:ncol(Go)
training<-sample(1:nrow(Go), size= floor(0.85 * nrow(Go)))
test<-c(1:nrow(Go))[- training]

Ymean=My(Y$rFitness, Y$row)

### Marginal
bm<-BMgwa(A=God@address,y=(Ymean-1)[training] ,vars=topcols-1, training = training-1,type=1)
newtopcols<- which(rank(bm) < 1001)

### Conditionals

# bm<-BMgwa(A=God@address,y=(Ymean-1)[training] ,vars=newtopcols-1, training = training-1,type=1)
brid<-BMgwa(A=God@address,y=(Ymean-1)[training] ,vars=newtopcols-1, training = training-1,type=3)
blas<-BMgwa(A=God@address,y=(Ymean-1)[training] ,vars=newtopcols-1, training = training-1,type=4)

ym<-BMprod(God@address,bm,test-1,newtopcols-1)+1
yrid<-BMprod(God@address,brid,test-1,newtopcols-1)+1
ylas<-BMprod(God@address,blas,test-1,newtopcols-1)+1

#### predictions
pa<-qplot(y=ym , x= Ymean[test], geom='point', plot=FALSE,
	ylab='Inferred fitness (# offspring x survival)',
	xlab='True fitness (# offspring x survival)')+
	geom_point(shape=19, alpha=0.2)+
	geom_hline(yintercept = 0,lty='dotted')+
	geom_vline(xintercept = 0,lty='dotted')+
  stat_smooth(method = 'glm',se=F,col='darkgrey',lty='dashed')
pb<-qplot(y=yrid , x= Ymean[test], geom='point', plot=FALSE,
	ylab='Inferred fitness (# offspring x survival)',
	xlab='True fitness (# offspring x survival)')+
	geom_point(shape=19, alpha=0.2)+
	geom_hline(yintercept = 0,lty='dotted')+
	geom_vline(xintercept = 0,lty='dotted')+
  stat_smooth(method = 'glm',se=F,col='darkgrey',lty='dashed')
pc<-qplot(y=ylas , x= Ymean[test], geom='point', plot=FALSE,
	ylab='Inferred fitness (# offspring x survival)',
	xlab='True fitness (# offspring x survival)')+
	geom_point(shape=19, alpha=0.2)+
	geom_hline(yintercept = 0,lty='dotted')+
	geom_vline(xintercept = 0,lty='dotted')+
  stat_smooth(method = 'glm',se=F,col='darkgrey',lty='dashed')

pb<-
  (ggdotscolor(y=yrid , x= Ymean[test],
	ylab='Inferred fitness (# offspring x survival)',
	xlab='True fitness (# offspring x survival)')+
	geom_point(shape=19, alpha=0.2)+
	geom_hline(yintercept = 0,lty='dotted')+
	geom_vline(xintercept = 0,lty='dotted') )%>% addggregression(se=F,col='darkgrey',lty='dashed')
pb
save_plot(filename = 'figs/fitprediction_gwa_real.pdf',plot=pb,base_height = 4,base_width = 4)


#### Manhattan

res<-cbind(map,bm)
# res$bc<-NA
res$blas<-NA
res$brid<-NA
res[newtopcols,'blas']<-blas
res[newtopcols,'brid']<-brid

head(res)
mhdata<-ggmanhattan_preparedata(res,stat.col='bm')

p<- ggplot(data=filter(mhdata$data, !is.na(brid)) )+
					 #, aes(y=brid,x=cumpos,group=factor(chr),color=factor(chr),size=0.1)
					 # ) +
  geom_point(aes(y=stat,x=cumpos,size=0.1),
  					 alpha=0.8,shape=16,color='black')+
# 	geom_point(aes(y=bc,x=cumpos,group=factor(chr),size=0.1),
#   					 alpha=0.8,shape=19,color='darkred')+
	geom_point(aes(y=brid,x=cumpos,group=factor(chr),size=0.1),
  					 alpha=0.8,shape=16,color='blue')+
	geom_point(aes(y=blas,x=cumpos,group=factor(chr),size=0.1),
  					 alpha=0.8,shape=16,color='darkgreen')+
  ylab(TeX('$\\beta $ (# offspring x survival)'))+
  xlab('Position (Mb)')+
	geom_hline(yintercept=0)+
	# ylim(c(-3,+3))+
  scale_x_continuous(breaks=mhdata$breaks, labels= mhdata$labels )+
  scale_color_manual(values= c("black","darkgrey", "black", "darkgrey" ,"black"),guide=FALSE)+
  scale_size_continuous(range=c(0,2), guide=FALSE)
p


p2<- ggplot(data=filter(mhdata$data, !is.na(brid)) )+
	    geom_density(aes(x=stat),fill=transparent('black'),bw=0.05)+
      geom_density(aes(x=brid),fill=transparent('blue'))+
      geom_density(aes(x=blas),fill=transparent('darkgreen'))+
      # xlab(TeX('$\\beta $ (# offspring x survival)'))+
      xlab("")+
      coord_flip()
p2

panel2<-plot_grid(p,p2,rel_widths = c(4,1))
save_plot(filename = 'figs/manhattan_gwa_real.pdf',plot=panel2,base_height = 3,base_width = 10)

# geom_point(aes(y=stat,x=cumpos,size=0.1),
#   					 alpha=0.8,shape=19,color='black')+
# # 	geom_point(aes(y=bc,x=cumpos,group=factor(chr),size=0.1),
# #   					 alpha=0.8,shape=19,color='darkred')+
# 	geom_point(aes(y=brid,x=cumpos,group=factor(chr),size=0.1),
#   					 alpha=0.8,shape=19,color='blue')+
# 	geom_point(aes(y=blas,x=cumpos,group=factor(chr),size=0.1),
#   					 alpha=0.8,shape=19,color='darkgreen')+
#   ylab("beta")+ xlab('Position (Mb)')+
# 	geom_hline(yintercept=0)+
# 	# ylim(c(-3,+3))+
#   scale_x_continuous(breaks=mhdata$breaks, labels= mhdata$labels )+
#   scale_color_manual(values= c("black","darkgrey", "black", "darkgrey" ,"black"),guide=FALSE)+
#   scale_size_continuous(range=c(0,2), guide=FALSE)
# p





# Go[,1:10]
# betam = BMmgwa(A = Go@address,
# 							 Y = Y$rFitness,
# 							 h = Y$row,
# 							 vars = topcols-1,
# 							 debug=TRUE)
# hist(betam)
# betam_dou = BMmgwa(A = God@address,
# 							 Y = Y$rFitness,
# 							 h = Y$row,
# 							 vars = topcols-1,
# 							 debug=TRUE)

# betam_dou = lm(My(Y$rFitness,Y$row) ~ Go[,topcols])


# fitm=BMprod(Go@address,betam, topcols-1)

# png('figs/gwa_marginal.png',height=3,width=5,useDingbats=FALSE)
# 	plot(betam, ylab=TeX('$\\beta_m$ (rel. fitness)'),
# 				xlab='SNP position',
# 				 pch=16)
# dev.off()


# pdf('figs/gwa_marginal.pdf',height=3,width=5,useDingbats=FALSE)
# plot(betam, ylab=TeX('$\\beta_m$ (rel. fitness)'),
# 			xlab='SNP position',cex=0,
# 			 pch=16)
# dev.off()

### Full Conditional
# betac= - BMcgwa(Go@address,Ymean-1,topcols-1)
# fitc=BMprod(Go@address,betac, topcols-1)

### Ridge
## (1) Get the best lambda with 100 SNPs
# betaridge= BMridge(Go@address,(Ymean-1) *(-1),topcols-1,l=seq(0.001,0.1,0.001))
# fitr=apply(betaridge$coef,2,function(i) BMprod(Go@address,betaridge$coef[,i], topcols-1) )
# ssr<-sapply(1:ncol(fitr),FUN = function(i) sum( (fitr[,i] - Ymean)^2) )
# plot(ssr, betaridge$lambda)

# (2) Compute coefficients for all SNPs

# betaridge= BMridge(Go@address,(Ymean-1) *(-1),topcols-1,l=c(1.0))
# betar<-betaridge$coef

# sourceCpp("gwa.cpp")
# betar= BMsimridge(Go@address,(Ymean-1) *(-1),topcols-1,l=c(1.0))
#
# fitr=BMprod(Go@address,betar, topcols-1)
#
# plot(Ymean,fitridge)


###############################################################################
## plots
################################################################################

p0<-ggdotscolor(
            betam ,x= betar,
						xlab=TeX("Conditional $\\beta$"),
						ylab=TeX("Marginal $\\beta$")
						) %>% addggregression(se=F,lty='dashed') #+ xlim(c(0,1)) +ylim(c(0,1))
p1<-ggdotscolor(fitm ,x= Ymean,
						xlab='True fitness',
						ylab='Inferred fitness (marginal)') %>%
						addggregression(se=F,lty='dashed')
p2<-ggdotscolor(fitr ,x= Ymean,
						xlab='True fitness',
						ylab='Inferred fitness (conditional)') %>%
						addggregression(se=F,lty='dashed')
save_plot(plot=plot_grid(p0,p1,p2,ncol=2), filename = 'figs/gwa_fitness_oriented.pdf',base_height = 8,base_width = 8)

saveRDS(file='figsfiles/gwa.rda',list(p0,p1,p2))

# save('gwasession.rda')
# saveRDS(file='databig/betam.rda',betam)
# saveRDS(file='databig/betac.rda',betac)
# saveRDS(file='databig/betar.rda',betar)
