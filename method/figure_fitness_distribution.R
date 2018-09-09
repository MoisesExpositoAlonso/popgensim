#### real data ####
## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(moiR)
load_all('.')


Y<-readRDS('data/Y.rds')

qplot(Y$rFitness,
	xlab='Relative fitness'
	)

# Calculate mean per site
means<-field.s %>%
  group_by(site, water,indpop) %>% summarize(meanfit= mean(Fitness))

ml<-filter(means,site=='madrid', water=='l',indpop=='p')$meanfit
mh<-filter(means,site=='madrid', water=='h',indpop=='p')$meanfit

# Make relative
ff<-field.s %>%
  filter(site=='madrid', indpop=='p') %>%
  mutate(rFitness= ifelse(water =='l',
                          Fitness/ml,
                          Fitness/mh
                          ))

head(ff$rFitness)

fitplot<-ff %>%
  ggplot(.) +
  geom_histogram(aes(x = rFitness,fill=water),bins = 50)+
  ylab('Number of genotypes') +
  xlab('Mean-scaled fitness (survival x # seeds)')+
  # facet_wrap(~site) +
  scale_fill_manual(values = transparent(watercolors(),alpha=0.8 ) )+
  geom_vline(xintercept = 1,lty='dashed')
fitplot

save_plot(filename = 'figs/fitness_distribution_real_waterdrought_madrid.pdf',
          plot = fitplot,base_width=4.2,base_height = 3.9)

####************************************************************************####
#### with simulated data ####

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


genomes<-readRDS('databig/genome.rda')
G <- attachgenomes(genomes)

Go<-bigmemory::attach.big.matrix('databig/go.desc')

m=100
# topcols=1:10
top=sample(1:ncol(Go),m)
# topcols<-top_cols(genomes,p = m)
## topcols=1:ncol(Go) ##< careful, this is genome-wide


###########################
## Example of distribution of selection coefficients
library(Rcpp)
sourceCpp("MCMCmean.cpp")

slist<- cbind(
  fn(PropoS(1e5,svar = 0.05)),
  fn(PropoS(1e5,svar = 0.1)),
  fn(PropoS(1e5,svar = .3))
) %>% as.data.frame()

head(slist)
p1<-ggplot(slist) +
  # geom_histogram(aes(x=V3),fill=transparent("grey80"))+
  # geom_histogram(aes(x=V2),fill=transparent("grey50")) +
  # geom_histogram(aes(x=V1),fill=transparent("grey20"))
  geom_histogram(aes(x=V3),fill=transparent("#e41a1c"))+
  geom_histogram(aes(x=V2),fill=transparent("#4daf4a")) +
  geom_histogram(aes(x=V1),fill=transparent("#377eb8"))
p1

p2<-ggplot(slist) +
  # geom_histogram(aes(x=log(1+V3)),fill=transparent("grey80"))+
  # geom_histogram(aes(x=log(1+V2)),fill=transparent("grey50"))+
  # geom_histogram(aes(x=log(1+V1)),fill=transparent("grey20"))
  geom_histogram(aes(x=log(1+V3)),fill=transparent("#e41a1c"))+
  geom_histogram(aes(x=log(1+V2)),fill=transparent("#4daf4a"))+
  geom_histogram(aes(x=log(1+V1)),fill=transparent("#377eb8"))
p2


## Example of distribution of fitness for 10 and 100 SNPs under selection

genomes<-readRDS('databig/genome.rda')
# G <- attachgenomes(genomes)

Go<-bigmemory::attach.big.matrix('databig/gc.desc')

# ## Preparate all parameters
PI=0.5
BE=0.1
AI=0

sdat1<-simulate_data(Go=Go,p=0.2,b=0.3,a=0,ss=0,m=100,svar = 0.05)
sdat2<-simulate_data(Go=Go,p=0.2,b=0.3,a=0,ss=0,m=100,svar = 0.1)
sdat3<-simulate_data(Go=Go,p=0.2,b=0.3,a=0,ss=0,m=100,svar = 0.3)
simunum=2

toplot<- c(
  sdat1$Y$Y3,
  # sdat2$Y$Y3,
  sdat3$Y$Y3
)
toplot<-data.frame(x=toplot)
dim(toplot)
head(toplot)
toplot$group<-sort(rep(1:simunum,515*5))

# Plot te distributions
p3<-ggplot(toplot)+
  geom_histogram(aes(x=x,group=group,fill=factor(group)),
                 color=transparent("white",alpha = 0), alpha=0.5,
                 bins=60)+
  scale_fill_manual(values=c(
                             "#377eb8",
                             # "#4daf4a",
                             "#e41a1c"
                             ))+
  xlim(c(-0.1,5))+
  guides(fill=FALSE)
p3


# add example of 3 genotype means under the 3 models
genos<-c(446 ,510, 120)
p3<-p3 +
  geom_vline(xintercept = sdat1$Ey3[genos], col="#4daf4a") +
  geom_vline(xintercept = sdat2$Ey3[genos], col="#e41a1c") +
  geom_vline(xintercept = sdat3$Ey3[genos], col="#377eb8")


panel<-plot_grid(plot_grid(p1,p2,nrow=1),p3, nrow=2)

save_plot(filename = 'figs/fitness_distribution_example_simulation.pdf',
          plot = panel,
          base_width=6,base_height = 6)


# sourceCpp('MCMCmean.cpp')
#
# Genos<-Ey_go(Go[,top], s= runif(length(top),-1,+1) ,1)
# Genos
#
# G1<-rnorm(100,Genos[1],sd = Genos[1] *0.3)
# G2<-rnorm(100,Genos[513],sd = Genos[513] *0.3)
# G3<-rnorm(100,Genos[500],sd = Genos[500] *0.3)
#
# G1[sample(1:length(G1), 20)]<-0.
# G2[sample(1:length(G1), 25)]<-0.
# G3[sample(1:length(G1), 10)]<-0.
#
# hist(G1)
# hist(G2)
# hist(G3)
#
# example<-ggplot() +
#   geom_histogram(data=data.frame(x=G1),aes(x ),bins = 100,fill=transparent('#7fc97f'))+
#   geom_histogram(data=data.frame(x=G2),aes(x ),bins = 100,fill=transparent('#beaed4'))+
#   geom_histogram(data=data.frame(x=G3),aes(x ),bins = 100,fill=transparent('#fdc086'))+
#   ylab('Number of individuals') +
#   xlab('Mean-scaled fitness (survival x # seeds)')+
#   geom_vline(xintercept = 1,lty='dashed')
#
#
# save_plot(filename = 'figs/fitness_distribution_example_model.pdf',
#           plot = example,
#           base_width=4,base_height = 4)

