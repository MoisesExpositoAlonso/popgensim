
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
gwares<-readRDS('dataint/Fitness_mlp.rda')
gwares$order<-1:nrow(gwares)
maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))

subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[1000]),]
head(subg)
subg$rs<-paste(subg$chr,subg$pos,sep='_')

###############################################################################
## Subset of matrix
#data("genomes")
G <- genomes$genotypes
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
# Calculate change in LD in the area of interest

V1=ldC(meanvarcent.mat(X))
hist(V1)
V2=ldC(meanvarcent.mat(X_fit(X,Yabs)))

diag(V1) = NA
diag(V2) = NA

## Observe changes
# theselected2<-V2[100,]
# theselected1<-V1[100,]

ldchange<-qplot(fn(V1),fn(V2)) + geom_abline(intercept = 0,slope = 1) +
  stat_smooth()+
  stat_smooth(method = 'glm')+
  xlim(-1,1) + ylim(-1,1)+
  ylab(TeX("r^2_{t1}"))+
  xlab(TeX("r^2_{t0}"))
ldchange

save_plot(file="figs/ld_change.pdf",ldchange,base_height = 4,base_width = 4)


t.test(fn(V2)-fn(V1))
lm(fn(V2)~fn(V1)) %>% summary

################################################################################
# Generate the random LD change distribution

# V2emp=sapply(1:100,function(i){
#   Vtemp=ldC(meanvarcent.mat(X_fit(X.,sample(absfit))))
#   diag(Vtemp)<-NA
#   Vtemp
# })

#
# ldchange<-qplot(fn(V1),V2emp[,2],geom = "hex") + geom_abline(intercept = 0,slope = 1) +
#   stat_smooth()+
#   stat_smooth(method = 'glm')+
#   xlim(-1,1) + ylim(-1,1)+
#   ylab(TeX("r^2_{t+1}"))+
#   xlab(TeX("r^2_t"))
# ldchange
#
#
# regs<-sapply(1:100,function(i){
#   lm(V2emp[,i] ~fn(V1) ) %>% coefficients %>% tail(1)
# })
# thereg<-lm(fn(V2) ~fn(V1) ) %>% coefficients %>% tail(1)
#
# regs<-ggplot(data.frame(regs))+geom_density(fill=transparent("black"),col="white",aes(x=regs)) + geom_vline(xintercept = thereg,col="blue") + xlab("Regression coefficients")
# save_plot(file="figs/regs.pdf",base_height = 4,base_width = 4,regs)
#
# table(regs<thereg) / sum(table(regs<thereg))
#
# save_plot(file="figs/regs_ldchange.pdf",base_height = 4,base_width = 8,plot_grid(ldchange,regs,align = "h"))


# ## Visualize
# library(gaston)
# V1b<-V1
# V2b<-V2
# diag(V1b) = 1
# diag(V2b) = 1
# V2b[V2b>1]<-1
# V1b[V1b>1]<-1
#
# LD.plot( apply((V2b-V1b),2,abs),  write.ld = NULL ,pdf.file = 'figs/ld_matrixchange2.pdf',
#          color.scheme = function(ld) rgb(0.9 - abs(ld), 0.9 - abs(ld), 0.9 - abs(ld))
#         )

################################################################################
### V empirical based on linear combinations
sourceCpp('src/ld.cpp')

is.na(V1) %>% table
is.na(V) %>% table


cor.test(V,V1)
plot(V,V1)

Vt1<-ldCnow(X)
hist(Vt1)
# V[V>10]<-10
# V[V< (-10)]<-(-10)
# V[V == Inf]<-(10)
# V[V == -Inf]<-(-10)
Vt1[is.na(Vt1)]<-1
# V<-V[lower.tri(V)]

Vt12<-ldCnext(X,Y)
Vt12[is.na(Vt12)]<-1


rdif<-V2-V1
Rdif<-(Vt1*Vt12)-Vt1


cor(fn(rdif),fn(Rdif))
plot(fn(rdif),log10(fn(Rdif)+0.0001))


# V<-ldCnextmat(X,Y)


# Vemp=sapply(1:100,function(i){
#   Vtemp=ldC(meanvarcent.mat(X_fit(X.,sample(absfit))))
#   diag(Vtemp)<-NA
#   Vtemp
# })

