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
gwares<-readRDS('dataint/rFitness_mlp.rda')
maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))

subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[501]),]
head(subg)
subg$rs<-paste(subg$chr,subg$pos,sep='_')

###############################################################################
## Subset of Genomematrix
genomes<-readRDS('databig/genome.rda')
G <- attachgenomes(genomes)
G[1:5,1:5]
Map <- genomes$map
head(Map)
Map$rs<-paste(Map$chr,Map$physical.pos,sep='_')
Map$order<-1:nrow(Map)

# Merge with top hits
matched<-merge(Map, subg, by='rs')
dim(matched)
head(matched)


X=change01(inputena.mat(G[,matched$order],value = 0))
X=change11neg(inputena.mat(G[,matched$order],value = 0))
X[1:5,1:5]
dim(X)

################################################################################
## Get phenotypes
Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Fitness_mlp")], by.y="id",all.x=T)[,7]
# Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Fitness_mli")], by.y="id",all.x=T)[,7]
Yabs<-Y
Y=relative(Y)
# Y=meanvarcent(Y)
Y[is.na(Y)]<-1
Yrel<-Y

hist(Yabs)
hist(Yrel)
hist(rexp(500,0.04))

1/250

myprod=1
for(i in 1:10){
  myprod= myprod * (1+0.7)
}
myprod

250^(1/10)

rgamma(5,0.1,shape = 0.5)

################################################################################
### Run the likelihood optimization
library(Rcpp)
sourceCpp('lgws.cpp')

# lambda(X[,1:2],c(.1,.2))
# lambda(X[,1:4],c(.1,.2,-.3,.5))
#
# lnLse(Y,X[,1:4],c(.1,.2,-.3,.5))
# lnLse(Y,X[,1:4],c(.1,.2,-.3,-100))

# sigma(
# as.matrix(expand.grid(c(-1,1),c(-1,1),c(-1,1))),
# rnorm(3,0.5,1)
# )




Xcopy<-X
X<-Xcopy[,1:100]
X<-Xcopy[,1:50]
X<-Xcopy[,1:10]
X<-X* (-1)

mean(X)

Yl<-simupheno(X,meanN = 0,sdN = 5,heritability = 0.4)
Y<-Yl$Y
true<-Yl$eff
hist(Y)
hist(Yrel)
MASS::fitdistr(Yrel,densfun = 'exponential')
hist(rexp(100,1))

lambda(X[1,],s = matched$beta[1:ncol(X)])
lambda(X[2,],s = matched$beta[1:ncol(X)])
lambda(X[3,],s = matched$beta[1:ncol(X)])



startpar<-rep(0.001,ncol(X))
startpar<-rnorm(ncol(X),1,5)
startpar<-matched$beta[1:ncol(X)]
system.time(
  sel<-optim(par=startpar,
                 fn = lik.gws.exp,
                 # gr = gr.gws.exp,
                 datalist=list(Y=Y,X=X),
                 method = "L-BFGS-B",
                 upper = rep(5,ncol(X)),
                 lower = rep(-5,ncol(X))
  )
)
sel
# sourceCpp('lgws.cpp')
lik.gws.exp(startpar,list(Y=Y,X=X))
lik.gws.exp(true,list(Y=Y,X=X))

lik.gws.exp(matched$beta[1:ncol(X)],list(Y=Y,X=X))
# lambda(X[1,], startpar)

plot(sel$par~ startpar)
plot(sel$par~ true)
plot(sel$par~ matched$beta[1:ncol(X)], xlab='GWA',ylab='Likelihood infer')



system.time(
  soursel<-optim(par=startpar,
                 fn = lik.gws.exp,
                 # gr = gr.gws.exp,
                 datalist=list(Y=Y,X=X),
                 method = "L-BFGS-B",
                 upper = rep(5,ncol(X)),
                 lower = rep(-5,ncol(X))
  )
)



# res<-sigma(X[,1:2],c(.1,.2),2)
# (1+c(.1)* X[,1]) *  (1+c(.2)* X[,2])


#
# lnLse(Y,X,rnorm(ncol(X),0.5,0.1),2)
# lnLse(Y,X,rnorm(ncol(X),0.5,0.1),1)
# lnLse(Y,X,rnorm(ncol(X),0.2,0.1),2)
#
# lnLse(Y,X,rnorm(ncol(X),0,0.1),2)
# lnLse(Y,X,rnorm(ncol(X),2,0.1),2)
# lnLse(Y,X,matched$beta,2)


# lik.gws<-function (par=c(2,rep(0.001,ncol(X) )),datalist) {
#  lnLse(datalist$Y,datalist$X,par[-1],par[1],FALSE)
# }

# Xcopy<-X
# X<-Xcopy[,1:50]
# X<-Xcopy[,1:500]

# Real data small
# startpar<-c(2,rep(0.001,ncol(X)))
# system.time(
# soursel<-optim(par=startpar,
#                fn = lik.gws,datalist=list(Y=Y,X=X),
#                method = "L-BFGS-B"
# )
# )


# pstart<-ggdotscolor(soursel$par[-1], x=startpar[-1],xlab="Starting values",ylab="Lgws inferred")%>%addggregression()
# pLik<-ggdotscolor(soursel$par[-1], matched$beta,ylab="GWA",xlab="Lgws")%>%addggregression()
# save_plot("figs/LGWS_inference.pdf",plot = pLik)
#
# 1/sigma(X,soursel$par[-1],soursel$par[1])
#
# rgamma(500,1/sigma(X,soursel$par[-1],soursel$par[1]),1/sigma(X,soursel$par[-1],soursel$par[1])) %>% qplot
# rgamma(500,.3,.3) %>% qplot()
# qplot(Y, xlim=c(-1,10))
# rgamma(500,1/sigma(X,soursel$par[-1],soursel$par[1]),1/sigma(X,soursel$par[-1],soursel$par[1])) %>% qplot(.,xlim=c(-1,10))
#
#
#
# # Simulated data
# Ydata<-simupheno(X,sparsity = 0,type = 'epistatic',epistasis = 1)
# hist(Ydata$Y)
# Ydata$Y <- relative(Ydata$Y)
# hist(Ydata$Y)
# mean(Ydata$Y)
#
# # pdf('likelihood_inference.pdf')
# ressimu<-optim(par=c(1,rep(0.001,ncol(X))),
#                fn = lik.gws,datalist=list(Y=Ydata$Y,X=X),
#                method="L-BFGS-B"
#         )
# ggdotscolor(y=ressimu$par[-1], Ydata$eff, ylab = 'Inferred',xlab = 'True') %>% addggregression() +ggtitle('unguided')
# ressimu<-optim(par=c(1,Ydata$eff),
#                fn = lik.gws,datalist=list(Y=Ydata$Y,X=X),
#                method="L-BFGS-B"
# )
# ggdotscolor(y=ressimu$par[-1], Ydata$eff, ylab = 'Inferred',xlab = 'True') %>% addggregression() +ggtitle('guided start')
# ressimu<-optim(par=c(1,addnoise(Ydata$eff,0.2)),
#                fn = lik.gws,datalist=list(Y=Ydata$Y,tX=t(X)),
#                method="L-BFGS-B"
# )
# ggdotscolor(y=ressimu$par[-1], Ydata$eff, ylab = 'Inferred',xlab = 'True') %>% addggregression() +ggtitle('guided with 20% sd. noise')
# ressimu<-optim(par=c(1,addnoise(Ydata$eff,0.8)),
#                fn = lik.gws,datalist=list(Y=Ydata$Y,tX=t(X)),
#                method="L-BFGS-B"
# )
# ggdotscolor(y=ressimu$par[-1], Ydata$eff, ylab = 'Inferred',xlab = 'True') %>% addggregression() +ggtitle('guided with 80% sd. noise')
# dev.off()
