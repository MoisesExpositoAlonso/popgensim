################################################################################
### Real data associate + inference from LD
################################################################################

################################################################################
message("Note: library GWS, run devtools::load_all if code has been modified")

# devtools::load_all('.')
suppressMessages(library(ggplot2))
suppressMessages(library(latex2exp))
suppressMessages(library(cowplot))
library(gws)
################################################################################
## Define genome matrix
message("Loading or simulating a genome matrix")

X<-gws::get_topSNPs(500,navalue=0.5)
dim(X)
X[1:5,1:10]

################################################################################
# Phenotype

Yl<-list(
  Y=get_realfit(relative = F),
  eff=1:ncol(X)
)



### How good the selection estimates predict the individual fitness
message("Estimating SNP effects using different models...")

gwatime<-system.time(
  res<-data.frame( true=Yl$eff,
                   bgwa=gwa(X,Yl$Y),
                   b=cgwa(X,Yl$Y),
                   bridge=ridgegwa(X,Yl$Y,type='penalized'),
                   blasso=lassogwa(X,Yl$Y)
  )
)
message("\t... running time was ",round(gwatime[3],digits = 3), " seconds or ", round(gwatime[3]/60),digits = 3, " minutes or ", round(gwatime[3]/60/60,digits = 3), " hours.")

p1<-(ggdotscolor(x=res$true,y=res$bgwa, xlab='True effect',ylab='Inferred effect',color='grey70')  )+ggtitle("Marginal")
p2<-(ggdotscolor(x=res$true,y=res$b, xlab='True effect',ylab='Inferred effect',color='grey70')  )+ggtitle("Conditional")
p3<-(ggdotscolor(x=res$true,y=res$bridge, xlab='True effect',ylab='Inferred effect',color='grey70'))+ggtitle("Ridge")
p4<-(ggdotscolor(x=res$true,y=res$blasso, xlab='True effect',ylab='Inferred effect',color='grey70')  )+ggtitle("Lasso")

panel<-plot_grid(p1,p2,p3,p4)
panel
save_plot(filename = paste('method/figs/assoc_','real',".pdf"),plot = panel,
          base_height = 7,base_width = 7)


# mean((abs(res$bgwa)-abs(res$bridge)) / (abs(res$bgwa)))
# sum(abs(res$bgwa-res$true) / sum(abs(res$true)) )
# sum(abs(res$bgwa) / sum(abs(res$true)) )
# sum(abs(res$blasso) / sum(abs(res$true)) )
# sum(abs(res$bridge) / sum(abs(res$true)) )

# Characterize somehow the difference between linked and nonlined
linkedeffect <- mean((abs(res$bgwa)-abs(res$blasso)) / (abs(res$bgwa))) %>% format(.,digit=3)

plinked<-(ggdotscolor(y=res$bgwa,x=res$bridge, ylab='Marginal effect',xlab='Conditional effect',color='grey70') %>% addggregression(se=F,lty='dashed') )+
  ggtitle(paste("mean(b-bc / b) = ",linkedeffect)) +
  geom_abline(slope=1,intercept = 0,col=('grey70'))+
  # geom_vline(xintercept = 0,col=('grey70'))+
  # geom_hline(yintercept = 0,col=('grey70'))+
  coord_cartesian(xlim=range(c(res$bgwa,res$bridge)),ylim=range(c(res$bgwa,res$bridge)) )
plinked
save_plot(filename = paste('method/figs/assoc_linkedeffect_','real',".pdf"),plot = plinked,
          base_height = 5,base_width = 5)

indirecteffect<-(ggdotscolor(x=res$true,y=res$bgwa-res$bridge, xlab='Position in genome',ylab=TeX('\beta_{i} x r^2_{ij}'),color='grey70'))+
  geom_line(col='grey50')+geom_hline(yintercept = 0)
  # ylab('difference'))
  # ylab()
indirecteffect

save_plot(filename = paste('method/figs/assoc_linkedeffect_difference','real',".pdf"),plot = indirecteffect,
          base_height = 2,base_width = 5)


################################################################################
### LD
################################################################################

Rcpp::sourceCpp('src/ld.cpp')

Yrel<-relative(Yl$Y)
Yrel[is.na(Yrel)]<-1
hist(Yrel)

## Visualization plot of LD change style old R2

pre_D<-ldCnext_D(X,Yrel,Dprior = TRUE)
pre_r2<-ldC(meanvarcent.mat(X))
obs_Dprime<-ldCnext_D(X,Yrel,Dprime = TRUE)
obs_Ddiff<-ldCnext_D(X,Yrel,Ddiff = TRUE)
obs_r2<-ldC(meanvarcent.mat(X_fit(X ,round(Yl$Y)))) # slow

lm(fn(obs_r2) ~ fn(pre_D)) %>% summary
summary(obs_r2)

pdf(paste('method/figs/ldobs_','r2',".pdf"),height = 5,width = 5,useDingbats = FALSE)
ggdotscolor(y=upperTmat(obs_r2) , x= upperTmat(pre_r2),color='grey70',
            xlab= TeX('r^{2} before selection'),
            ylab= TeX('r^{2} after selection')
            ) %>% addggregression(lty='dashed',col='darkred') +
  geom_abline(intercept = 0,slope = 1)+
  xlim(c(-0.2,+1))+
  ylim(c(-0.2,+1))
dev.off()

ggdotscolor(y=upperTmat(obs_Dprime) , x= upperTmat(pre_D),color='grey70',
            xlab= TeX('D before selection'),
            ylab= TeX('D after selection')
) %>% addggregression(lty='dashed',col='darkred') +
  geom_abline(intercept = 0,slope = 1)

## Visualization of R from Felsenstein with random

obs<-ldCnext(X,Yrel,R=TRUE,dolog = FALSE) %>% upperTmat()
# random<-ldCnext(X,sample(Yrel),R=TRUE,dolog = FALSE) %>% upperTmat()
# hist(obs-random)
random<-ldCnext_perm(X,sample(Yrel),R=TRUE,diff = FALSE,dolog = FALSE)
random.mean<-apply(random,1,mean)
random.95<-apply(random,1,function(i) quantile(i,p=0.95))
hist(random.mean)
summary(random.mean)
hist(obs)

ggplot(data.frame(Rexcess=obs - random.mean),aes(x=Rexcess))+geom_histogram(fill='grey30',color='white') +
  xlab( 'R excess')+
  geom_vline(xintercept = 0,col='skyblue3',lwd=1.5,lty='dashed')
summary(obs - random.95)
summary(obs - random.mean)
summary(obs - random.mean)

pdf(paste('method/figs/ldobs_','real',".pdf"),height = 5,width = 5,useDingbats = FALSE)
qplot(x=obs,main='',xlab='R after selection (w1w4+1/w2w3+1)')
qplot(x=random.mean,main='',xlab='R after selection random')
ggplot(data.frame(Rexcess=obs - random.mean),aes(x=Rexcess))+geom_histogram(fill='grey30',color='white') +
  xlab( 'R excess')+
  geom_vline(xintercept = 0,col='skyblue3',lwd=1.5,lty='dashed')
dev.off()

mean(obs)
mean(random)

## Expected

qplot(x=obs,main='',xlab='R after selection (w1w4+1/w2w3+1)')


# Visualize the matrix
require(RColorBrewer)
tiny<-ldCnext(X,Yrel,R=TRUE,dolog = FALSE)
# tiny<-tiny[1:5,1:5]

matrixchange<-ggplot(data.frame(x=fn(row(tiny)),y=fn(col(tiny)),varcol = fn(tiny)))+
  geom_tile(aes(x,y,fill=varcol))+
  scale_fill_gradientn("",colours=palette(jet.colors())) #+
  # theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #       axis.text.y=element_blank(),axis.ticks=element_blank(),
  #       axis.title.x=element_blank(),
  #       axis.title.y=element_blank(),legend.position="none",
  #       panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
  #       panel.grid.minor=element_blank(),plot.background=element_blank())
matrixchange
save_plot(filename=paste0('method/figs/ldchange_matrix','real',".pdf"),
          base_height =  5,base_width = 5,useDingbats = FALSE,
          matrixchange )

matrixchange<-ggplot(data.frame(x=fn(row(pre_r2)),y=fn(col(pre_r2)),varcol = fn(pre_r2)))+
  geom_tile(aes(x,y,fill=varcol))+
  scale_fill_gradientn("",colours=palette(jet.colors()))

### Model of power epistasis
sumdifsq<-function(e,y){
 sum( (y - upperTmat(ldCexpect(s = res$bridge, e)) )^2 )
}
res_power<-optim(par=c(e=1),fn =sumdifsq,y=obs,method = 'BFGS')
res_power


# # power with a bit of noise
# sumdifsq_noise<-function(pars,y){
#   sum( (y - (upperTmat(ldCexpect(s = res$bridge+ rnorm(length(res$bridge),pars[2],pars[3]), pars[1])) ) )^2 )
# }
# res_power_noise<-optim(par=c(e=1,mu=.1,sd=.1),fn =sumdifsq_noise,y=obs,method = 'BFGS')
# res_power_noise
# ggplot(data=
#          # data.frame(exp=upperTmat(ldCexpect(s = res$bridge, e=res_power_noise$par[1]))+ rnorm(length(obs),res_power_noise$par[2],res_power_noise$par[3]) ,
#           data.frame(exp=upperTmat(ldCexpect(s = res$bridge *rnorm(length(res$bridge),res_power_noise$par[2],res_power_noise$par[3]) , e=res_power_noise$par[1])) ,
#                                obs=obs
#          )) +
#   geom_histogram(aes(x=obs),fill=transparent('black',0.8))+
#   geom_histogram(aes(x=exp),fill=transparent('navy',0.8))



### Model of scalar epistasis
sumdifsq_scalar<-function(e,y){
  sum( (y - upperTmat(ldCexpectscalar(s = res$bridge, e)) )^2 )
}
res_scalar<-optim(par=c(e=1),fn =sumdifsq_scalar,y=obs,method = 'BFGS')
res_scalar


compare<-ggplot(data=
              data.frame(exp=upperTmat(ldCexpectscalar(s = res$bridge, e=res_scalar$par)),
                         exppow=upperTmat(ldCexpect(s = res$bridge, e=res_power$par)),
                         obs=obs
              )) +
  geom_histogram(aes(x=obs),fill=transparent('black',0.8),col=transparent('white'))+
  geom_histogram(aes(x=exp),fill=transparent('navy',0.4),col=transparent('white'))+
  geom_histogram(aes(x=exppow),fill=transparent('darkred',0.4),col=transparent('white'))+
  xlim(c(0.8,+1.2))+
  coord_cartesian(ylim=c(0,1500))
compare
save_plot(filename=paste('method/figs/ld_e_heuristix_','real',".pdf"),base_height =  5,base_width = 5,useDingbats = FALSE,
          compare )

################################################################################
## Inference
message("Estimating SNP effects using LD matrix change...")

message("\t... starting optimization ...")

# optimtime<-system.time(
#   opres<-optim(par = c(1.5,rnorm(ncol(X),0,1)),
#                fn = ld_Gauss_lik,
#                y=obs,
#                method="L-BFGS-B",
#                lower =c(-2, rep(-1,ncol(X)) ),
#                upper =c(+2, rep(+1,ncol(X)) )
#   )
# )


message("\t... optimization time was ",round(optimtime[3]), " seconds or ", round(optimtime[3])/60, " minutes or ", round(optimtime[3])/60/60, " hours.")

#
# ldres<-qplot(opres$par[-1],x=Yl$eff, xlab='True effects',ylab='Inferred effects from LD',main=paste('Epistatic term = ',round(opres$par[1],3)))
# save_plot(filename = paste('method/figs/ldinf_',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),
#           plot = ldres,
#           base_height = 7,base_width = 7)



message("\t...done")


# sink() # end log file
# close(logfile)
