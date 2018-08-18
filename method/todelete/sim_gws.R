################################################################################
### Given genetic architecture: simulate + associate + inference from LD
################################################################################

suppressMessages(library(argparser) )

cat("\n")
p <- arg_parser("Simulation of fitness effects and inference via association and LD change")
p <- add_argument(p, "--genome", help="The subset of SNPs should be analysed",default="top")
p <- add_argument(p, "--sparsity", help="The degree of sparsity of selection, from 0 to 1",default=0.95)
p <- add_argument(p, "--phenotype", help="Additive or epistatic phenotype generation from positions",default="additive")
p <- add_argument(p, "--heritability", help="The degree of heretability of selection, from 0 to 1",default=1)
p <- add_argument(p, "--e", help="The epistatic term in case the phenotype is epistatic",default=1)
argv<-parse_args(p)

# logfile<-file(paste0('method/logs/',argv$genome,"_sparse",argv$sparsity,'_',argv$phenotype ,"_h",argv$heritability,".log"),open="wt")
# sink(logfile)
# sink(logfile, type = "message") # opening a file for log


cat(paste("\n ----------------------------------------------------- \n"))
cat("\n Interpreting arguments:\n")
  cat(paste("\n  --> SNPs subset: ",argv$genome,"\n"))
  cat(paste("\n  --> Sparcity: ",argv$sparsity,"\n"))
  cat(paste("\n  --> Additivity/epistasis: ",argv$phenotype,"\n"))
  cat(paste("\n  --> Heritability: ",argv$heritability,"\n"))
cat(paste("\n ----------------------------------------------------- \n"))


################################################################################
# devtools::load_all('.')
message("Note: library GWS, run devtools::load_all if code has been modified")
suppressMessages( # to avoid all the loading messages
  library('gws')

  )
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

################################################################################
## Define genome matrix
message("Loading or simulating a genome matrix")
X<-get_genomematrix(type =argv$genome)
dim(X)

################################################################################
# How well models predict true effects
### sparse vs gaussian & additive vs epistatic
message("Simulating SNP effects and phenotypes")
Yl<-simupheno(X,
              sparsity = argv$sparsity,
              type = argv$phenotype,
              epistasis = argv$e,
              heritability = argv$heritability
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

p1<-(ggdotscolor(x=res$true,y=res$bgwa, xlab='True effect',ylab='Inferred effect',color='grey70') %>% addggregression(se=F,lty='dashed') )+ggtitle("Marginal")
p2<-(ggdotscolor(x=res$true,y=res$b, xlab='True effect',ylab='Inferred effect',color='grey70') %>% addggregression(se=F,lty='dashed') )+ggtitle("Conditional")
p3<-(ggdotscolor(x=res$true,y=res$bridge, xlab='True effect',ylab='Inferred effect',color='grey70') %>% addggregression(se=F,lty='dashed') )+ggtitle("Ridge")
p4<-(ggdotscolor(x=res$true,y=res$blasso, xlab='True effect',ylab='Inferred effect',color='grey70') %>% addggregression(se=F,lty='dashed') )+ggtitle("Lasso")

panel<-plot_grid(p1,p2,p3,p4)
panel
save_plot(filename = paste0('method/figs/assoc_',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),plot = panel,
          base_height = 7,base_width = 7)


mean((abs(res$bgwa)-abs(res$bridge)) / (abs(res$bgwa)))
sum(abs(res$bgwa-res$true) / sum(abs(res$true)) )
sum(abs(res$bgwa) / sum(abs(res$true)) )
sum(abs(res$blasso) / sum(abs(res$true)) )
sum(abs(res$bridge) / sum(abs(res$true)) )

# Characterize somehow the difference between linked and nonlined
linkedeffect <- mean(abs(res$bgwa-res$true) / abs(res$bgwa) )

plinked<-(ggdotscolor(x=res$bgwa,y=res$bridge, xlab='Marginal effect',ylab='Conditional effect',color='grey70') %>% addggregression(se=F,lty='dashed') )+
  ggtitle(paste("mean(b-bc / b) = ",linkedeffect)) +
  geom_abline(slope=1,intercept = 0,col=('grey70'))+
  coord_cartesian(xlim=range(c(res$bgwa,res$bridge)),ylim=range(c(res$bgwa,res$bridge)) )
plinked
save_plot(filename = paste0('method/figs/assoc_linkedeffect_',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),plot = plinked,
          base_height = 5,base_width = 5)



################################################################################
## Visualization plot of LD change


Yrel<-relative(Yl$Y)
Yrel[is.na(Yrel)]<-1
hist(Yrel)


obs<-ldCnext(X,Yrel,R=TRUE,dolog = FALSE) %>% upperTmat()
# random<-ldCnext(X,sample(Yrel),R=TRUE,dolog = FALSE) %>% upperTmat()
random<-ldCnext_perm(X,sample(Yrel),R=TRUE,dolog = FALSE)
random.mean<-apply(random,1,mean)

pdf(paste0('method/figs/ldobs_',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),height = 7,width = 7,useDingbats = FALSE)
qplot(x=obs,main='',xlab='R after selection (w1w4+1/w2w3+1)')
qplot(x=random.mean,main='',xlab='R after selection random')
qplot(obs - random.mean,xlab='R excess')
dev.off()

# Visualize the matrix
require(RColorBrewer)
tiny<-ldCnext(X,Yrel,R=TRUE,dolog = FALSE)
# tiny<-tiny[1:5,1:5]

matrixchange<-ggplot(data.frame(x=fn(row(tiny)),y=fn(col(tiny)),varcol = fn(tiny)))+
  geom_tile(aes(x,y,fill=varcol))+
  scale_fill_gradientn("",colours=palette(jet.colors()))+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
save_plot(filename=paste0('method/figs/ldchange_matrix',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),
          base_height =  5,base_width = 5,useDingbats = FALSE,
          matrixchange )


## Inference



### Model of power epistasis
sumdifsq<-function(e,y){
  sum( (y - upperTmat(ldCexpect(s = res$bridge, e)) )^2 )
}
res_power<-optim(par=c(e=1),fn =sumdifsq,y=obs,method = 'BFGS')
res_power


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
save_plot(filename=paste0('method/figs/ldobs_',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),
          base_height =  5,base_width = 5,useDingbats = FALSE,
          compare )

#
# message("Estimating SNP effects using LD matrix change...")
#
# message("\t... starting optimization ...")
#
# optimtime<-system.time(
#                         opres<-optim(par = c(1.5,rnorm(ncol(X),0,1)),
#                               fn = ld_Gauss_lik,
#                               y=obs,
#                               method="L-BFGS-B",
#                               lower =c(-2, rep(-1,ncol(X)) ),
#                               upper =c(+2, rep(+1,ncol(X)) )
#                         )
# )
#
#
# message("\t... optimization time was ",round(optimtime[3]), " seconds or ", round(optimtime[3])/60, " minutes or ", round(optimtime[3])/60/60, " hours.")
#
#
# ldres<-qplot(opres$par[-1],x=Yl$eff, xlab='True effects',ylab='Inferred effects from LD',main=paste('Epistatic term = ',round(opres$par[1],3)))
# save_plot(filename = paste('method/figs/ldinf_',argv$genome,"_",argv$sparsity,'_',argv$phenotype ,"_",argv$heritability,".pdf"),
#           plot = ldres,
#           base_height = 7,base_width = 7)
#
#

message("\t...done")


# sink() # end log file
# close(logfile)
