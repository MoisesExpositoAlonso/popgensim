
############ funcitons
mafsim<-function(p, type='uniform',rate=1, cutoff=1){
  stopifnot(type %in% c('uniform','exponential'))

  if(type=='uniform'){
    runif(p,0,0.5)
  }else if(type=='exponential'){
    es<-rexp(p,rate = 1)
    # es/max(es)
    0.5* es/ (cutoff*max(es))
  }
}

Xsim<-function(n=100,p=1000, maf){
  # X<-cbind(sapply(1:p,function(i) rbinom(n = n,size=1,prob = maf[i])))
  X<-cbind(sapply(1:p,function(i) sample(c(-1,+1),size=n,replace=T,prob = c(1-maf[i],maf[i] ) )))
return(X)
}


ssim<-function(p,svar=0.1){
  PropoS(1000,svar)
}

wsim<-function(X,s,mode=1, epi=1,mu=1){
  wC(X,s,mode,epi,mu)
}


sampleEys<-function(Eys,a,b,p,rep=1){
  Yobs<-c()
  for(i in 1:length(Eys)){
    Yobs<-c(Yobs,rnorm(rep,Eys[i], abs(a+(Eys[i]*b)) ))
  }
  Yobs[Yobs<0] <-0
  if(p!=0){
    Yobs[sample(1:length(Yobs),size = ceiling(p*length(Yobs)) ) ]<-0
  }
return(Yobs)
}




#####**********************************************************************#####
multievo<-function(p,
                   n,
                   rate,
                   svar,
                   a,
                   b,
                   mu,
                   mode,
                   No,
                   Nmax,
                   d,
                   tmax,
                   Replicates=2){

    pop<-lapply(1:Replicates,function(i){
                maf=mafsim(p,'exponential',1)
                X=Xsim(n,p,maf)
                s=PropoS(p,svar)
                w=  wC(X,s,mode,1,mu)

                po<-multigenpopsimC(fitness=w,
                                     No=No,
                                     Nmax=Nmax,
                                     a=a,
                                     b=b,
                                     t=tmax,
                                     d=d
                                     )

                return(po)
      })

  return(pop)
}




#####**********************************************************************#####

popsizeplot<-function(obj){

  dat_<-sapply(obj,function(x) apply(as.matrix(x),2,sum) )
  dat<-data.frame(popsize=as.numeric(dat_), generations=1:nrow(dat_), replicate=sort(rep(1:ncol(dat_), nrow(dat_)))  )

spline_int <- as.data.frame(spline(dat$popsize~ dat$generations))

p<-ggplot(dat) +
  geom_line(aes(y=popsize,x=generations,group=replicate), color='grey')+
  geom_line(data = spline_int, aes(y=y,x=x),lwd=2)+
      ylab('N(t)') +
      xlab("Generations")
  return(p)
}


####  Allele frequencies
allelefreqplot<-function(obj){

  dat_<-sapply(obj,function(x) apply(as.matrix(x),2,sum) )
  dat<-data.frame(popsize=as.numeric(dat_), generations=1:nrow(dat_), replicate=sort(rep(1:ncol(dat_), nrow(dat_)))  )

spline_int <- as.data.frame(spline(dat$popsize~ dat$generations))

p<-ggplot(dat) +
  geom_line(aes(y=popsize,x=generations,group=replicate), color='grey')+
  geom_line(data = spline_int, aes(y=y,x=x),lwd=2)
  # stat_summary(fun.y = mean, geom="line", aes(y=popsize,x=generations))
  # stat_summary(aes(y=popsize,x=generations))
  # stat_smooth( method = lm, formula = y ~ poly(x, nrow(dat_)-1 ), aes(y=popsize,x=generations))
  # stat_smooth(aes(y=popsize,x=generations), se=F,col='black',span=10)

  # p<-ggplot() +
  #     ylim(c(0,max(as.numeric(as.matrix(dat))))) +
  #     xlim(c(1, input$tmax))+
  #     ylab('N(t)') +
  #     xlab("Generations")
  # for(rep in 1:ncol(dat)){
  #   newdat<-data.frame(y=dat[,rep],
  #                      x=1:input$tmax)
  #   p <- p+ geom_line(data=newdat,aes(y=y,x=x))
  # }
  return(p)
}
