#### UTILITIES  ####

sampleEys<-function(Eys,a,b,p,rep=5){
  Yobs<-c( sapply(1:nrow(Go), function(i) rnorm(rep,Eys[i],a+(Eys[i]*b))))
  Yobs[Yobs<0] <-0
  if(p!=0) Yobs[sample(1:length(Yobs),size = ceiling(p*length(Yobs)) ) ]<-0
return(Yobs)
}

####************************************************************************####

##### CLASS #####
napMCMCres<-function(b,a,p,mu,epi,svar,ss,
                     newtopcols,training,testing,Y,Go,Ey1,Ey2,Ey3,s,ori){

  ### define objects
  value<-list(
              b=b,
              a=a,
              p=p,
              mu= mu,
              epi=epi,
              svar=svar,
              ss=ss,

              newtopcols=newtopcols,
              training=training,
              testing=testing,
              Y=Y,
              Go=Go,
              Ey1=Ey1,
              Ey2=Ey2,
              Ey3=Ey3,
              s=s,
              ori=ori,
              results=list()
            )
  class(value) <- append(class(value),"napMCMCres")

  return(value)
}

##### METHODS #####

###### (1) napMCMC #####
  runMCMC <- function(sdat,iterations=5e4,burnin,Simumode,Fitnessmode,Proposalmode, Priormode,doplot){
      UseMethod("runMCMC")
  }

  runMCMC.napMCMCres <- function(sdat,iterations=5e4,burnin=0.1,Simumode=1,Fitnessmode=1,Proposalmode=1, Priormode=1,doplot=TRUE){
    Si<-fc(Simumode)
    Fi<-fc(Fitnessmode)
    Pri<-fc(Priormode)
    Pro<-fc(Proposalmode)

    rows<-sdat$Y$row
    Ys<-sdat$Y[,Fitnessmode]

    # Run
sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["run"]]<-
          napMCMC(y = Ys ,
                  h = rows,
                  A=Go,
                  s=rep(0,m),
                  svar=0.1,
                  ss=0.1,
                  m=sdat$newtopcols-1,
                  n=sdat$training-1,
                  iterations=iterations,
                  Fitnessmode=Fitnessmode,
                  Priormode=Priormode,
                  Proposalmode=Proposalmode)

    # Summarize run
    sdat<-MCMCsummary(sdat,burnin=burnin,Simumode=Simumode,Fitnessmode=Fitnessmode,Priormode=Priormode,Proposalmode=Proposalmode)
    # Plot
    if(doplot) sdat<-splot(sdat,Simumode=Simumode,Fitnessmode=Fitnessmode,Priormode=Priormode,Proposalmode=Proposalmode)
    if(doplot) sdat<-iplot(sdat,Simumode=Simumode,Fitnessmode=Fitnessmode,Priormode=Priormode,Proposalmode=Proposalmode)
    if(doplot) sdat<-parameterplot(sdat,Simumode=Simumode,Fitnessmode=Fitnessmode,Priormode=Priormode,Proposalmode=Proposalmode)
    return(sdat)
  }


  # (1b) napMCMC clean
  MCMCsummary <- function(sdat,burnin,Simumode,Fitnessmode,Priormode, Proposalmode){
      UseMethod("MCMCsummary")
  }

  MCMCsummary.napMCMCres <- function(sdat,burnin=0.1,Simumode=1,Fitnessmode=1,Priormode=1, Proposalmode=1){
    Si<-fc(Simumode)
    Fi<-fc(Fitnessmode)
    Pri<-fc(Priormode)
    Pro<-fc(Proposalmode)

    schain<- sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["run"]]$chain
    burn= floor(burnin  * nrow(schain))

    # sinf<-apply(schain,2,median)
    sinf<-MCMCglmm::posterior.mode(as.mcmc(schain))
    sinf_range<-HPDinterval(as.mcmc(schain))
    # sinf<-apply(as.mcmc(schain[-c(1:burn),]),2,median)
    # sinf_range<-HPDinterval(as.mcmc(schain[-c(1:burn),]))

    sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["sinf"]] <-sinf
    sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["sinf_range"]] <-sinf_range

    return(sdat)
  }

##### # (2) runGWA #####
  runGWA <- function(sdat,gwatype,Simumode, doplot){
      UseMethod("runGWA")
  }
  runGWA.napMCMCres <- function(sdat,gwatype,Simumode=1,doplot=TRUE){

    rows<-sdat$Y$row
    Ys<-sdat$Y[,Simumode]

    sdat$results[[fc(Simumode)]][[paste0("GWA",gwatype)]][["run"]]  <-

      BMgwa(X0=sdat$Go, training=sdat$training,vars=sdat$newtopcols,
             yraw = My(Ys,rows),type = gwatype)

    # Plot
    if(doplot) sdat<-splot(sdat,Simumode=Simumode,gwatype=gwatype)
    if(doplot) sdat<-iplot(sdat,Simumode=Simumode,gwatype=gwatype)
    return(sdat)
  }

###### (3) plot selection coefficients #####
splot_ <-function(s,sinf,sinf_range=NULL){

    ## Bias and accuracy at the selsection level
    lmobj<-summary(lm(s~sinf))
    # lmobj<-summary(lm(sinf ~s))
    accuracy<-lmobj$r.squared %>% round(.,digits=3)
    if(dim(lmobj$coefficients)[1] ==1){
      bias<-"na"
    }else{
      bias<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
    }

    if( is.null(sinf_range)) sinf_range =cbind(sinf,sinf)

    ## PLot selection coefficients
psel<-qplot(x=s,y=sinf,
xlim=c(range(c(s,sinf_range))),
ylim=c(range(c(s,sinf_range)))
) +
  # geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )+
      geom_abline(intercept = 0,slope = 1,lty="dotted")+
      ylab("Inferred selection coefficients")+
      xlab("True selection coefficients")+
        ggtitle(TeX(paste("$R^2$ = ",accuracy, ", $\\beta = $",bias)))
    psel
    if(!is.null(sinf_range)) {
      psel<-psel+geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )
    }
  return(list(psel=psel,accuracy=accuracy,bias=bias))
}

splot <- function(sdat,Simumode, Fitnessmode,Priormode,Proposalmode,gwatype){
    UseMethod("splot")
}

splot.napMCMCres <- function(sdat,Simumode=1,
                               Fitnessmode=1,
                               Priormode=1,
                              Proposalmode=1,
                             gwatype=NULL){
    Si<-fc(Simumode)
    Fi<-fc(Fitnessmode)
    Pri<-fc(Priormode)
    Pro<-fc(Proposalmode)

  if(!is.null(gwatype)){ ## Plot GWA
    ## Extract
    s<-sdat$s * sdat$ori
    sinf<- sdat$results[[Si]][[paste0("GWA",gwatype)]][["run"]]$coefficients
    se<-sdat$results[[Si]][[paste0("GWA",gwatype)]][["run"]]$stderr
    sinf_upp<- sinf + (se*1.96)
    sinf_low<- sinf - (se*1.96)
    sinf_range<-cbind(sinf_low,sinf_upp)

    res<-splot_(s,sinf,sinf_range)

    sdat$results[[Si]] [[paste0("GWA",gwatype)]][["plot_s"]]<-res$psel
    sdat$results[[Si]] [[paste0("GWA",gwatype)]][["R2_B_s"]]<-c(res$accuracy, res$bias)

  }else{ ## Plot MCMC
    ## Extract
    s<-sdat$s * sdat$ori
    sinf<-  sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["sinf"]]
    sinf_range<- sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["sinf_range"]]

    res<-splot_(s,sinf,sinf_range)

   sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["plot_s"]]<-res$psel
   sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["R2_B_s"]]<-c(res$accuracy, res$bias)
  }
  return(sdat)
}


##### (4) plot individual fit #####

iplotreal<-function(true, inf){

  toplot<-data.frame(x=true, y=inf)
  toplot$xnozero <- toplot$x
  toplot$xnozero[toplot$xnozero==0 ] <- NA

  return(iplot_(toplot))
}

iplot_<-function(toplot){
    ## Bias and accuracy at the individual level
    lmobj2<-summary(lm(toplot$y~toplot$xnozero))
    # lmobj2<-summary(lm(toplot$xnozero ~toplot$y))
    accuracy2<-lmobj2$r.squared %>% round(.,digits=3)
    if(dim(lmobj2$coefficients)[1] ==1){
    bias2<-"na"
    }else{
    bias2<-coefficients(lmobj2)[2,1]  %>% round(.,digits=3)
    }


    pind<-ggplot(data=toplot,aes(y=y,x=x)) +
      geom_point(col="grey") +
      stat_smooth(aes(y=y , x=xnozero), se=F,
                  method="glm",lty="dashed",col="black")+
            ylim(range(c(toplot$x,toplot$y))) +
            xlim(range(c(toplot$x,toplot$y)))+
            ylab("Inferred individual fitness")+
            xlab("True individual fitness")+
        geom_abline(intercept = 0,slope = 1,lty="dotted")+
        geom_hline(yintercept = 0,lty="dotted")+
        geom_vline(xintercept = 0,lty="dotted")+
        ggtitle(TeX(paste("$R^2$ = ",accuracy2, ", $\\beta = $",bias2)))
  return(list(pind=pind,accuracy2=accuracy2,bias2=bias2))
}

iplot <- function(sdat,Simumode, Fitnessmode,Priormode,Proposalmode,gwatype){
    UseMethod("iplot")
}
iplot.napMCMCres <- function(sdat,Simumode=1,
                               Fitnessmode=1,
                              Priormode=1,Proposalmode=1,
                              gwatype=NULL){
    Si<-fc(Simumode)
    Fi<-fc(Fitnessmode)
    Pri<-fc(Priormode)
    Pro<-fc(Proposalmode)


  if(!is.null(gwatype)){ ## Plot GWA
    ## Extract
    trueY<-sdat[[paste0("Ey",Simumode)]] [sdat$testing,]
    sinf<- sdat$results[[Si]][[paste0("GWA",gwatype)]][["run"]]$coefficients
    intercept<- sdat$results[[Si]][[paste0("GWA",gwatype)]][["run"]]$meanasintercept

    infY<- ((sdat$Go[sdat$testing,sdat$newtopcols] %*%  sinf ) +1 ) *intercept
    # infY<- BMpred(sdat$Go@address,
    #               sinf,
    #               sdat$testing-1,
    #               sdat$newtopcols-1,
    #               intercept
    #               )

    ## PLot  individuals
    toplot<-data.frame(y= infY, x=trueY )
    toplot$xnozero<-toplot$x
    toplot$xnozero[toplot$xnozero==0]<-NA

    res<-iplot_(toplot)
    res
    print(res)

    sdat$results[[Si]][[paste0("GWA",gwatype)]][["plot_i"]]<-res$pind
    sdat$results[[Si]][[paste0("GWA",gwatype)]][["R2_B_i"]]<-c(res$accuracy2, res$bias2)

  }else{ ## Plot MCMC
    ## Extract
    trueY<-sdat[[paste0("Ey",Simumode)]] [sdat$testing,]

    sinf<- sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["sinf"]]
    sinf_range<-sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["sinf_range"]]

    infY<-Ey_go(sdat$Go[sdat$testing,sdat$newtopcols],
                s = sinf,mode = Fitnessmode)

    toplot<-data.frame(y= infY, x=trueY )
    toplot$xnozero<-toplot$x
    toplot$xnozero[toplot$xnozero==0]<-NA

    # print((toplot))
    res<-iplot_(toplot)
    res
    print(res)
    sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["plot_i"]]<-res$pind
    sdat$results[[Si]][[Fi]][[Pri]][[Pro]] [["R2_B_i"]]<-c(res$accuracy2, res$bias2)
  }

  return(sdat)
}

#### (0) Plot global parameters #####

parametersummary<-function(parchain,parnames){
  r<-as.mcmc(parchain)
  colnames(r) <- parnames
  summary(r)
}
parameterplot_<-function(parchain,pnames){

  colnames(parchain)<-pnames

  plotlist<-list()

  for(i in pnames){
    panel<-plot_grid(
      qplot(parchain[,i],
            x=1:nrow(parchain),
            # ylim=c(0,1),
            geom='line',xlab='iterations',ylab=i), #+geom_hline(yintercept = sdat[[i]], col='grey',lty='dashed')   ,
      qplot(parchain[,i],geom="density",
            xlab=i,
            # xlim=c(0,1),
            fill=I(transparent("black"))) #+geom_vline(xintercept = sdat[[i]], lty="dashed",col='grey')

    )
    plotlist[[i]]<-panel
  }

  bigpanel<-plot_grid(plotlist=plotlist,ncol=1)

  return(bigpanel)
}

parameterplot <- function(sdat,Simumode, Fitnessmode,Priormode,Proposalmode,gwatype){
    UseMethod("parameterplot")
}

parameterplot.napMCMCres <- function(sdat,
                              Simumode=1,
                               Fitnessmode=1,
                               Priormode=1,
                              Proposalmode=1
                               ){

#   gNAP<-expand.grid(simumodes,fitnessmodes,priormodes,proposalmodes)
#
#
#   for(i in 1:nrow(gNAP)){
#     Si<-fc(gNAP[i,1])
#     Fi<-fc(gNAP[i,1])
#     Pri<-fc(gNAP[i,1])
#     Pro<-fc(gNAP[i,1])
  Si<-fc(Simumode)
  Fi<-fc(Fitnessmode)
  Pri<-fc(Priormode)
  Pro<-fc(Proposalmode)

  parchain<- sdat$results[[Si]][[Fi]][[Pri]][[Pro]]$run$parchain
  parnames<- sdat$results[[Si]][[Fi]][[Pri]][[Pro]]$run$parnames

  bigpanel<-parameterplot_(parchain, parnames)

  sdat$results[[Si]][[Fi]][[Pri]][[Pro]][["parplot"]] <- bigpanel

  return(sdat)
}


#### (5) Extract information ####
extract <- function(sdat,Simumode, Fitnessmode,Priormode,Proposalmode,gwatype){
    UseMethod("extract")
}
extract.napMCMCres <- function(sdat,Simumode=1,
                               Fitnessmode=1,
                              Priormode=1,Proposalmode=1,
                              gwatype=NULL){
    Si<-fc(Simumode)
    Fi<-fc(Fitnessmode)
    Pri<-fc(Priormode)
    Pro<-fc(Proposalmode)

  if(!is.null(gwatype)){
    return( sdat$results[[Si]][[paste0("GWA",gwatype)]] )
  }else{
    return( sdat$results[[Si]][[Fi]][[Pri]][[Pro]] )
  }
}

#### (6) plot grid ####
plotgrid <- function(sdat,simumodes,fitnessmodes,priormodes,proposalmodes,gwatypes,
                     plottype, pdfname, pdf,printpdf){
    UseMethod("plotgrid")
}
plotgrid.napMCMCres <- function(sdat,

                              simumodes=c(1,2,3),
                              fitnessmodes=c(1,2,3),
                              priormodes=c(1,3),
                              proposalmodes=c(3),
                              gwatypes=c(1,3),

                              plottype="s" ,
                              pdfname="out",
                              pdf=paste0(pdfname,"_",plottype,".pdf"),
                              printpdf=TRUE
                              ){

    gNAP<-expand.grid(simumodes,fitnessmodes,priormodes,proposalmodes)
    head(gNAP)
    gGWA<-expand.grid(simumodes,gwatypes)
    head(gGWA)


    ncols= length(simumodes)
    nrows= ceiling( (nrow(gNAP) + nrow(gGWA)) / ncols )

if(plottype=='all'){
    # message("plotting all")
    grid1<-plotgrid(sdat,simumodes,fitnessmodes,priormodes,proposalmodes,gwatypes,plottype='s',spdfname,printpdf=FALSE)
    grid2<-plotgrid(sdat,simumodes,fitnessmodes,priormodes,proposalmodes,gwatypes,plottype='i',pdfname,printpdf=FALSE)
    allgrids<-plot_grid( grid1,grid2,align = c("nv"))
    if(printpdf){ save_plot(filename=pdf,plot=allgrids,base_height = 3*nrows,base_width = 3*ncols*2, limitsize=FALSE) }
    return(allgrids)
}else{

  if(plottype=="s"){
    # message("plotting s")
    plotlist= lapply(1:nrow(gNAP), function(g){ extract(sdat,gNAP[g,1],gNAP[g,2],gNAP[g,3],gNAP[g,4])$plot_s })
    plotlist= append(plotlist, lapply(1:nrow(gGWA), function(g){ extract(sdat,gGWA[g,1],gwatype = gGWA[g,2])$plot_s }) )

  }else if(plottype=="i"){
    # message("plotting i")
    plotlist= lapply(1:nrow(gNAP), function(g){ extract(sdat,gNAP[g,1],gNAP[g,2],gNAP[g,3],gNAP[g,4])$plot_i })
    plotlist= append(plotlist, lapply(1:nrow(gGWA), function(g){ extract(sdat,gGWA[g,1],gwatype = gGWA[g,2])$plot_i }) )
  }else{
    stop("Not a valid plottype, only s or i or all")
  }

    labs_a<-apply(gNAP,1,function(i) paste(i,collapse=""))
    labs_b<-paste0(apply(gGWA,1,function(i) paste(i,collapse="")),"GWA")
    sapply(1:length(plotlist),
           function(i){
              if(is.null(plotlist[[i]])){
                message("Null plot")
                message(c(labs_a,labs_b)[i] )
                plotlist[[i]]<- ggplot()
              }
           }
      )

panel<-plot_grid(plotlist = plotlist,ncol = ncols,nrow = nrows,
                # labels = c(labs_a,labs_b),
                align = c("nv"))

  if(printpdf)  {  save_plot(filename=pdf,plot=panel,base_height = 3*nrows,base_width = 3*ncols, limitsize=FALSE) }
return(  panel)
}
}


####************************************************************************####

#' Simulate
#'
#' @param Go
#' @param p
#' @param b
#' @param bmin
#' @param bmax
#' @param a
#' @param amin
#' @param amax
#' @param m
#' @param prophidden
#' @param ss
#' @param epi
#'
#' @return
#' @export
#'
#' @examples
simulate_data<-function(Go,
      p=0.1,
      b=0.1,
      a=0,
      mu=1,
      epi=1,
      svar=0.5,
      ss=0.8,
      m=10,
      trainingp=1,
      prophidden=0.0
      ){

  newtopcols<-sort(sample(1:ncol(Go),m,replace = FALSE))
  training<-sort(sample(1:nrow(Go),ceiling(trainingp*nrow(Go)), replace = FALSE))

  if(trainingp ==1){
    testing=training
  }else{
    testing=(1:nrow(Go))[-training]
  }

  s = PropoS(m,svar)
  s= s* rbinom(n = length(s),size = 1,p=1-ss)

  ## Simulate hidden selection coefficients
  if(prophidden!=0){
    mhid=floor(m*prophidden)
    hiddencols<-sample(1:ncol(Go),size = mhid,replace = FALSE)
    # s_hid = rexp(mhid,0.1);
    # s_hid=s_hid/max(s_hid,na.rm = T);
    # s_hid= s_hid* sample(x = c(-1,+1),size = mhid,replace = TRUE)
    s_hid = PropoS(mhid,svar)
  }else{
    hiddencols<-c()
    s_hid<-c()
  }

 ## Genome matrix
  X<- Go[,newtopcols]
 ## ORIENTATION!
  ori<-sample(x = c(-1,1),size = m,replace = TRUE)
  X<- t(apply(X,1,function(i){ i* ori} ) )

   ## Create expectations
  Ey1= Ey_go(X,s = c(s,s_hid),mode=1) # WITH HIDDEN COLUMNS
  Ey2= Ey_go(X,s = c(s,s_hid),mode=2) # WITH HIDDEN COLUMNS
  Ey3= Ey_go(X,s = c(s,s_hid),mode=3) # WITH HIDDEN COLUMNS

  plot_grid(
    qplot(s, main="selection coefficients"),
    qplot(Ey1,xlab="True genotype fitness 1"),
    qplot(Ey2,xlab="True genotype fitness 2"),
    qplot(Ey3,xlab="True genotype fitness 3")
  ) %>% print()


  ## Sample 5 repliates from those expectations

  Y1<-sampleEys(Eys = fn(Ey1),a = a,b = b,p = p,rep=5)
  Y2<-sampleEys(Eys = fn(Ey2),a = a,b = b,p = p,rep=5)
  Y3<-sampleEys(Eys = fn(Ey3),a = a,b = b,p = p,rep=5)


  # The ids of genotypes

  Hids<-sort(rep(1:nrow(Go), 5))


  ## Buid data
  Y<-data.frame(Y1,Y2,Y3,row=Hids)

  ## Filter data to the training ones
  # Y<-filter(Y, row %in% training)

  simu<-napMCMCres(b=b,
                    a=a,
                    p=p,
                    mu=mu,
                    epi=epi,
                    svar=svar,
                    ss=ss,
                    newtopcols=newtopcols,
                     training=training,
                     testing=testing,
                    Go=Go,
                    Y=Y,
                    Ey1=Ey1,
                    Ey2=Ey2,
                    Ey3=Ey3,
                    s=s,
                   ori=ori
                    )

  return(simu)
}


#' run
#'
#' @param sdat
#' @param Go
#' @param epi
#' @param simumodes
#' @param fitnessmodes
#' @param priormodes
#' @param proposalmodes
#' @param gwatypes
#' @param iterations
#' @param burnin
#' @param donap
#' @param dogwa
#' @param doplot
#'
#' @return
#' @export
#'
#' @examples
run_MCMC_GWA<-function(sdat,
                    Go,
                    epi=1,
                    simumodes,
                    fitnessmodes,
                    priormodes,
                    proposalmodes,
                    gwatypes,
                    iterations=5e2,
                    burnin=0.1,
                    donap=TRUE,
                    dogwa=TRUE,
                    doplot=TRUE,
                    pdfname=NULL
                    ){

gNAP<-expand.grid(simumodes,fitnessmodes,priormodes,proposalmodes)
gGWA<-expand.grid(simumodes,gwatypes)
i=1

if(donap){
    for(i in 1:nrow(gNAP)){
      message(gNAP[i,])
      sdat<-runMCMC(sdat,
                    Simumode=gNAP[i,1],
                    Fitnessmode=gNAP[i,2],
                    Priormode=gNAP[i,3],
                    Proposalmode=gNAP[i,4],
                    iterations=iterations,
                    burnin=burnin,
                    doplot=doplot
                    )

    }
}

if(dogwa){
    for(i in 1:nrow(gGWA)){
      sdat<-runGWA(sdat,
                   Simumode=gGWA[i,1],
                   gwatype=gGWA[i,2],
                   doplot=doplot)

    sdat}
}

if(!is.null(pdfname)){
  # plotgrid(sdat,simumodes,fitnessmodes,priormodes,proposalmodes,gwatypes,
  #          plottype="s",pdfname = pdfname)
  # plotgrid(sdat,simumodes,fitnessmodes,priormodes,proposalmodes,gwatypes,
  #          plottype="i",pdfname = pdfname)
plotgrid(sdat,simumodes,fitnessmodes,priormodes,proposalmodes,gwatypes,
plottype="all",pdfname = pdfname)
}
return(sdat)
}

