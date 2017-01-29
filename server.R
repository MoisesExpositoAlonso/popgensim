library(shiny)
library(ggplot2)
library(cowplot)

# source("global.R")

# Define server
function(input, output) {

alleleFreq <- function(mu, nu, m, wAA, wAa, waa, p0, pc, tmax, d, N,rep) {
    p <- numeric(tmax)
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      # mutation first
      pp <- (1-mu)*p[t] + (1-p[t])*nu

      # next, migration
      ppp <- (1-m)*pp + m*pc

      # then selection
      if ( ( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 ) > 0) {
        p[t+1] <- ( wAA*ppp^2 + wAa*ppp*(1-ppp) ) / ( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 )
      } else {
        p[t+1] <- NA
      }

      # then drift
      NAA <- round(N*p[t+1]^2)
      NAa <- round(N*2*p[t+1]*(1-p[t+1]))
      Naa <- round(N*(1-p[t+1])^2)

      if (NAA <= 0) {
        NAA <- 0
        NAA.prime <- 0
      } else {
        NAA.prime <- sum(rbinom(n=NAA, size=1, prob=d))
      }
      if (NAa <= 0) {
        NAa <- 0
        NAa.prime <- 0
      } else {
        NAa.prime <- sum(rbinom(n=NAa, size=1, prob=d))
      }
      if (Naa <= 0) {
        Naa <- 0
        Naa.prime <- 0
      } else {
        Naa.prime <- sum(rbinom(n=Naa, size=1, prob=d))
      }
      N.prime <- NAA.prime + NAa.prime + Naa.prime

      if (N.prime <= 0) {
        p[t+1] <- NA
      } else {
        p[t+1] <- (2*NAA.prime + NAa.prime) / (2*N.prime)
      }
    } #end t loop


    return(p)
  }

  p <- reactive({ alleleFreq(mu=as.numeric(input$mu),
                             nu=as.numeric(input$nu),
                             m=as.numeric(input$m),
                             wAA=as.numeric(input$wAA),
                             wAa=as.numeric(input$wAa),
                             waa=as.numeric(input$waa),
                             p0=as.numeric(input$p0),
                             pc=as.numeric(input$pc),
                             tmax=as.numeric(input$tmax),
                             d=as.numeric(input$d),
                             N=as.numeric(input$N),
                             rep=as.numeric(input$rep)
                             ) })

  colours <- c("blue", "purple", "green")


  output$allelePlot <- renderPlot({
  toplot=data.frame(Generation=1:input$tmax,
                    A=p(),
                    a=1-p())
  p<-ggplot(data=toplot) + geom_line(aes(y=A,x=Generation),color=colours[1]) + geom_line(aes(y=a,x=Generation),color=colours[2]) + ylab("Allele frequency") + xlab("Generations")

  if(input$rep >1){
  for(rep in 1:input$rep){
   newtoplot=data.frame(Generation=1:input$tmax,
                    A=p(),
                    a=1-p())
   p<-p+ geom_line(data=newtoplot,aes(y=A,x=Generation),color=colours[1])
  }
  }
  print(p)
  })

  output$genoPlot <- renderPlot({
  # toplot2=data.frame(Generation=1:input$tmax,
  #                   AA=p()^2,
  #                   Aa=2*p()*(1-p()),
  #                   aa=(1-p())^2
  #                   )
  # ggplot(data=toplot2) + geom_line(aes(y=AA,x=Generation),color=colours[1]) + geom_line(aes(y=aa,x=Generation),color=colours[2]) + geom_line(aes(y=Aa,x=Generation),color=colours[3]) + ylab("Genotype frequency") + xlab("Generations")
  #
  })

  # debug only
  output$debug <- renderText({ p() })
}
