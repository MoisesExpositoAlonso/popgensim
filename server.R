library(shiny)
library(ggplot2)
library(cowplot)


# Define server
function(input, output) {


alleleFreq <- function(mu, nu, m, wAA, wAa, waa, p0, psource, tmax, d, N,rep) {
sapply( 1:rep , FUN=function(rep){
    p <- c()
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      # mutation first
      pp <- (1-mu)*p[t] + (1-p[t])*nu

      # next, migration
      ppp <- (1-m)*pp + m*psource

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
}) #end sapply
}


p <- reactive({ alleleFreq(mu=as.numeric(input$mu),
                             nu=as.numeric(input$nu),
                             m=as.numeric(input$m),
                             wAA=as.numeric(input$wAA),
                             wAa=as.numeric(input$wAa),
                             waa=as.numeric(input$waa),
                             p0=as.numeric(input$p0),
                             psource=as.numeric(input$psource),
                             tmax=as.numeric(input$tmax),
                             d=as.numeric(input$d),
                             N=as.numeric(input$N),
                             rep=as.numeric(input$rep)
                             ) })

  colours <- c("red", "blue", "green")


  output$allelePlot <- renderPlot({
  p<-ggplot() +ylim(c(0,1)) + xlim(c(0, input$tmax))+ ylab("Allele frequency") + xlab("Generations") +  scale_colour_manual("",labels = c("A", "a"),values = colours[1:2])

  for(rep in 1:input$rep){
   newtoplot=data.frame(Generation=1:input$tmax,
                    A=p()[,rep],
                    a=1-p()[,rep])
   p<-p+ geom_line(data=newtoplot,aes(y=A,x=Generation,color="A")) +
     geom_line(data=newtoplot,aes(y=a,x=Generation,color="a"))
  }
  print(p)
})

  output$genoPlot <- renderPlot({
  p<-ggplot() +ylim(c(0,1)) + xlim(c(0, input$tmax))+ ylab("Genotype frequency") + xlab("Generations") +  scale_colour_manual("",breaks = c("AA", "aa","Aa"),labels = c("AA", "aa","Aa"),values = colours)

  for(rep in 1:input$rep){

  newtoplot=data.frame(Generation=1:input$tmax,
                    AA=p()[,rep]^2,
                    aa=(1-p()[,rep])^2,
                    Aa=2*p()[,rep]*(1-p()[,rep])
                    )
  p=p+ geom_line(data=newtoplot,aes(y=AA,x=Generation,color="AA")) +
    geom_line(data=newtoplot,aes(y=aa,x=Generation,color="aa")) +
    geom_line(data=newtoplot,aes(y=Aa,x=Generation,color="Aa"))

  }
print(p)
})

  # output$histogram <- renderPlot({
  # newtoplot=data.frame(Generation=1:input$tmax,
  #                 AA=p()[input$tmax,]^2,
  #                 aa=(1-p()[input$tmax,])^2,
  #                 Aa=2*p()[input$tmax,]*(1-p()[,rep])
  #                 )
  # p <- ggplot(data=newtoplot)+geom_density
  # print(p)
  #
  # })

}
