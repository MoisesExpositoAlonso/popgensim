library(shiny)
library(ggplot2)
library(cowplot)

# source("global.R")

# Define server
function(input, output) {

alleleFreq <- function(mu, nu, m, wAA, wAa, waa, p0, pc, tmax) {
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
    }
      # then drift

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
                             tmax=as.numeric(input$tmax)) })

  colours <- c("darkgreen", "purple", "blue")

  output$allelePlot <- renderPlot({
    # par(mar=c(5,5,1,1))
    # plot(1:input$tmax, p(), type="l", xlim=c(0,input$tmax), ylim=c(0,1), xlab="Generation", ylab="Allele frequency", col=colours[1], lty=1, lwd=3, cex.lab=2)
    # points(1:input$tmax, 1-p(), type="l", col=colours[3], lty=3, lwd=3)
    # abline(h=input$pc, lty=3, lwd=3)
    # legend("topright", legend=c("A", "a"), col=colours[c(1,3)], lty=c(1,3), lwd=3, cex=2)
  toplot=data.frame(Generation=1:input$tmax,
                    A=p(),
                    a=1-p())
  ggplot(data=toplot) + geom_line(aes(y=A,x=Generation))
  })

  output$genoPlot <- renderPlot({
    par(mar=c(5,5,1,1))
    plot(1:input$tmax, p()^2, type="l", xlim=c(0,input$tmax), ylim=c(0,1), xlab="Generation", ylab="Genotype frequency", col=colours[1], lty=1, lwd=3, cex.lab=2)
    points(1:input$tmax, (1-p())^2, type="l", col=colours[3], lty=3, lwd=3)
    points(1:input$tmax, 2*p()*(1-p()), type="l", col=colours[2], lty=2, lwd=3)
    legend("topright", legend=c("AA", "Aa", "aa"), col=colours, lty=c(1,2,3), lwd=3, cex=2)
  })

  # debug only
  output$debug <- renderText({ p() })
}
