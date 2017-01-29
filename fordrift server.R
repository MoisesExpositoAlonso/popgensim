library(shiny)

# Define server
shinyServer(function(input, output) {

  alleleFreq <- function(d, N, p0, tmax) {
    p <- numeric(tmax)
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      NAA <- round(N*p[t]^2)
      NAa <- round(N*2*p[t]*(1-p[t]))
      Naa <- round(N*(1-p[t])^2)

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
    }

    return(p)
  }

  p1 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p2 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p3 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p4 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p5 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p6 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p7 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p8 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p9 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })
  p10 <- reactive({ alleleFreq(d=as.numeric(input$d), N=as.numeric(input$N), p0=as.numeric(input$p0), tmax=as.numeric(input$tmax)) })

  colours <- c("darkred", "purple", "blue")

  output$allelePlot <- renderPlot({
    par(mar=c(5,5,1,1))
    plot(1:input$tmax, p1(), type="l", xlim=c(0,input$tmax), ylim=c(0,1), xlab="Generation", ylab="Allele frequency", col=colours[1], lty=1, lwd=3, cex.lab=2)
    points(1:input$tmax, p2(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p3(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p4(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p5(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p6(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p7(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p8(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p9(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p10(), type="l", col="lightgrey", lty=1, lwd=3)
    points(1:input$tmax, p1(), type="l", col=colours[1], lty=1, lwd=3)
    abline(h=input$pc, lty=3, lwd=3)
    legend("topright", legend=c("A"), lty=1, lwd=3, cex=2)
  })

  output$genoPlot <- renderPlot({
    par(mar=c(5,5,1,1))
    plot(1:input$tmax, p1()^2, type="l", xlim=c(0,input$tmax), ylim=c(0,1), xlab="Generation", ylab="Genotype frequency", col=colours[1], lty=1, lwd=3, , cex.lab=2)
    points(1:input$tmax, (1-p1())^2, type="l", col=colours[3], lty=3, lwd=3)
    points(1:input$tmax, 2*p1()*(1-p1()), type="l", col=colours[2], lty=2, lwd=3)
    legend("topright", legend=c("AA", "Aa", "aa"), col=colours, lty=c(1,2,3), lwd=3, cex=2)
  })

  # debug only
  output$debug <- renderText({ p() })
})
