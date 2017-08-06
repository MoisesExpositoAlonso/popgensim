library(shiny)
library(ggplot2)
library(cowplot)
library(devtools)
library(Rcpp)
devtools::load_all(".")
# install(".")

# Define server
function(input, output) {


p <- reactive({
  # popgensim::allelesim(
  popgensim::allelesimCmat( # this is the C implementation
    mu=as.numeric(input$mu),
                             nu=as.numeric(input$nu),
                             m=as.numeric(input$m),
                             wAA=as.numeric(input$wAA),
                             wAa=as.numeric(input$wAa),
                             waa=as.numeric(input$waa),
                             p0=as.numeric(input$p0),
                             psource=as.numeric(input$psource),
                             tmax=as.numeric(input$tmax),
                             Fi=as.numeric(input$Fi),
                             d=as.numeric(input$d),
                             N=as.numeric(input$N),
                             rep=as.numeric(input$rep)
                             ) })

  colours <- c("#fb8072", "#80b1d3", "#7fc97f")


  output$allelePlot <- renderPlot({

  p<-ggplot() +ylim(c(0,1)) + xlim(c(0, input$tmax))+ ylab("Allele frequency") + xlab("Generations") +  scale_colour_manual("",breaks=c("A","a"),labels = c("A", "a"),values = colours[1:2])

  for(rep in 1:input$rep){
   newtoplot=data.frame(Generation=1:input$tmax,
                    A=p()[rep,],
                    a=1-p()[rep,])
   p<-p+ geom_line(data=newtoplot,aes(y=A,x=Generation,color="A")) +
     geom_line(data=newtoplot,aes(y=a,x=Generation,color="a"))
  }
  print(p)

})

  output$genoPlot <- renderPlot({

   p<-ggplot() +ylim(c(0,1)) + xlim(c(0, input$tmax))+ ylab("Genotype frequency") + xlab("Generations") +  scale_colour_manual("",breaks = c("AA", "Aa","aa"),labels = c("AA", "Aa","aa"),values = colours[c(1,3,2)])

  for(rep in 1:input$rep){

  newtoplot=data.frame(Generation=1:input$tmax,
                    AA=p()[rep,]^2,
                    aa=(1-p()[rep,])^2,
                    Aa=2*p()[rep,]*(1-p()[rep,])
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
