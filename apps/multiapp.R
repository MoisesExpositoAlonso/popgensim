library(shiny)
library(ggplot2)
library(cowplot)
library(devtools)
library(Rcpp)
devtools::load_all(".")
# install(".")

####************************************************************************####
# Define server
myserver<-function(input, output) {

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
}

####************************************************************************####
# define ui

myui<- function(){
fluidPage(
titlePanel("One biallelic locus evolution of diploid populations: mutation, migration, selection, and drift"),
  # Application title
  # headerPanel("One locus evolution of diploid populations: mutation, migration, selection, and drift"),

sidebarPanel(
    "The frequencies of A (p) and a (q=1-p) depend on the mutation rate in both directions: p = (1-mu)*p + q*nu ",
    # mutation rate A->a
    sliderInput("mu", "Mutation rate A->a, mu:", min=0, max=0.001, value=0.0001, step=1*10^-5),

    # mutation rate a->A
    sliderInput("nu", "Retro-mutation rate a->A, nu:", min=0, max=0.001, value=0.0001, step=1*10^-5)#,

),sidebarPanel(
   "The frequencies of A (p) can increase as a funciton of the proportion of new migrants comprising the population (m) and the frequency of A in the source population (psource):\n (1-m)*p + m*psource \n",
    # migration rate
    sliderInput("m", "Migration rate of A from source population, m:", min=0, max=1.0, value=0, step=0.001),

    # starting A allele frequency on continent (mainland)
    sliderInput("psource", "Allele frequency of A in source population, psource", min=0, max=1, value=0.90, step=0.01) #,

),sidebarPanel(
    " Allele frequency of A based on relative fitness of the three genotypes. The equation is basically the relative fitness of all individuals carrying at least one A, divided by the total population fitness:
    ( wAA*p^2 + wAa*p*(1-p) ) / ( wAA*p^2 + wAa*2*p*(1-p) + waa*(1-p)^2 )",
    # relatvie fitnesses
    sliderInput("wAA", "Relative fitness of AA: wAA", min=0, max=1.0, value=1.00, step=0.01),
    sliderInput("wAa", "Relative fitness of Aa: wAa", min=0, max=1.0, value=0.9, step=0.01),
    sliderInput("waa", "Relative fitness of aa: waa", min=0, max=1.0, value=0.85, step=0.01) #,

),sidebarPanel(
  "Drift is expressed as a proportion of individuals of all in the populations that survive the 'extinction by chance': N1 = d * N0",
    # drift coefficient
    sliderInput("d", "Drift coefficient: d", min=0, max=1, value=0.50, step=0.01),

    # starting population size
    sliderInput("N", "Initial population size: N0", min=1, max=1000, value=100, step=1) #,

),sidebarPanel(
  "Inbreeding is expressed as the proportion of genotypes that reproduce with themselves: example fAA = p^2 * (1-Fi) + p*Fi",
    # inbreeding coefficient
    sliderInput("Fi", "Inbreeding coefficient: Fi", min=0, max=1, value=0, step=0.01)

),sidebarPanel(
    # starting A allele frequency on focal population
    sliderInput("p0", "Initial allele frequency of A, p0:", min=0, max=1, value=0.5, step=0.01),
    # generation time
    sliderInput("tmax", "Number of generations:", min=1, max=1000, value=50, step=1),

    # replicates
    sliderInput("rep", "Number of independent loci evolving (replicates):", min=1, max=50, value=10, step=1)

  ),

  mainPanel(
    plotOutput("allelePlot") ,
    plotOutput("genoPlot")
  )
)
}
####************************************************************************####
# run app

app <- shinyApp(
    ui = myui,
    server = myserver
  )
runApp(app)
