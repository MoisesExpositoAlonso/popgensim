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
    multievo(
                      p=as.numeric(input$p),
                      n=as.numeric(input$n),
                      rate=as.numeric(input$rate),
                      svar=as.numeric(input$svar),
                      a=as.numeric(input$a),
                      b=as.numeric(input$b),
                      mu=as.numeric(input$mu),
                      d=as.numeric(input$d),
                      No=as.numeric(input$No),
                      Nmax=as.numeric(input$Nmax),
                      tmax=as.numeric(input$tmax),
                      mode=as.numeric(input$mode),
                      Replicates=as.numeric(input$Replicates)
                    )
})

  output$Popsize <- renderPlot({
                      obj<-p()
                      plotmy<-popsizeplot(obj)
                      print(plotmy)
                      })



}

####************************************************************************####
# define ui

myui<- function(){
fluidPage(
titlePanel("Multilocus evolution of diploid clonal populations: mutation, selection, and drift"),

sidebarPanel(
    "Number of loci",
    sliderInput("p", "p:", min=1, max=10000, value=1000, step=100),
    "Number of different genomes",
    sliderInput("n", "n:", min=1, max=1000, value=100, step=10)
),sidebarPanel(
   "Site Frequency Spectrum (SFS) shape",
    sliderInput("rate", "Rate in exponential distribution:", min=1e-6, max=1e6, value=1, step=100     ),
    "Distribution of Fitness Effects (DFE) shape",
    sliderInput("svar", "Variance of the logNormal distribution:", min=0, max=1, value=0.1, step=0.1),
   "Pseudoheritability",
    sliderInput("a", "Baseline variance of Normal fitness distribution", min=0, max=1, value=0.5, step=0.1),
    sliderInput("b", "Increment in variance of Normal fitness distribution", min=0, max=1, value=1, step=0.1),
    sliderInput("mu", "Relative to absolute fitness transformation, Rel x mu = Abs", min=1e-3, max=10, value=0.8, step=0.2),

    selectInput("mode",
                  label = "Select selection mode",
                  choices = list("Multiplicative"=1,
                                 "Additive"=2,
                                 "Inverse multiplicative"=3),
                  selected = "Multiplicative")
),sidebarPanel(
  "Drift",
    sliderInput("d", "Drift coefficient: d", min=0, max=1, value=0.7, step=0.01),
    sliderInput("No", "Initial population size: N0", min=1, max=1e6, value=100000, step=10),
    sliderInput("Nmax", "Carrying capacity: Nmax", min=1, max=1e6, value=1e6, step=100)
),sidebarPanel(
    sliderInput("tmax", "Number of generations:", min=1, max=50, value=5, step=1),
    sliderInput("Replicates", "Number of replicates:", min=1, max=100, value=50, step=1)
),

  mainPanel(
    plotOutput("Popsize")
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
