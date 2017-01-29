library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Population Genetics Models: Drift"),

  sidebarPanel(
    # navigation
    p(a("HOME", href="../index.html")),

    # Descriptive text
    p("Use the sliders below to explore the effects of changing parameters associated with drift on the evolution of a population."),

    # drift coefficient
    sliderInput("d", "Drift coefficient:", min=0, max=1, value=0.70, step=0.01),

    # starting population size
    sliderInput("N", "Initial population size:", min=1, max=1000, value=100, step=1),

    # starting A allele frequency on island
    sliderInput("p0", "Initial allele frequency of A:", min=0, max=1, value=0.50, step=0.01),

    # generation time
    sliderInput("tmax", "Number of generations:", min=1, max=1000, value=100, step=1)
  ),

  mainPanel(
    #textOutput("debug"),
    plotOutput("allelePlot"),
    plotOutput("genoPlot")
  )
))
