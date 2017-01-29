library(shiny)

# Define UI for miles per gallon application
fluidPage(

  # Application title
  headerPanel("Evolution of diploid populations: popgen equations for mutation, migration, selection, and drift"),

  sidebarPanel(
    # navigation
    # p(a("HOME", href="../index.html")),

    # Descriptive text
    p("Use the sliders below to explore the effects of changing parameters associated with a combination of evolutionary forces on the evolution of a population."),


    # mutation rate A->a
    sliderInput("mu", "Mutation rate, mu:", min=0, max=0.001, value=0.0001, step=1*10^-5),

    # mutation rate a->A
    sliderInput("nu", "Mutation rate, nu:", min=0, max=0.001, value=0.00001, step=1*10^-5),

    # migration rate
    sliderInput("m", "Migration rate, m:", min=0, max=1.0, value=0.05, step=0.001),

    # starting A allele frequency on continent (mainland)
    sliderInput("pc", "Allele frequency on continent, pC:", min=0, max=1, value=0.90, step=0.01),

    # relatvie fitnesses
    sliderInput("wAA", "Relative fitness of AA:", min=0, max=1.0, value=1.00, step=0.01),
    sliderInput("wAa", "Relative fitness of Aa:", min=0, max=1.0, value=1.00, step=0.01),
    sliderInput("waa", "Relative fitness of aa:", min=0, max=1.0, value=0.90, step=0.01),

    # starting A allele frequency on island
    sliderInput("p0", "Initial allele frequency of A on thi island, pI:", min=0, max=1, value=0.10, step=0.01),

    # generation time
    sliderInput("tmax", "Number of generations:", min=1, max=1000, value=100, step=1)
  ),

  mainPanel(
    # textOutput("debug"),
    plotOutput("allelePlot"),
    plotOutput("genoPlot")
  )
)
