library(shiny)

# Define UI for miles per gallon application
fluidPage(
titlePanel("One biallelic locus evolution of diploid populations: mutation, migration, selection, and drift"),
  # Application title
  # headerPanel("One locus evolution of diploid populations: mutation, migration, selection, and drift"),

sidebarPanel(
    "The frequencies of A (p) and a (q=1-p) depen on the mutation rate in both directions: \n p = (1-mu)*p + q*nu ",
    # mutation rate A->a
    sliderInput("mu", "Mutation rate A->a, mu:", min=0, max=0.001, value=0.0001, step=1*10^-5),

    # mutation rate a->A
    sliderInput("nu", "Retro-mutation ratea->A, nu:", min=0, max=0.001, value=0.00001, step=1*10^-5)#,

),sidebarPanel(
   "The frequencies of A (p) can increase as a funciton of the proportion of new migrants comprising the population (m) and the frequency of A in the source population (psource):\n (1-m)*p + m*psource \n",
    # migration rate
    sliderInput("m", "Migration rate of A from source population, m:", min=0, max=1.0, value=0.05, step=0.001),

    # starting A allele frequency on continent (mainland)
    sliderInput("pc", "Allele frequency of A in source population, psource", min=0, max=1, value=0.90, step=0.01) #,

),sidebarPanel(
    " Allele frequency of A based on fitness of the three genotypes:
    ( wAA*p^2 + wAa*p*(1-p) ) / ( wAA*p^2 + wAa*2*p*(1-p) + waa*(1-p)^2 )",
    # relatvie fitnesses
    sliderInput("wAA", "Relative fitness of AA:", min=0, max=1.0, value=1.00, step=0.01),
    sliderInput("wAa", "Relative fitness of Aa:", min=0, max=1.0, value=1.00, step=0.01),
    sliderInput("waa", "Relative fitness of aa:", min=0, max=1.0, value=0.90, step=0.01) #,

),sidebarPanel(
    # drift coefficient
    sliderInput("d", "Proportion of individuals that randomly survive to next generation, drift coefficient :", min=0, max=1, value=0.70, step=0.01),

    # starting population size
    sliderInput("N", "Initial population size:", min=1, max=1000, value=100, step=1) #,

),sidebarPanel(
    # starting A allele frequency on focal population
    sliderInput("p0", "Initial allele frequency of A on local population, p0:", min=0, max=1, value=0.10, step=0.01),
    # generation time
    sliderInput("tmax", "Number of generations:", min=1, max=1000, value=100, step=1),

    # replicates
    sliderInput("rep", "Number of evolution replicates:", min=1, max=50, value=1, step=1)

  ),

  mainPanel(
    # textOutput("debug"),
    plotOutput("allelePlot") ,
    plotOutput("genoPlot")
  )
)
