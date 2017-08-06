# Interactive shiny R plot of allele frequencies under a Wright-Fisher population

Tune the forces of evolution: mutation, migration, selection, and drift (including also non-random mating); to simulate allele frequencies under the classic population genetics equations.


# Get it

## Option 1 directly from github

``` sh
# Packages required
library(ggplot)
library(cowplot
library(shiny) 
library(Rcpp)

# To start the app just run:
runGitHub("popgensim","MoisesExpositoAlonso")

``` 


### Option 2 Download package
``` sh
git clone https://github.com/MoisesExpositoAlonso/popgensim.git

``` 


# Run simulations independently of the Shinny app
Note: You need to download the package before

```
# Packages required
library(devtools)
library(ggplot) 
library(cowplot)
library(shiny) 
library(Rcpp)

install(".") # install the popgensim package


# Run a simulations defining each population parameter manually:

allelesimCmat(
    mu=0.001,
    nu=0.001,
    m=0,
    wAA=0.5,
    wAa=0.5,
    waa=0.5,
    p0=0.5,
    psource=0.5,
    N=1000,
    tmax=100,
    rep=50)

# To compare the speed of the R and C++ implementation:

library(microbenchmark) 
microbenchmark(
  Rimplement=allelesim(rep = 10,tmax=50),
    Cimplement=allelesimCmat(rep = 10,tmax=50)
    )



```
