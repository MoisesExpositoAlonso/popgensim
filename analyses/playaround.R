library(ggplot2)
library(cowplot)

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


basefreqplot <- function(tmax){
    p<-ggplot() +ylim(c(0,1)) + xlim(c(0, tmax))+ ylab("Allele frequency") + xlab("Generations")
return(p)
}

addtrajectories <-function(ps,rep,tmax, thecol="grey"){

  for(nrep in 1:rep){
   newtoplot=data.frame(Generation=1:tmax,
                    A=ps[,nrep],
                    a=1-ps[,nrep])
   p<-p+ geom_line(data=newtoplot,aes(y=A,x=Generation),color=thecol) #+  geom_line(data=newtoplot,aes(y=a,x=Generation,color="a"))

}
  return(p)
}
# ---------------------------------------------------------------------



mu=0.001
nu=0.001
m=0
wAA=0.5
wAa=0.5
waa=0.5
p0=0.5
psource=0.5
tmax=100
d=0.3
N=1000
rep=50

pn<-alleleFreq(mu, nu, m, wAA, wAa, waa, p0, psource, tmax, d, N,rep)

ps<-alleleFreq(mu, nu, m, wAA=0.55, wAa, waa, p0, psource, tmax, d, N,rep=10)

library(moiR)
library(ggplot2); library(cowplot)
p <- basefreqplot(tmax)
p<-addtrajectories(pn,rep,tmax,thecol=transparent("grey") )
p<-addtrajectories(ps,rep=10,tmax,thecol="green4")
p
