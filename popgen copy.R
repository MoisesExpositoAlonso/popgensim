  alleleFreq <- function(d=0.5, N=100, p0=0.5, tmax=50) {
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


dev.off()

plot.alleleFreq <- function(rep, ...){

sapply(1:rep,  FUN = function(i){
        if(i==1){ plot(alleleFreq(...),type = "l",axes = F, xlab = "",ylab = "")  }
        else{  points(alleleFreq(...),type = "l",axes = F,xlab = "",ylab = "") }
})
}

plot.alleleFreq(3,d=0.9,tmax=1000,N=10000)
