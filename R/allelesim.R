

allelesim <- function(
mu=0,
nu=0,
m=0,
wAA=1,
wAa=0.75,
waa=0.5,
Fi=0,
p0=0.5,
psource=0.1,
tmax=100,
d=0,
N=1000,
rep=50
  ) {

  sapply( 1:rep , FUN=function(rep){
    p <- c()
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      # mutation first
      pp <- (1-mu)*p[t] + (1-p[t])*nu

      # next, migration
      ppp <- (1-m)*pp + m*psource

      # then selection
      w_hat=( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 )
      print(w_hat)
      if ( w_hat > 0 & !is.na(w_hat) ) {
        p[t+1] <- ( wAA*ppp^2 + wAa*ppp*(1-ppp) ) / w_hat
      } else {
        p[t+1] <- NA
      }
      # then imbreeding (this equation is general)
      fAA <- (p[t+1]^2 * (1-Fi)) + p[t+1]*Fi
      fAa <- 2*p[t+1]*(1-p[t+1]) * (1-Fi)
      faa <- (1-p[t+1])^2 * (1-Fi) + (1-p[t+1])* Fi
      # no imbreeding
      # fAA <- p[t+1]^2
      # fAa <- 2*p[t+1]*(1-p[t+1])
      # faa <- (1-p[t+1])^2

      # then drift
      NAA <- round(N*fAA)
      NAa <- round(N*fAa)
      Naa <- round(N*faa)

      # if (NAA <= 0) {
      #   NAA <- 0
      #   NAA.prime <- 0
      # } else {
      #   NAA.prime <- sum(rbinom(n=NAA, size=1, prob=d))
      # }
      # if (NAa <= 0) {
      #   NAa <- 0
      #   NAa.prime <- 0
      # } else {
      #   NAa.prime <- sum(rbinom(n=NAa, size=1, prob=d))
      # }
      # if (Naa <= 0) {
      #   Naa <- 0
      #   Naa.prime <- 0
      # } else {
      #   Naa.prime <- sum(rbinom(n=Naa, size=1, prob=d))
      # }
      # N.prime <- NAA.prime + NAa.prime + Naa.prime
      #
      # if (N.prime <= 0) {
      #   p[t+1] <- NA
      # } else {
      # p[t+1] <- (2*NAA.prime + NAa.prime) / (2*N.prime)
      # }

      p[t+1] <- (2*NAA + NAa) / (2*Naa)

    } #end t loop


    return(p)
}) #end sapply
}
