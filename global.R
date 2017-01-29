
alleleFreq <- function(mu, nu, m, wAA, wAa, waa, p0, pc, tmax) {
    p <- numeric(tmax)
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      # mutation first
      pp <- (1-mu)*p[t] + (1-p[t])*nu

      # next, migration
      ppp <- (1-m)*pp + m*pc

      # then selection
      if ( ( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 ) > 0) {
        p[t+1] <- ( wAA*ppp^2 + wAa*ppp*(1-ppp) ) / ( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 )
      } else {
        p[t+1] <- NA
      }
    }

    return(p)
  }
