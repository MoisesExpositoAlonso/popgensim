alleleFreq <- function(mu, nu, m, wAA, wAa, waa, p0, psource, tmax, d, Fi, N,rep) {
    p <- c()
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      # mutation first
      pp <- (1-mu)*p[t] + (1-p[t])*nu

      # next, migration
      if(m!=0) ppp <- (1-m)*pp + m*psource
      else ppp=pp

      # then selection
      meanfit=( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 )

      if (is.na(meanfit)){  p[t+1] <- NA
      }else if (meanfit  > 0) {
        p[t+1] <- ( wAA*ppp^2 + wAa*ppp*(1-ppp) ) / ( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 )
      }
      else {
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

      if(d!=0){
      # then drift
      NAA <- round(N*fAA)
      NAa <- round(N*fAa)
      Naa <- round(N*faa)

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
    } #end t loop


    return(p)
}) #end sapply
}

tmax=50
mu= 7e-9
nu=1e-20
m=0; psource=0
Fi=0.99
N=1000
rep=1
d=0 #; d=0.5; d=0.99
d=0.9

sel=0.1
wAA=1; waa= wAA-sel; wAa= (wAA+waa)/2


p0=seq(0,1,0.01) [2]


sim<- alleleFreq(mu, nu, m, wAA, wAa, waa, p0, psource, tmax, d, Fi, N,rep)

#### one selection and frequency array ####
present_freq[,"tuebingen"]

farray <- sapply(present_freq[,"tuebingen"],FUN = function(p0){
  alleleFreq(mu, nu, m, wAA, wAa, waa, p0, psource, tmax, d, Fi, N,rep)
})
meltsimarray<- reshape2::melt(farray)
boxplot(meltsimarray$value ~meltsimarray$Var1, color="black")

ggplot(data=meltsimarray)+geom_boxplot(aes(y=value,x=Var1, group=Var1),outlier.size = 0,outlier.colour = "white",fill="black")

ggplot(data=meltsimarray)+geom_violin(aes(y=value,x=Var1, group=Var1),fill="black", trim = T,draw_quantiles = T)

#### selection array and frequency array ####
sel=seq(.001,.2,length.out = 100)
s2array <-sapply(sel, FUN=function(s){
  sel=s
  wAA=1; waa= wAA-sel; wAa= (wAA+waa)/2

  farray <- sapply(present_freq[,"tuebingen"],FUN = function(p0){
    alleleFreq(mu, nu, m, wAA, wAa, waa, p0, psource, tmax, d, Fi, N,rep)
  })[tmax,]

})
length(s2array)
dim(s2array)
class(s2array)
str(s2array)

qplot

ms2array<- reshape2::melt(s2array)
ms2array[,2] <- sel[ms2array[,2]]

qplot(x=ms2array[,2],y=ms2array[,3],alpha=0.1)
qplot(x=ms2array[,2],y=ms2array[,3],geom="bin2d")

table(ms2array[,3] >.9 , round(ms2array[,2],digits=3 ) )

# possel<- sort(rep(1:50,151) )
# length(possel)
# ms2array$sel=possel
# head(ms2array)
#
# meltsimarray
#
# # selection - allele freq distribution
#
#
# ggplot()+
#   # geom_boxplot(data=meltsimarray,aes(y=value,x=Var1, group=Var1),outlier.size = 0,outlier.colour = "white",fill="black")+
#   geom_boxplot(data=dplyr::filter(ms2array,sel==1),aes(y=value,x=Var2, group=Var),outlier.size = 0,outlier.colour = "white",fill="black")
#
# # meltsimarray<- reshape2::melt(s2array[[1]])
# # ggplot(data=meltsimarray)+geom_boxplot(aes(y=value,x=Var1, group=Var1),outlier.size = 0,outlier.colour = "white",fill="black")
# #
# # meltsimarray<- reshape2::melt(s2array[[50]])
# # ggplot(data=meltsimarray)+geom_boxplot(aes(y=value,x=Var1, group=Var1),outlier.size = 0,outlier.colour = "white",fill="black")
# #
# # meltsimarray<- reshape2::melt(s2array[[100]])
# # ggplot(data=meltsimarray)+geom_boxplot(aes(y=value,x=Var1, group=Var1),outlier.size = 0,outlier.colour = "white",fill="black")
#
#
#
# # library(plotly)
# #
# # library(plotly)
# # mydata = read.csv("density_plot.txt")
# # df = as.data.frame(mydata)
# # plot_ly(df, x = Y, y = X, z = Z, group = X, type = "scatter3d", mode = "lines")
# #
# #
# # plot_ly(df, x = Y, y = X, z = Z, group = X, type = "scatter3d", mode = "lines")
# #
# #
# # ggplot(all.complete, aes(x=humid_temp)) +
# #   geom_density(aes(group=height, colour=height, fill=height.f, alpha=0.1)) +
# #   guides(fill = guide_legend(override.aes = list(colour = NULL))) +
# #   labs(main="Temperature by Height", x="Temperature", y="Density")
