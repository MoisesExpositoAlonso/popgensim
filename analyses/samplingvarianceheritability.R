# heritability based on sampling variance?


w=rnorm(100,1,2)
# w=1:100
hist(w)

y=sapply(1:length(w),function(i)rnorm(1,w[i],0*abs(w[i])))

plot(w,y)
cor.test(w,y)

var(w,na.rm = T) / var(y,na.rm = T)


results<-c()

for( i in 1:100){

frac=seq(0,1,length.out = 100)[i]

y=sapply(1:length(w),function(x)rnorm(1,w[x],(1-frac)*abs(w[x])))
# y=sapply(1:length(w),function(x)rnorm(1,w[x],(1-frac)))

h<-var(w,na.rm = T) / var(y,na.rm = T)
h<-summary(lm(w~y))$coefficients[2,1]
h<-abs(cor(w,y))

results<-c(results, h)


}

plot(results~1-seq(0,1,length.out = 100))


cor(results,seq(0,1,length.out = 100))

i=0

i=1
hreal=0

moiR::randomvariance
