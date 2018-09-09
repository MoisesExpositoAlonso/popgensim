# http://stat.ethz.ch/~meier/teaching/cheming/5_regression_large.pdf


G <- genomes$genotypes
X. = G[,c(1:1000)] # Instead of sampling, I get them consecutive
X. = G[,sample(1:ncol(G),200)] # Instead of sampling, I get them consecutive
X. = G[,c(1:600)] # Instead of sampling, I get them consecutive
N=nrow(X.)
M=ncol(X.)

genomes$map[1000,"physical.pos"] - genomes$map[1,"physical.pos"]


### 1.2 Mean center and variance scale
# Center and scale genome matrix
X. = apply(X., 2, function(x) { x[ x== (-1)] <- 1 ; x})
X. = apply(X., 2, function(x) { x[ is.na(x)] <- 1 ; x})
# X. = apply(X., 2, function(x) { x - 1 }) ## this centers the matrix in 0
X = apply(X.,2, function(x) {
    mu=mean(x)
    sig=sd(x)+1e-10
    (x-mu)/sig
  } )

### 1.3 Calculation of LD based on raw covariation

# Calculate
V= t(X) %*% X
V= 1/N * t(X) %*% X
D = diag(V)
qplot(x = fn(V),geom="histogram", xlab="LD distribution")

Vinv=MASS::ginv(V)

### Get the phenotypes

Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Survival_fruit_mli")], by.y="id",all.x=T)[,7]
# Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Flowering_time_mli")], by.y="id",all.x=T)[,7]
Y=relative(Y)
Y = meanvarcent(Y)
qplot(Y,geom="histogram",main="Survival to reproduction (Madrid+drought)")
Y.bad=Y
Y.bad[is.na(Y.bad)]<-mean(Y.bad,na.rm=TRUE) # Mean impute

# Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Survival_fruit_thi")], by.y="id",all.x=T)[,7]
# Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id","Flowering_time_thi")], by.y="id",all.x=T)[,7]
# Y=relative(Y)
# Y = meanvarcent(Y)
# qplot(Y,geom="histogram",main="Survival to reproduction (Tübingen+water)")
# Y.good=Y
# Y.good[is.na(Y.good)]<-mean(Y.good,na.rm=TRUE) # Mean impute

### Calculate coefficients

bgwa=sapply(1:M,function(m) solve(t(X[,m])%*%X[,m]) %*% t(X[,m]) %*% Y.bad)
# b= 1/N * Vinv %*% t(X) %*% Y.bad
lambda=0.0001
lambda=0.1
b= Vinv %*% bgwa

lambda=0
b.ridge= solve(t(X)%*%X + lambda*diag(M))  %*% t(X) %*% Y.bad

cgwa

# b.ridge= solve(V + lambda*V)  %*% t(X) %*% Y.bad
b.ridge= MASS::ginv(V + lambda*V)  %*% t(X) %*% Y.bad

b.ridge= MASS::ginv(t(X)%*%X + lambda *diag(M))  %*% t(X) %*% Y.bad
b.ridge= MASS::ginv(t(X)%*%X + seq(0.0001,0.1,length.out = M) *diag(M))  %*% t(X) %*% Y.bad

b.ridge= MASS::ginv(t(X)%*%X + normalize(abs(bgwa-b)) *diag(M))  %*% t(X) %*% Y.bad

b= MASS::ginv(t(X)%*%X)  %*% t(X) %*% Y.bad #Vinv %*% bgwa
plot(b,b.ridge)
plot(bgwa,b.ridge)
plot(bgwa,b)

qplot(y=b,x=b.ridge) + stat_smooth(method="glm")
cor.test(b.ridge,b)


MSE= 1/N * sum( (Y.bad - (X %*% b))^2 )
MSEgwa= 1/N * t(Y.bad - (X %*% bgwa)) %*% (Y.bad - (X %*% bgwa))
print(MSEgwa)
print(MSE)

b.bad=b
bgwa.bad=bgwa

# b= 1/N * Vinv %*% t(X) %*% Y.good
bgwa=sapply(1:M,function(m) solve(t(X[,m])%*%X[,m]) %*% t(X[,m]) %*% Y.good)
b= Vinv %*% bgwa

MSE= 1/N * sum( (Y.good - (X %*% b))^2 )
MSEgwa= 1/N * t(Y.good - (X %*% bgwa)) %*% (Y.good - (X %*% bgwa))
print(MSEgwa)
print(MSE)

b.good=b
bgwa.good=bgwa


library(glmnet)
# load data
data(longley)
x <- as.matrix(longley[,1:6])
y <- as.matrix(longley[,7])
# fit model
fit <- glmnet(x, y, family="gaussian", alpha=0, lambda=0.001,standardize = FALSE)
# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, X, type="link")
# summarize accuracy
mse <- mean((Y.bad - predictions)^2)
predictions

################################################################################
### See how the phenotype predicts
d=data.frame(y=Y.bad, y.cond= (X %*% b.bad), y.marg=(X %*% bgwa.bad),y.ridge=(X %*% b.ridge),y.lasso=fn(predictions))
head(d)

ggplot(data=d,
         # aes(y=y,x=y.ridge)
         aes(y=y,x=y.lasso)
         )+
  geom_point(color="red",alpha=0.2)+
         # geom_point(aes(x=Y,y=Y.cond),color="red",alpha=0.2)+
         # geom_point(aes(x=Y,y=Y.marg),color="black",alpha=0.2)+
         # coord_cartesian(ylim=c(range(Y.bad)),xlim=c(range(Y.bad)))+
         geom_abline(slope = 1,intercept = 0)+
  ylab("Real phenotype")+
  xlab("Predicted phenotype") #+
  # stat_smooth(aes(y=Y,x=Y.cond),formula = y ~ poly(x, 3),method="glm")+
  # stat_smooth(aes(y=Y,x=Y.cond),formula = y ~ poly(x, 2),method="glm")+
  # stat_smooth(aes(y=Y,x=Y.cond),formula = y ~ poly(x, 1),method="glm")+
  # stat_smooth(aes(y=Y,x=Y.cond),formula = y ~ exp(x) ,method="glm")+
  # stat_smooth(aes(y=Y,x=Y.cond),formula = y ~ exp(x) ,method="glm")


