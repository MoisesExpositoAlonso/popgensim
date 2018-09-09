BMgwa<-function(X0, yraw,training=NULL, vars=NULL, type=1, lambda=1, max_iter=1000L, tol=1e-5){
 
  if(class(X0) == "externalptr"){
  	return(
    	BMgwa1(X0,yraw,vars,training,type,lambda,max_iter,tol)
  	)
  }else if(class(X0) == "matrix"){
    return(
  		BMgwa2(X0,yraw,type,lambda,max_iter,tol)
    )
  }else if(class(X0) == "big.matrix"){
    if(is.null(training)){training=1:nrow(X0)}
    if(is.null(vars)){vars=1:ncol(X0)}
    return(
  		BMgwa2(X0[training,vars],yraw,type,lambda,max_iter,tol)
    )
  }else{
    stop('Unknown data type')
  }

}

# napMCMC<-function (y, h, A, s, m, n ,
#                    b = 0.5, bmin = 0, bmax = 1,
#                   a = 0.4, amin = 0, amax = 1,
#                   p = 0.5,
#                   mu = 1, mumin = 0,mumax = 10,
#                   epi = 1, epimin = 1, epimax = 1,
#                   svar = .1, svarmin = 0, svarmax = 2,
#                   ss = .1, ssmin = 0, ssmax = 1,
# 
#                   bw =  ifelse(Fitnessmode ==3L, 1, 0.1),
#                   nupdates = 1L,
# 
#                   min = ifelse(Fitnessmode ==3L, 0, -1),
#                   max = ifelse(Fitnessmode ==3L, +2.5, 1),
# 
#                   iterations = 10000,
#                   TEST = FALSE,
#                   verbose = FALSE,
#                   debug = FALSE,
#                   Fitnessmode = 1L,
#                   Priormode = 1L,
#                   Proposalmode = ifelse(Fitnessmode ==3L, 2L, 1L),
#                   file2sink = "output.log", sink2file = FALSE){
# 
# 
#   BMgwa2(y, h, A, s, m, n,
#     b, bmin, bmax, a, amin, amax, p, mu, mumin, mumax, epi, epimin,
#     epimax, svar,svarmin,svarmax,ss,ssmin,ssmax, bw, nupdates, min, max, iterations, TEST, verbose,
#     debug, fn(Fitnessmode), fn(Priormode), fn(Proposalmode), file2sink, sink2file)
# 
# 
# }



# gwa<-function(X,Y,M=ncol(X),N=nrow(X)){
#   X=meanvarcent.mat(X)
#   bgwa=sapply(1:M,function(m) solve(t(X[,m])%*%X[,m]) %*% t(X[,m]) %*% Y)
#   MSEgwa= 1/N * t(Y - (X %*% bgwa)) %*% (Y - (X %*% bgwa))
#   return(bgwa)
# }
#
# logistic<-function(x,L=1,x0=0.05,k=10){
#   L / (1 + exp(-k*(x-x0)))
# }
#
# ridgegwa<-function(X,Y,M=ncol(X),N=nrow(X), lambda=1, type='regular',k=100,x0=0.90){
#   # Maybe implement an empirical approach
#
#   if(type=='regular'){
#     b.ridge= MASS::ginv(t(X)%*%X + lambda *diag(M))  %*% t(X) %*% Y
#   }else if(type=='ld'){
#     b=cgwa(X,Y)
#     bgwa=gwa(X,Y)
#     b.ridge= MASS::ginv(t(X)%*%X + lambda*(normalize(abs(bgwa-b))) *diag(M))  %*% t(X) %*% Y
#   }else if(type=='ldpenalized'){
#     b=cgwa(X,Y)
#     bgwa=gwa(X,Y)
#     b.ridge= MASS::ginv(t(X)%*%X + lambda*(1-normalize(logistic(normalize(abs(bgwa-b)),k=k,x0 =x0))) *diag(M))  %*% t(X) %*% Y
#   }else if(type=='penalized'){
#     b=cgwa(X,Y)
#     # b.ridge= MASS::ginv(t(X)%*%X + lambda*(1-normalize(logistic(normalize(abs(b)),k=k,x0 = x0))) *diag(M))  %*% t(X) %*% Y
#     b.ridge= MASS::ginv(t(X)%*%X + lambda*(normalize(logistic(-normalize(abs(b)),k=k,x0 = x0))) *diag(M))  %*% t(X) %*% Y
#   }else if(type=='threshold'){
#     b=cgwa(X,Y)
#     b.ridge= MASS::ginv(t(X)%*%X + lambda*(1-normalize(logistic(normalize(abs(b)),k=10e6,x0 =x0))) *diag(M))  %*% t(X) %*% Y
#   }else{
#     stop('Did not recognize the type. Choose among: regular, ld, penalized, threshold')
#   }
#   return(b.ridge)
# }
#
#
# cgwa<-function(X,Y,Vinv=NULL,M=ncol(X),N=nrow(X)){
#
#   X=meanvarcent.mat(X)
#
#   # b= 1/N * Vinv %*% t(X) %*% Y
#   b= MASS::ginv(t(X)%*%X)  %*% t(X) %*% Y
#
#   # MSE= 1/N * sum( (Y - (X %*% b))^2 )
#   # message("The residual variance is: ",MSE)
#   return(b)
# }
#
# lassogwa<-function(X,Y,M=ncol(X),N=nrow(X)){
#   suppressMessages(require(glmnet))
#
#   X=meanvarcent.mat(X)
#
#   cvfit = cv.glmnet(X, Y)
#   # plot(cvfit)
#   # cvfit$lambda.min
#
#
#   # fit <- glmnet(X,Y, family="gaussian", alpha=0, lambda=0.001)
#   fit <- glmnet(X,Y, family="gaussian",lambda = cvfit$lambda.min)
#   fit
#
#   # predictions <- predict(fit, X, type="link")
#   #
#   # mse <- mean((Y - predictions)^2)
#   # message("The mean square error is ", mse)
#
#   return(fn(fit$beta))
# }
#
#
# getVinv<-function(X,M=ncol(X),N=nrow(X)){
#   V= 1/N * t(X) %*% X
#   Vinv=MASS::ginv(V)
#   return(Vinv)
# }
