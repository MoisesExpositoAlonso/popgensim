run_MCMCglmm=function(
                            dat,
                            var,
                            fixformula,
                            family,
                            randomformula=NULL,
                            strongerprior=FALSE ,
                            ...
                            )
{


# Format the formulas
message("parsing formulas")
fixformula=formula( paste(var,fixformula ))
randomformula=formula(randomformula)


# Define prior if needed
if(strongerprior == TRUE & randomformula != TRUE){
message("attemption a non-improper prior")

  randomnum=length (strsplit(as.character(randomformula),split = "+" ,fixed=TRUE)[[2]])
  Rlist= list(V = 1, nu = 0.002)
  Glist=list(V = 1, nu =0.002)
  Glists= lapply(1:randomnum, function(i) Glist )
  names(Glists) <- paste0("G",1:randomnum)
prior <- list(R = Rlist, G = Glists)

}else{ prior=NULL}

# Run the model
message("runing MCMC glmm ... ")
message(paste("  -> response variable: ",var) )
message(paste("  -> fixed effects formula: ",as.character(fixformula) ))
message(paste("  -> random effects formula: ", as.character(randomformula)) )
message(paste("  -> distribution family: ",family) )

fitmod<-MCMCglmm(data=dat,
                 family=family,
                 fixed= fixformula,
                 random= randomformula,
                 verbose=TRUE,
                 pr=TRUE,
                 pl=TRUE,
                 prior=prior, ...
                 )

return(fitmod)
}

get_effects<-function(){}
