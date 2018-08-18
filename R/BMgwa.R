GWA<-function(Go,
								y,
	  						vars,
	  						training=(1:nrow(Go) )-1,
	  						type=1,
					      lambda =1,
						  	max_iter=1000,
						  	tol = 1e-5
	  						){
BMgwa(Go@address,y,vars,training,type,lambda,max_iter,tol)
}