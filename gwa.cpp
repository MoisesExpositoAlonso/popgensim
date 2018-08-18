#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <string>

#define ARMA_64BIT_WORD 1
//// https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define

// when armadillo is loaded, remove this below
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>



////////////////////////////////////////////////////////////////////////////////

// // [[Rcpp::export]]
// arma::vec gwa(arma::mat X, arma::vec Y,M=ncol(X),N=nrow(X)){
//   X=meanvarcent.mat(X)
//   bgwa=sapply(1:M,function(m) solve(t(X[,m])%*%X[,m]) %*% t(X[,m]) %*% Y)
//   MSEgwa= 1/N * t(Y - (X %*% bgwa)) %*% (Y - (X %*% bgwa))
//   return(bgwa)
// }



// [[Rcpp::export]]
arma::vec polypred(arma::vec y, arma::vec x, arma::vec xnew,int order=2){
	arma::vec p = arma::polyfit(x,y,10);
	arma::vec ynew = arma::polyval(p,x);
	return(ynew);
}




// arma::Mat<int> BMsubset(SEXP & A, const arma::uvec & mycols){
//       Rcpp::XPtr<BigMatrix> bigMat(A);
//       arma::Mat<int> X((int*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
//                                       // consider saying true, perhaps is faster
//       return(X.cols(mycols));
// }
arma::Mat<double> BMread(SEXP A){
      Rcpp::XPtr<BigMatrix> bigMat(A);
			arma::Mat<double> X;

			if( bigMat->matrix_type() == 4){
 	      arma::Mat<int> X0((int*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
	    	X = arma::conv_to< arma::Mat<double> >::from(X0);
				return(X);
			}else{
 	      arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
				X=X0;
			  return(X0);
			}
}

// [[Rcpp::export]]
arma::Mat<double> BMsubset(SEXP A, const arma::uvec & myrows, const arma::uvec & mycols ){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
                                      // consider saying true, perhaps is faster
      // Subset matrix
    	if(myrows.n_elem == X0.n_rows){
    		X0=X0.cols(mycols);
    	}else if(mycols.n_elem == X0.n_rows){
    	  X0=X0.rows(myrows);
    	}else{
    		X0=X0.submat(myrows,mycols);
    	}
      return(X0);
}
// [[Rcpp::export]]
arma::Mat<double> BMcolsubset(SEXP A, const arma::uvec & mycols ){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
      // Subset matrix
    		X0=X0.cols(mycols);
      return(X0);
}

// [[Rcpp::export]]
NumericVector BMprod(XPtr<BigMatrix> bMPtr, const NumericVector& x,
                     const arma::uvec & myrows,const arma::uvec & mycols ) {

  MatrixAccessor<double> macc(*bMPtr);

  NumericVector res(myrows.n_elem);
  int i, j;
  for (j = 0; j <mycols.n_elem; j++) {
    for (i = 0; i < myrows.n_elem; i++) {
      res[i] += macc[mycols(j)][myrows(i)] * x[j];
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericVector BMpred(XPtr<BigMatrix> bMPtr, const NumericVector& x,
                     const arma::uvec & myrows,const arma::uvec & mycols,
                     double intercept){
  return(BMprod(bMPtr, x,myrows,mycols) + intercept);
}


// [[Rcpp::export]]
mat Xmcenter(mat X){
  mat newX(X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols; j++){
   newX.col(j) = (X.col(j) - arma::mean( X.col(j))) ;
  }
  return(newX);
}

// [[Rcpp::export]]
arma::colvec My(const arma::colvec & y, const arma::colvec & h){
  /*
  * Mean trait per genotype
  */

  // Declarations
  arma::colvec hunique = unique(h);
  arma::colvec m(hunique.n_elem);
  // cout << hunique << endl;


  for(int i=0; i< hunique.n_elem; i++ ){
    // Create temporal vector
    arma::colvec ytmp;
    // Fill with all values corresponding to the same genotype
    for(int j=0; j<y.n_elem;j++){
      if(h(j) == hunique(i)) {
        ytmp.resize(ytmp.size()+1);
        ytmp(ytmp.size()-1) = y(j);
      }
    }
    // Compute variance
      if(ytmp.n_elem ==1){
       // v(i)=0;
       	m(i)=ytmp(0);
      }else{
      	m(i)=arma::mean( ytmp );
      }
  }
  return(m);
}


////////////////////////////////////////////////////////////////
///  GWA
////////////////////////////////////////////////////////////////



// [[Rcpp::export]]
List lmC(const arma::vec & y, const arma::mat & X) {

    int n = X.n_rows, k = X.n_cols;

    // Centering
    arma::mat Xc=Xmcenter(X);
    arma::vec yc=y-arma::mean(y);

    arma::colvec coef = arma::solve(Xc, yc);
    arma::colvec resid = yc - Xc*coef;

    double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
    arma::colvec stderrest =
        arma::sqrt(sig2 * arma::diagvec( arma::pinv(arma::trans(Xc)*Xc)) );

   return List::create(Named("coefficients") = coef,
                       Named("stderr")       = stderrest);
}

// [[Rcpp::export]]
arma::colvec BMcgwa(const SEXP A,const arma::vec & y, const arma::uvec & vars) {
	// Read pointer
	Rcpp::XPtr<BigMatrix> bigMat(A);

	// Map to matrix
    arma::Mat<double> X0((double*) bigMat->matrix(), bigMat -> nrow(), bigMat -> ncol(), false);

	// Subset matrix
	X0=X0.cols(vars);

	// Centering
    arma::mat X=Xmcenter(X0);

    // Ordineary Least Squares
	// arma::colvec coef = arma::pinv(arma::trans(X)*X) * arma::trans(X) * arma::colvec(y);
	arma::colvec coef = solve(X,y);

 return(coef);
}

// [[Rcpp::export]]
arma::colvec BMmgwa(const SEXP A,
                    const arma::colvec & y,
                    const arma::uvec & vars,
                    bool debug=false) {


	// Read pointer
	Rcpp::XPtr<BigMatrix> bigMat(A);

	// Map to matrix
    arma::Mat<double> X0((double*) bigMat->matrix(), bigMat -> nrow(), bigMat -> ncol(), false);
	X0=X0.cols(vars);

	// Centering
    arma::mat X=Xmcenter(X0);

	// format change
	// y= arma::mat(y);

	// Ordineary Least Squares
  	arma::colvec coef(X.n_cols);
	arma::colvec val;


	// cout << "calculate effects" << endl;
	for(int i=0; i< X.n_cols; i++){
		arma::mat Xsub= X.col(i);
		// val = 1/arma::accu(Xsub) * arma::trans(Xsub) * (arma::mat(y));
		// val = 1/arma::accu(Xsub) * arma::trans(Xsub) * (arma::mat(y));
		// coef(i) = val(0);
		arma::vec val = solve(Xsub,arma::mat(y));
		if(debug) cout << val << endl;
		coef(i) = val(0);
	}

 return(coef);
}

// [[Rcpp::export]]
vec softmax_cpp(vec x, vec y) { // for utils in BMlasso
  return sign(x) % max(abs(x) - y, zeros(x.n_elem));
}

// [[Rcpp::export]]
arma::vec BMlasso(const SEXP A,const arma::colvec & y, const arma::uvec & vars,
             double lambda=1,double tol = 1e-5, int max_iter = 100){
	/*
	* Ridge regression
	*
	*  adapted from 	https://github.com/fditraglia/econ722/blob/master/RcppArmadillo/
	*/

	// Read pointer
	Rcpp::XPtr<BigMatrix> bigMat(A);

	// Map to matrix
	arma::Mat<double> X((double*) bigMat->matrix(),  bigMat -> nrow(), bigMat -> ncol(), true);

	// Subset matrix
	X=X.cols(vars);

	// Precompute some values
	int p = X.n_cols;
	mat XX = X.t() * X;
	vec Xy = X.t() * y;
	vec Xy2 = 2 * Xy;
	mat XX2 = 2 * XX;

	// Solver
	vec beta = solve(XX + diagmat(lambda * ones(p)), Xy);

	bool converged = false;
	int iteration = 0;
	vec beta_prev, aj, cj;

	while (!converged && (iteration < max_iter)){
		beta_prev = beta;
		for (int j = 0; j < p; j++){
			aj = XX2(j,j);
			cj = Xy2(j) - dot(XX2.row(j), beta) + beta(j) * XX2(j,j);
			beta(j) = as_scalar(softmax_cpp(cj / aj, lambda / aj));
		}
		iteration = iteration + 1;
		converged =  norm(beta_prev - beta, 1) < tol;
	}

	// return List::create(Named("beta") = beta,
	//             Named("n_iter") = iteration,
	//             Named("converged") = converged);
	return(beta);
}


// [[Rcpp::export]]
arma::vec BMsimridge(const SEXP A,const arma::colvec & y, const arma::uvec & vars,double lambda=1){
	/*
	* Ridge regression
	*
	*  adapted from 	https://github.com/fditraglia/econ722/blob/master/RcppArmadillo/
	*/

	// Read pointer
	Rcpp::XPtr<BigMatrix> bigMat(A);

	// Map to matrix
	arma::Mat<double> X0((double*) bigMat->matrix(),  bigMat -> nrow(), bigMat -> ncol(), true);

	// Subset matrix
	X0=X0.cols(vars);

	// center matrix
	arma::mat X=Xmcenter(X0);

	// Precompute some values
	int p = X.n_cols;
	mat XX = X.t() * X;
	vec Xy = X.t() * y;
	vec Xy2 = 2 * Xy;
	mat XX2 = 2 * XX;

	// Solver
	vec beta = solve(XX + diagmat(lambda * ones(p)), Xy);

	return(beta);
}

// [[Rcpp::export]]
List BMridge(const SEXP A,const arma::colvec & y, const arma::uvec & vars,
             const colvec & lambda) {
	/*
	* Ridge regression
	*
	* adapted from https://github.com/fditraglia/econ722/blob/master/RcppArmadillo/
	*/

	// Read pointer
	Rcpp::XPtr<BigMatrix> bigMat(A);

	// Map to matrix
	arma::Mat<double> X((double*) bigMat->matrix(),  bigMat -> nrow(), bigMat -> ncol(), true);

	// Subset matrix
	X=X.cols(vars);

	// Precompute some values
	int n_lam = lambda.n_elem;
	// int n = X.n_rows;
	int p = X.n_cols;

	mat coef(p, n_lam, fill::zeros);
	colvec y_tilde = join_cols(y, zeros<colvec>(p));


	for(int i = 0; i < n_lam; i++){
		mat X_tilde = join_cols(X, sqrt(lambda(i)) * eye(p, p));
		mat Q, R;
		qr_econ(Q, R, X_tilde);
		coef.col(i) = solve(R, Q.t() * y_tilde);
	}

	return List::create(Named("coef") = coef,
	          					Named("lambda") = lambda);
}

// [[Rcpp::export]]
List BMgwa1(const SEXP A,const arma::colvec & yraw,
	  const arma::uvec & vars,
	  const arma::uvec & training ,
	  int type=1,
      double lambda =1,
	  int max_iter=1000,
	  double tol = 1e-5
	  ){


   // Reading and subsetting
	 // arma::Mat<double> X0 = BMread(A);
   Rcpp::XPtr<BigMatrix> bigMat(A);
   arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);

   if(training.n_elem == X0.n_rows){
    X0= X0.cols(vars);
   }else{
    X0= X0.submat(training, vars);
   }
	// Centering
  arma::mat X=Xmcenter(X0);
  arma::vec y = (yraw/arma::mean(yraw) ) - 1;

	// Ordineary Least Squares
    int nsnps = X.n_cols;

  	arma::colvec coef(nsnps);

	switch(type){
		case 1:{
			cout << "Marginal GWA" << endl;
			arma::colvec val;
			for(int i=0; i< nsnps; i++){
				arma::mat Xsub= X.col(i);
				arma::vec val = solve(Xsub,arma::mat(y));
				coef(i) = val(0);
			}
			break;}
		case 2:{
			cout << "Conditional GWA" << endl;
			coef = solve(X,y);
			break;}
		case 3:{
			cout << "Conditional ridge GWA" << endl;
			mat XX = X.t() * X;
			vec Xy = X.t() * y;
			coef = solve(XX + diagmat(lambda * ones( nsnps)), Xy);
			break;}
		case 4:{
			cout << "Conditional lasso GWA" << endl;
			mat XX = X.t() * X;
			vec Xy = X.t() * y;
			vec Xy2 = 2 * Xy;
			mat XX2 = 2 * XX;

			vec beta = solve(XX + diagmat(lambda * ones( nsnps)), Xy);

			bool converged = false;
			int iteration = 0;
			vec beta_prev, aj, cj;

			while (!converged && (iteration < max_iter)){
				beta_prev = beta;
				for (int j = 0; j < nsnps; j++){
					aj = XX2(j,j);
					cj = Xy2(j) - dot(XX2.row(j), beta) + beta(j) * XX2(j,j);
					beta(j) = as_scalar(softmax_cpp(cj / aj, lambda / aj));
				}
				iteration = iteration + 1;
				converged =  norm(beta_prev - beta, 1) < tol;
			}
			coef= beta;
			break;}
	} // end switch

 // Get residuals
    int n = X.n_rows, k = X.n_cols;
    arma::colvec resid = y - X*coef;
    double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
    arma::colvec stderrest =
        arma::sqrt(sig2 * arma::diagvec( arma::pinv(arma::trans(X)*X)) );

   return List::create(Named("coefficients")    = coef,
                       Named("stderr")          = stderrest,
                       Named("meanasintercept") = arma::mean(yraw)
                       );
	// return(coef);
} // end gwa

// [[Rcpp::export]]
List BMgwa2(arma::mat X0,
            const arma::colvec & yraw,
        	  int type=1,
            double lambda =1,
        	  int max_iter=1000,
        	  double tol = 1e-5
	  ){

	// Centering
  arma::mat X=Xmcenter(X0);
  arma::vec y = (yraw/arma::mean(yraw) ) - 1;

	// Ordineary Least Squares
    int nsnps = X.n_cols;

  	arma::colvec coef(nsnps);

	switch(type){
		case 1:{
			cout << "Marginal GWA" << endl;
			arma::colvec val;
			for(int i=0; i< nsnps; i++){
				arma::mat Xsub= X.col(i);
				arma::vec val = solve(Xsub,arma::mat(y));
				coef(i) = val(0);
			}
			break;}
		case 2:{
			cout << "Conditional GWA" << endl;
			coef = solve(X,y);
			break;}
		case 3:{
			cout << "Conditional ridge GWA" << endl;
			mat XX = X.t() * X;
			vec Xy = X.t() * y;
			coef = solve(XX + diagmat(lambda * ones( nsnps)), Xy);
			break;}
		case 4:{
			cout << "Conditional lasso GWA" << endl;
			mat XX = X.t() * X;
			vec Xy = X.t() * y;
			vec Xy2 = 2 * Xy;
			mat XX2 = 2 * XX;

			vec beta = solve(XX + diagmat(lambda * ones( nsnps)), Xy);

			bool converged = false;
			int iteration = 0;
			vec beta_prev, aj, cj;

			while (!converged && (iteration < max_iter)){
				beta_prev = beta;
				for (int j = 0; j < nsnps; j++){
					aj = XX2(j,j);
					cj = Xy2(j) - dot(XX2.row(j), beta) + beta(j) * XX2(j,j);
					beta(j) = as_scalar(softmax_cpp(cj / aj, lambda / aj));
				}
				iteration = iteration + 1;
				converged =  norm(beta_prev - beta, 1) < tol;
			}
			coef= beta;
			break;}
	} // end switch

 // Get residuals
    int n = X.n_rows, k = X.n_cols;
    arma::colvec resid = y - X*coef;
    double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
    arma::colvec stderrest =
        arma::sqrt(sig2 * arma::diagvec( arma::pinv(arma::trans(X)*X)) );

   return List::create(Named("coefficients")    = coef,
                       Named("stderr")          = stderrest,
                       Named("meanasintercept") = arma::mean(yraw)
                       );
	// return(coef);
} // end gwa


// List BMgwa2(const SEXP A,const arma::colvec & yraw,
// 	  const arma::uvec & vars,
// 	  const arma::uvec & training ,
// 	  int type=1,
//       double lambda =1,
// 	  int max_iter=1000,
// 	  double tol = 1e-5
// ){
//
// }


// // [[Rcpp::export]]
// arma::vec BMgwa(const SEXP A,const arma::colvec & yraw,
// 	  const arma::uvec & vars,
// 	  int type=1,
//       double lambda =1,
// 	  int max_iter=1000,
// 	  double tol = 1e-5
// 	  ){
//
// 	// Reading and subsetting
//   arma::Mat<double> X0= BMcolsubset(A,vars);
//
// 	// Centering
//   arma::mat X=Xmcenter(X0);
//   arma::vec y = (yraw/arma::mean(yraw) ) - 1;
//
//
// 	// Ordineary Least Squares
//   int nsnps = X.n_cols;
//
//   	arma::colvec coef(nsnps);
//
// 	switch(type){
// 		case 1:{
// 			cout << "Marginal GWA" << endl;
// 			arma::colvec val;
// 			for(int i=0; i< nsnps; i++){
// 				arma::mat Xsub= X.col(i);
// 				arma::vec val = solve(Xsub,arma::mat(y));
// 				coef(i) = val(0);
// 			}
// 			break;}
// 		case 2:{
// 			cout << "Conditional GWA" << endl;
// 			coef = solve(X,y);
// 			break;}
// 		case 3:{
// 			cout << "Conditional ridge GWA" << endl;
// 			mat XX = X.t() * X;
// 			vec Xy = X.t() * y;
// 			coef = solve(XX + diagmat(lambda * ones( nsnps)), Xy);
// 			break;}
// 		case 4:{
// 			cout << "Conditional lasso GWA" << endl;
// 			mat XX = X.t() * X;
// 			vec Xy = X.t() * y;
// 			vec Xy2 = 2 * Xy;
// 			mat XX2 = 2 * XX;
//
// 			vec beta = solve(XX + diagmat(lambda * ones( nsnps)), Xy);
//
// 			bool converged = false;
// 			int iteration = 0;
// 			vec beta_prev, aj, cj;
//
// 			while (!converged && (iteration < max_iter)){
// 				beta_prev = beta;
// 				for (int j = 0; j < nsnps; j++){
// 					aj = XX2(j,j);
// 					cj = Xy2(j) - dot(XX2.row(j), beta) + beta(j) * XX2(j,j);
// 					beta(j) = as_scalar(softmax_cpp(cj / aj, lambda / aj));
// 				}
// 				iteration = iteration + 1;
// 				converged =  norm(beta_prev - beta, 1) < tol;
// 			}
// 			coef= beta;
// 			break;}
// 	} // end switch
//
// 	return(coef);
// } // end gwa



// arma::colvec BMmgwa(const SEXP A,const arma::colvec & y, const arma::uvec & vars, bool debug=false) {
// 	// Read pointer
// 	Rcpp::XPtr<BigMatrix> bigMat(A);
//
// 	if(bigMat->matrix_type() !=8) stop("Big matrix is not of type double");
//
//
// 	// Map to matrix
//     arma::Mat<double> X((double*) bigMat->matrix(),
//      bigMat -> nrow(),  bigMat -> ncol(), true);
//
// 	// Subset matrix
// 	X=X.cols(vars);
//
// 	// Mean center
// 	X=Xmcenter(X);
//
//     // Ordineary Least Squares
//     arma::colvec coef(vars.n_elem);
// 	arma::colvec val;
//
// 	for(int i=0; i< vars.n_elem; i++){
// 		arma::mat Xsub= X.col(i);
// 		val = 1/arma::accu(Xsub) * arma::trans(Xsub) * arma::mat(y);
// 		coef(i) = val(0);
// 	}
//
//  return(coef);
// }


// // [[Rcpp::export]]
// arma::colvec BMmgwa(const SEXP A,
//                     const arma::colvec & Y,
//                     const arma::colvec & h,
//                     const arma::uvec & vars,
//                     bool debug=false){
// 	// Read pointer
// 	Rcpp::XPtr<BigMatrix> bigMat(A);
// 	arma::mat X;
// 		switch(bigMat->matrix_type()){
// 			case 1:
// 		    throw Rcpp::exception("unknown type detected for big.matrix object!");
//   		case 2:
// 		    throw Rcpp::exception("unknown type detected for big.matrix object!");
// 		case 4:{
// 				cout << "Big matrix is not of type double. Attempting transformation..." << endl;
// 	   		arma::Mat<int> Xi((int*) bigMat->matrix(), bigMat -> nrow(),  bigMat -> ncol(), true);
// 	 	 		// Subset matrix
// 		 		Xi=Xi.cols(vars);
// 		 		X = arma::conv_to< arma::Mat<double> >::from(Xi);
// 		 		return(BMmgwa(X,Y,h,debug));
// 		  }
// 		case 8:{
// 				arma::Mat<double> X1((double*) bigMat->matrix(), bigMat -> nrow(),  bigMat -> ncol(), true);
// 	 	 		// Subset matrix
// 		 		X=X1.cols(vars);
// 		 		return(BMmgwa(X,Y,h,debug));
// 			}
// 		}
// }


// // [[Rcpp::export]]
// arma::colvec BMcgwa_mm(const SEXP A,const arma::vec & y, const arma::uvec & m) {
// 		// Read pointer
// 		Rcpp::XPtr<BigMatrix> bigMat(A);
//
// 		int nrows = bigMat -> nrow();
// 		int ncols = bigMat -> ncol();
//
// 		// Map to matrix
//     arma::Mat<int> Am((int*) bigMat->matrix(), nrows, ncols, true); // slower when it is true
//
// 		// Subset matrix
// 		Am=Am.cols(m);
//
//     // Convert to double for linear algebra operations
//     arma::Mat<double> X = arma::conv_to< arma::Mat<double> >::from(Am);
//
//     // Ordineary Least Squares
// 		arma::colvec coef = arma::inv(arma::trans(X)*X) * arma::trans(X) * arma::colvec(y);
//
//  return(coef);
// }

//

  // arma::colvec coef = arma::solve(X,y);
  // arma::vec coef = arma::inv(X.t()*X) * X.t() * y;
  // arma::vec coef = (X.t()*X) * X.t() * y;
  // arma::colvec coef = arma::inv(arma::trans(X)*X) * arma::trans(X) * arma::colvec(y);
  // int j=0;
  // arma::colvec coef = 1/(arma::trans((X.col(j)))*(X.col(j))) * arma::trans((X.col(j))) * arma::colvec(y);
  // arma::colvec coef;
  // double val;
  // cout << X.col(j) << endl;
  // cout << (arma::trans((X.col(j)))*(X.col(j)) << endl;
  // cout << arma::colvec(y)  << endl;
  // cout << arma::trans((X.col(j))) * arma::colvec(y) << endl;
  // arma::colvec val = (1/(arma::trans((X.col(j)))*(X.col(j))) * arma::trans((X.col(j))) * arma::colvec(y));
  // coef = join_rows( coef, val);


//
// // [[Rcpp::export]]
// arma::Mat<int> BMsub(SEXP A, arma::vec , arma::uvec mycols){
//       Rcpp::XPtr<BigMatrix> bigMat(A);
//       arma::Mat<int> X((int*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
//
//       return(X.cols(mycols));
// }







// // [[Rcpp::export]]
// arma::colvec BMmgwa(SEXP A, arma::vec y) {

//   Rcpp::XPtr<BigMatrix> bigMat(A);
//   MatrixAccessor<int> Am(*bigMat);
//   int nrows = bigMat->nrow();
//   int ncolumns = bigMat->ncol();

//   arma::Mat<int> X((int*) bigMat->matrix(), bigMat->nrow(),bigMat->ncol(), false, false);

//   arma::colvec coef(bigMat->ncol());
//   double val;
//   int j=0;
//   // for (int j = 0; j < ncolumns; j++){

//      val = 1/(arma::trans((X.col(j)))*(X.col(j))) * arma::trans((X.col(j))) * arma::colvec(y);
//      // coef(j)=val;
//   // }

// // arma::colvec coef = arma::solve(X,y);
// // arma::vec coef = arma::inv(X.t()*X) * X.t() * y;
// // arma::vec coef = (X.t()*X) * X.t() * y;
// // arma::colvec coef = arma::inv(arma::trans(X)*X) * arma::trans(X) * arma::colvec(y);
// // int j=0;
// // arma::colvec coef = arma::inv(arma::trans(arma::mat(X.row(j)))*arma::mat(X.row(j))) * arma::trans(arma::mat(X.row(j))) * arma::colvec(y);


// // arma::colvec coef(ncolumns);
// // int j=0;
// // const arma::Row<uint32_t> j=0;
// // for (int j = 0; j < ncolumns; j++){
// //      for (int i = 0; i < nrows; i++){

//        // auto Xs=X.rows(j);
//        // auto Xst=X.rows(j).t();
//        // auto V = Xst * Xs;
//        // arma::mat Vinv = arma::pinv(V,0.01);
//        // // arma::mat V=X.rows(j).t() * X.rows(j);
//        // arma::rowvec Xst=X.rows(j);
//        // arma::mat V=X.row.t() * X.row(j);
//        // arma::mat Vinv=arma::pinv(V);
//        // coef(j) = Vinv *  X.row(j) * y;
//          // bgwa=sapply(1:M,function(m) solve(t(X[,m])%*%X[,m]) %*% t(X[,m]) %*% Y)

//  //     }
//  // }
//   return(val);
// }


// arma::vec prodArmaSub(XPtr<BigMatrix> xpA, const arma::vec& x,
//                       const arma::Row<uint32_t>& ind) {
//   arma::Mat<char> Am((char *) xpA->matrix(), xpA->nrow(), xpA->ncol(), false);

//   return Am.rows(ind) * x;


// // [[Rcpp::export]]
// arma::vec mgwa(arma::mat A, arma::vec y) {

  // arma::colvec coef(ncolumns);
  // int j=0;
  // for (int j = 0; j < ncolumns; j++){
       // for (int i = 0; i < nrows; i++){
         // arma::vec coef = arma::spsolve(A, y);  // solve one system
         // arma::vec coef = arma::solve(A, y,"lapack");
        // bgwa=sapply(1:M,function(m) solve(t(X[,m])%*%X[,m]) %*% t(X[,m]) %*% Y)

       // }
   // }
// sp_mat A = sprandu<sp_mat>(1000, 1000, 0.1);

// vec b = randu<vec>(1000);
// mat B = randu<mat>(1000, 5);

// vec x = spsolve(A, b);  // solve one system
// mat X = spsolve(A, B);  // solve several systems

// return(coef);
// }
