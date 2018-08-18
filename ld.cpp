#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
#include <chrono>
#include <ctime>

#include <cstdio>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <string>

#define ARMA_64BIT_WORD 1


// when armadillo is loaded, remove this below
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;


#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(bigmemory)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


////////////////////////////////////////////////////////////////////////////////
//// Declarations
////////////////////////////////////////////////////////////////////////////////

arma::vec upperTmat(const arma::mat mat);

////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat Xmvcenter(arma::mat X);

// [[Rcpp::export]]
arma::mat LDrelative(SEXP A, arma::uvec  m, bool debug = false){

  Rcpp::XPtr<BigMatrix> bigMat(A);
  if(bigMat->matrix_type() !=8) stop("Big matrix is not of type double");

  // Read the genome matrix from address
  arma::Mat<double> X((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
  X=X.cols(m);

  // mean and var center for LD calculation
  X=Xmvcenter(X);
  if(debug) cout << X << endl;

  // Get the relative LD for the proposals
  arma::mat R2 =  arma::trans(X)*X ;
  if(debug)  cout << R2 << endl;
  R2 = R2/ arma::sum(upperTmat(R2));
  if(debug)  cout << arma::sum(upperTmat(R2)) << endl;
  return(R2);
}

arma::mat LDrelative(arma::mat X, bool debug = false){

  // mean and var center for LD calculation
  X=Xmvcenter(X);
  if(debug) cout << X << endl;

  // Get the relative LD for the proposals
  arma::mat R2 =  arma::trans(X)*X ;
  if(debug)  cout << R2 << endl;
  R2 = R2/ arma::sum(upperTmat(R2));
  if(debug)  cout << arma::sum(upperTmat(R2)) << endl;
  return(R2);
}
