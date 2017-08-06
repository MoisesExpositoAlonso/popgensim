#include <stdio.h>
#include <iostream>
#include <Rcpp.h> // when armadillo is loaded, remove this
// #include <RcppArmadillo.h> // when armadillo is loaded, remove #include <Rcpp.h>

// Set seed
Rcpp::RNGScope scope;


using namespace Rcpp;

///// #include <RcppArmadillo.h>
///// [[Rcpp::depends(RcppArmadillo)]]


//[[Rcpp::export]]
NumericVector allelesimCvec(
    double mu=0,
    double nu=0,
    double m=0,
    double wAA=1,
    double wAa=0.75,
    double waa=0.5,
    double p0=0.5,
    double psource=0.1,
    double Fi=0,
    double d=0,
    int N=100,
    int tmax=10
  ){
  NumericVector p(tmax);

  p(0) = p0;

  for (int t = 0; t < tmax-1; t++) {

  // introduce mutations
  p(t+1) = (1-mu)*p(t) + (1-p(t))*nu;

  // introduce migrations
  p(t+1) = (1-m)*p(t+1) + m*psource;

  // selection
  double w_hat = ( wAA*pow(p(t+1),2.0)+ wAa*2*p(t+1)*(1-p(t+1)) + waa*pow(1-p(t+1),2.0));
  if ( std::isnan(w_hat) || w_hat==0 ) {
  p(t+1) =   NumericVector::get_na()  ;
  }
  else {
  p(t+1) =  ( wAA*pow(p(t+1),2.0)+ wAa*p(t+1)*(1-p(t+1)) ) / w_hat;
  }
  // p(t+1) =  ( wAA*pow(p(t+1),2.0)+ wAa*p(t+1)*(1-p(t+1)) ) / w_hat+0.00000000001;
          // the above code will never produce NA or 0

  // inbreeding

  double fAA = (pow(p(t+1),2.0) * (1-Fi)) + p(t+1)*Fi;
  double fAa = 2*p(t+1)*(1-p(t+1)) * (1-Fi);
  double faa = pow(1-p(t+1),2.0) * (1-Fi) + (1-p(t+1))* Fi;

  // finite populations
  double NAA = round(N*fAA);
  double NAa = round(N*fAa);
  double Naa = round(N*faa);

  // random death
  if(d != 0){
    NAA = R::rbinom(NAA, 1-d) ;
    NAa = R::rbinom(NAa, 1-d) ;
    Naa = R::rbinom(Naa, 1-d) ;
  }

  // output
  p(t+1) = (2*NAA + NAa) / (2*(NAA+NAa+Naa));

  }

  return p;

}



//[[Rcpp::export]]
NumericMatrix allelesimCmat(
    double mu=0,
    double nu=0,
    double m=0,
    double wAA=1,
    double wAa=0.75,
    double waa=0.5,
    double p0=0.5,
    double psource=0.1,
    double Fi=0,
    double d=0,
    int N=100,
    int tmax=10,
    int rep=10
  )
{

  // Define data structure objects
  NumericMatrix p(rep,tmax);

 for(int i=0;i<rep; i++){
  p(i,_) = allelesimCvec(mu,nu,m,wAA,wAa,waa,p0,psource,Fi,d,N, tmax);
 }

return p;
}





// //[[Rcpp::export]]
// NumericMatrix allelesimC(
//     double mu=0,
//     double nu=0,
//     double m=0,
//     double wAA=1,
//     double wAa=0.75,
//     double waa=0.5,
//     double p0=0.5,
//     double psource=0.1,
//     double Fi=0,
//     int N=100,
//     int tmax=10,
//     int rep=10
//   )
// {
//   NumericMatrix p(rep,tmax);
//
//   for(int i=0; i<rep;i++){
//     p(i,0) = p0;
//   }
// // `
//   for (int r=0; r<rep; r++){
//   for (int t = 1; t < tmax; t++) {
// //
// //   // introduce mutations
//   double porigin=p(r,t) ;
//   double pp = (1-mu)*porigin + (1-porigin)*nu;
// //
// //   // introduce migrations
// //   double ppp = (1-m)*pp + m*psource;
// //
// //   // selection
// //   double w_hat = ( wAA*pow(ppp,2.0)+ wAa*2*ppp*(1-ppp) + waa*pow(1-ppp,2.0));
// //     if ( std::isnan(w_hat) ) {
// //     p(r,t+1) =  NumericVector::get_na() ;
// //     }
// //     else {
// //     p(r,t+1) =  ( wAA*pow(ppp,2.0)+ wAa*ppp*(1-ppp) ) / w_hat;
// //     }
// //
// //   // inbreeding
// //   double pppp = p(r,t+1);
// //
// //   double fAA = (pow(pppp,2.0) * (1-Fi)) + pppp*Fi;
// //   double fAa = 2*pppp*(1-pppp) * (1-Fi);
// //   double faa = pow(1-pppp,2.0) * (1-Fi) + (1-pppp)* Fi;
// //
// //   // finite populations
// //   int NAA = round(N*fAA);
// //   int NAa = round(N*fAa);
// //   int Naa = round(N*faa);
// //
// //   // output
// // p(r,t+1) = (2*NAA + NAa) / (2*(NAA+NAa+Naa));
// // p(r,t+1) = (2*NAA + NAa) / (2*(NAA+NAa+Naa));
//    p(r,t+1) = pp;
//
//   }}
//
//   return p;
//
// }



// if (NAA <= 0) {  // This could help to kill part of the population
//   NAA <- 0
//   NAA.prime <- 0
// } else {
//   NAA.prime <- sum(rbinom(n=NAA, size=1, prob=d))
// }
// if (NAa <= 0) {
//   NAa <- 0
//   NAa.prime <- 0
// } else {
//   NAa.prime <- sum(rbinom(n=NAa, size=1, prob=d))
// }
// if (Naa <= 0) {
//   Naa <- 0
//   Naa.prime <- 0
// } else {
//   Naa.prime <- sum(rbinom(n=Naa, size=1, prob=d))
// }
// N.prime <- NAA.prime + NAa.prime + Naa.prime
//
// if (N.prime <= 0) {
//   p[t+1] <- NA
// } else {
//   p[t+1] <- (2*NAA.prime + NAa.prime) / (2*N.prime)
// }





// // [[Rcpp::export]]
// int main()
// {
//    // printf() displays the string inside quotation
//    printf("Hello, World!");
//    return 0;
// }

// int main() {
//   allelesimC();
//
//   return 0;
// }

//
// /*** R
// library(microbenchmark)
// x <- runif(1e5)
// microbenchmark(
//   mean(x),
//   meanC(x)
// )
// */
