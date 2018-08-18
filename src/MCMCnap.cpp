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
/// Profiling utilities
////////////////////////////////////////////////////////////////////////////////

// RcppExport SEXP start_profiler(SEXP str) {
//   ProfilerStart(as<const char*>(str));
//   return R_NilValue;
// }

// RcppExport SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }


////////////////////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////////////////////



#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}


template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

template <typename T>
arma::vec vec2arma(const T& x) {
  return  Rcpp::as<arma::vec>(x);
}

// [[Rcpp::export]]
NumericVector sample_num( NumericVector x,
                          int size,
                          bool replace,
                          NumericVector prob = NumericVector::create()
)
{
  NumericVector ret = RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}


// #define MIN_NUM = std::numeric_limits<float>::min(); // problem is that it does not know the type
const double MIN_NUM = std::numeric_limits<float>::min();
// #define PIraw = 3.14159265358979323846;
// const double PI= 3.14159265358979323846;
// # define PI 3.14159265358979323846  /* pi */
//const double MAX_NUM = std::numeric_limits<float>::max();


////////////////////////////////////////////////////////////////////////////////
/// Matrix and LD utilities
////////////////////////////////////////////////////////////////////////////////


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
arma::vec upperTmat(const arma::mat mat){
  arma::vec output((mat.n_cols*(mat.n_cols-1))/2);
  arma::mat::const_iterator it = mat.begin() + mat.n_rows; //iterator at rows to skip in every step, starts at second column
  long toSkipInVec = 0;
  for(int i = 1; i < mat.n_cols; i++) //Starts with 1 to skip the diagonal
  {
    std::copy(it, it + i, output.begin() + toSkipInVec);
    toSkipInVec += i;
    it += mat.n_rows;
  }
  return output;
}

// [[Rcpp::export]]
arma::mat Xmvcenter(arma::mat X){
  arma::mat newX(X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols; j++){
   newX.col(j) = (X.col(j) - arma::mean( X.col(j))) /arma::stddev( X.col(j));
  }
  return(newX);
}

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



////////////////////////////////////////////////////////////////////////////////
/// Likelihood, Probabilities, Proposals
////////////////////////////////////////////////////////////////////////////////


////// Utils proposals //////

// [[Rcpp::export]]
double runif_reflect(double minhere,double maxhere,double min,double max){
  // int counter=1;
  double newval;
  if( min == max){
    newval= min; // if values want to be kept constant
  }else{
    newval =Rcpp::runif(1,minhere,maxhere)(0);
    if(newval<min){newval = (min-newval) + min;}
    else if(newval>max){newval = max- (newval-max);}
  }
  if(newval < min || newval>max){   // Check if it is out of bounds
      newval=(max-min)/2;
  }
  return(newval);
}

////////////////////////////////////////////////////////////////////////////////
////// Hyperparameter Proposals //////

class GPROPOSAL{
  private:
    double b; double bmin; double bmax;
    double a; double amin; double amax;
    double p; double pmin=0; double pmax=1;
    double mu; double mumin; double mumax;
    double epi; double epimin; double epimax;
    double svar; double svarmin; double svarmax;
    double ss; double ssmin; double ssmax;
    int nupdates=1;
    bool verbose=false;
    double bw=0.1;
  public:
    GPROPOSAL(
            double b_=0.5,double bmin_=0,double bmax_=1,
            double a_=0.1,double amin_=0,double amax_=1,
            double p_=0.5,
            double mu_=1,double mumin_=0, double mumax_=50,
            double epi_=1,double epimin_=0, double epimax_=5,
            double svar_=0.5,double svarmin_=0, double svarmax_=5,
            double ss_=0.1,double ssmin_=0, double ssmax_=1
            ){
      b=b_;bmin=bmin_;bmax=bmax_;
      a=a_;amin=amin_;amax=amax_;
      p=p_;
      mu=mu_;mumin=mumin_;mumax=mumax_;
      epi=epi_;epimin=epimin_;epimax=epimax_;
      svar=svar_;svarmin=svarmin_;svarmax=svarmax_;
      ss=ss_;ssmin=ssmin_;ssmax=ssmax_;
    }
    void setbandwidth(double bandwidth){bw=bandwidth;}
    void setupdatesnum(int ups){nupdates=ups;}
    void setverbose(bool verbose_){verbose=verbose_;}
    void printatributes(){
      cout <<"bw = " << bw << endl;
      cout <<"b = " << b << " [" << bmin << " " << bmax << "]" << endl;
      cout <<"a = " << a << " [" << amin << " " << amax << "]" << endl;
      cout <<"p = " << p << " [" << pmin << " " << pmax << "]" << endl;
      cout <<"mu = " << mu << " [" << mumin << " " << mumax << "]" << endl;
      cout <<"epi = " << epi << " [" << epimin << " " << epimax << "]" << endl;
      cout <<"svar = " << svar << " [" << svarmin << " " << svarmax << "]" << endl;
      cout <<"ss = " << ss << " [" << ssmin << " " << ssmax << "]" << endl;
    }
    arma::vec fn(arma::vec g){
      // New proposal
      arma::vec news=g;
      // Update one position
      double minhere,maxhere;
      double newval;
      if(verbose) cout << "Loop to substitute position" << endl;
      for(int j=0; j< nupdates; j++){
        int randomIndex = rand() % g.size();
        switch(randomIndex){
            if(verbose) cout << g(randomIndex) << endl;
            case 0:
              minhere=g(randomIndex)- (bw *(bmax-bmin)) ;
              maxhere=g(randomIndex)+ (bw *(bmax-bmin));
              newval= runif_reflect(minhere,maxhere,bmin,bmax); break;
            case 1:
              minhere=g(randomIndex)- (bw *(amax-amin)) ;
              maxhere=g(randomIndex)+ (bw *(amax-amin));
              newval= runif_reflect(minhere,maxhere,amin,amax); break;
            case 2:
              minhere=g(randomIndex)- (bw *(pmax-pmin)) ;
              maxhere=g(randomIndex)+ (bw *(pmax-pmin));
              newval= runif_reflect(minhere,maxhere,pmin,pmax); break;
            case 3:
              minhere=g(randomIndex)- (bw *(mumax-mumin)) ;
              maxhere=g(randomIndex)+ (bw *(mumax-mumin));
              newval= runif_reflect(minhere,maxhere,mumin,mumax); break;
            case 4:
              minhere=g(randomIndex)- (bw *(epimax-epimin)) ;
              maxhere=g(randomIndex)+ (bw *(epimax-epimin));
              newval= runif_reflect(minhere,maxhere,epimin,epimax); break;
            case 5:
              minhere=g(randomIndex)- (bw *(svarmax-svarmin)) ;
              maxhere=g(randomIndex)+ (bw *(svarmax-svarmin));
              newval= runif_reflect(minhere,maxhere,svarmin,svarmax); break;
            case 6:
              minhere=g(randomIndex)- (bw *(ssmax-ssmin)) ;
              maxhere=g(randomIndex)+ (bw *(ssmax-ssmin));
              newval= runif_reflect(minhere,maxhere,ssmin,ssmax); break;
        }
        if(verbose) cout << newval << endl;
        news(randomIndex) = newval;
      }
      if(verbose) cout << "End loop" << endl;
      return(news);
  }
};


// [[Rcpp::export]]
arma::mat test_GPROPOSAL(double b=1,
                          double a=1,
                          double p=1,
                          double mu=1,
                          double epi=1,
                          double svar=1,
                          double ss=0.1,
                          double bw=0.5,
                          int iter=1000,
                          bool verbose = false,
                          int nupdates=7
                          ){

  arma::vec g(7);
  g(0)= b;
  g(1)= a;
  g(2)= p;
  g(3)= mu;
  g(4)= epi;
  g(5)= svar;
  g(6)= ss;
  GPROPOSAL GProp;

  GProp.setupdatesnum(nupdates);
  GProp.setbandwidth(bw);
  GProp.setverbose(verbose);
  GProp.printatributes();

  arma::mat res(7,iter);
  res.col(0)=g;

  cout << "Runing several proposals "<< endl;
  for(int i=1; i<iter; i++){
    g=GProp.fn(g);
    res.col(i)=g;
  }
  return(res);
}


////////////////////////////////////////////////////////////////////////////////
////// Hyperparameter Proposals //////

class SPROPOSAL {
  private:
    int nupdates;
    double min;
    double max;
    int mode;
    bool verbose;
    double bw=0.5;
    arma::mat R=arma::randu<arma::mat>(1,1);  // for LD implementatioon
  public:
  SPROPOSAL(
            int nupdates_,
            int mode_=1,
            double min_=-0.98039215686274505668,
            double max_=50,
            bool verbose_= false
            ){
   nupdates=nupdates_;min=min_;max=max_;mode=mode_;verbose=verbose_;
    }
  void setR2(arma::mat R2){R=R2;}
  void setbandwidth(double bandwidth){bw=bandwidth;}
  void printatributes(){
    cout <<"bw = " << bw << endl;
    cout <<"nupdates = " <<nupdates << endl;
    cout <<"min = " << min << endl;
    cout <<"max = " <<max << endl;
    cout <<"mode = " <<mode << endl;
    cout <<"verbose = " <<verbose << endl;
  }
  arma::vec fn(arma::vec s, double svar=0.1){
    switch(mode){
      case 1:
        return update(s);
        break;
      case 2:
        return updateLD(s);
      case 3:
        return updatelog(s,svar);
        break;
      default:
        return update(s);
        break;
    }
  }
   arma::vec update(arma::vec s){
      arma::colvec news=s;
      double minhere,maxhere,newval;
      if(verbose) cout << "Loop to substitute position" << endl;
      for(int j=0; j< nupdates; j++){
        int randomIndex = rand() % s.size();
        minhere=s(randomIndex)-bw;
        maxhere=s(randomIndex)+bw;
        newval = runif_reflect(minhere,maxhere,min,max);
        news(randomIndex) = newval;
      }
      if(verbose) cout << "End loop" << endl;
      return(news);
  }
    arma::vec updatelog(arma::vec s,double svar){
      arma::colvec news= s;
      double meanhere,newval;
      if(verbose) cout << "Loop to substitute position" << endl;
      for(int j=0; j< nupdates; j++){
        int randomIndex = rand() % s.size();
        meanhere=log(1+s(randomIndex)); //*** transform mean to log dimension
        // newval = Rcpp::rnorm(1,1,svar)(0); // this will always draw from the same distribution. The chain is not reversible anymore and the aceptance ratio has to be updated
        //*// if(std::isinf(meanhere)) meanhere= (max-min)/2; // check for infinity
        if(std::isinf(meanhere)) meanhere= 0; // check for infinity
        newval = Rcpp::rnorm(1,meanhere,svar)(0); // sample with the mean being the last value.
        //*// if(newval<log(1+min)) newval= log(1+min)+ (newval-log(1+min)); // refract in limits
        //*// if(newval>log(1+max)) newval= log(1+max)- (newval-log(1+max)); // refract in limits
        // NOTE: if any of the lines //*// the chain gets biased towards lower values
        news(randomIndex) = exp(newval)-1; //*** back transform to natural dimension
        if(verbose) cout << newval << endl;
        if(verbose) cout << news(randomIndex) << endl;
      }
      if(verbose) cout << "End loop" << endl;
      return(news);
  }
  arma::vec updateLD(arma::vec s){
      arma::colvec news=s;
      double minhere,maxhere,newval;
      if(verbose) cout << "Loop to substitute position" << endl;
      int randomIndex = rand() % s.size();
      minhere=s(randomIndex)-bw;
      maxhere=s(randomIndex)+bw;
      newval = runif_reflect(minhere,maxhere,min,max);
      news(randomIndex) = newval;
      double diff = news(randomIndex) - s(randomIndex);
      if(verbose) cout << "Difference with original value = " << diff << endl;
      for(int i=0; i<s.n_elem & i!= randomIndex; i++){
        news(i) -= diff * R(randomIndex,i) ;
        if(verbose) cout << news(i) << endl;
          if(news(i) < min) news(i) =min+abs(news(i) -min);
          if(news(i) > max) news(i) =max-abs(news(i) -max);
      }
      return(news);
  }
};

// [[Rcpp::export]]
arma::vec PropoS( // I think this function might be unnecessary
                 int nupdates=1,
                 double svar = 0.5){
  arma::vec s(nupdates);
  // s.fill(0);
  SPROPOSAL ps(nupdates,svar);
  return(ps.fn(s,svar));
}

// [[Rcpp::export]]
arma::mat test_SPROPOSAL(int nupdates=1,
                         double svar =0.5,
                         int m=3,
                         int iter=3000,
                         int mode=3
                         ){

  // initialize class
  arma::vec s(m);
  s.fill(0);
  SPROPOSAL PS(nupdates,mode);

  // create result matrix
  arma::mat res(m,iter);
  res.col(0)=s;

  cout << "Runing several proposals "<< endl;
  for(int i=1; i<iter; i++){
    s=PS.fn(s,svar);
    res.col(i)=s;
  }

  return(res);
}

////////////////////////////////////////////////////////////////////////////////
////// PRIOR

class PRIOR{
public:
  double min;
  double max;
  double mean;
  double svar;
  double ss;
  int mode;
  PRIOR(double par1=0,double par2=1, int mode_=1){
       mode=mode_;
      switch(mode){
        case 1: // moc mode, return 1
          break;
        case 2: // true uniform
          min=par1;
          max=par2;
          break;
        case 3: // log +1 normal
          mean=par1;
          svar=par2;
          break;
        case 4: // log +1 mixture normal with sparcity
          svar=par1;
          ss=par2;
          break;
        default:
          min=0;max=1;mean=0,svar=0.5;
          break;
    }
  }
  void printatributes(){
    cout <<"min = " << min << endl;
    cout <<"max = " <<max << endl;
    cout <<"s variance = " <<svar << endl;
    cout <<"mode = " <<mode << endl;
  }
// Prior probability functions
  double uniform(const arma::colvec & s){
    double L= 0;
    int N=s.n_elem;
    for(int i=0;i<N ;i++){
      L+= R::dunif(s(i),min,max,true);
    }
    return L;
  }
  double loggaussian(const arma::colvec & s, double svar){ // int n=s.n_elem;
    arma::vec x = log(1+s); // double L = -.5*n*log(2*PI) -.5*n*log(svar) -(1/(2*svar))*sum(arma::pow((x-mean),2));
    double L=0;
    for(int i = 0; i<s.n_elem; i++){
      L+= R::dnorm(x(i),0,svar,true);
    }
    return L;
  }
double logmixgaussian(const arma::colvec & s, double svar, double ss){
    arma::vec x = log(1+s);
    double L=0;
    for(int i = 0; i<s.n_elem; i++){
      if(x(i)==0){
        L += log(ss  + (1-ss) * R::dnorm(x(i),0,svar,false)) ;
      }else{
        L += (1-ss) * R::dnorm(x(i),0,svar,true);  // L+= R::dnorm(x(i),0,svar,true);
      }
    }
    return L;
  }
// Prior wrapper
  double fn(const arma::colvec & s){
    switch(mode){
      case 1: // moc mode, return 1
        return 1.0; break;
      case 2: // true uniform
        return uniform(s); break;
      default:
        return 1.0; break;
    }
  }
  double fn(const arma::colvec & s,const double & svar,const double & ss){
        switch(mode){
      case 1: // moc mode, return 1
        return 1.0; break;
      case 2: // true uniform
        return uniform(s); break;
      case 3: // logGaussian
        return loggaussian(s,svar); break;
      case 4: // mixture distribution point mass zero and logGaussian
        return logmixgaussian(s,svar,ss); break;
      default:
        return 1.0; break;
    }
  }
};


// [[Rcpp::export]]
arma::vec test_PRIOR(arma::vec s,
                      double min=-1,
                      double max=10,
                      double mean=0,
                      double svar=0.5,
                      double sparsity=0.1,
                      int mode=1
                    ){
  arma::vec res(4);

  PRIOR Pri; // mode 1 = uniform moc
  // cout << "Prior mode = 1" << endl;
  // cout << Pri.fn(s) << endl;
  res.row(0) =  Pri.fn(s);

  PRIOR Pri2(min,max,2); // mode 1 = uniform moc
  // cout << "Prior mode = 2" << endl;
  // cout << Pri2.fn(s) << endl;
  res.row(1) =  Pri2.fn(s);

  PRIOR Pri3(mean,svar,3); // mode 1 = uniform moc
  // cout << "Prior mode = 3" << endl;
  // cout << Pri3.fn(s,svar,sparsity) << endl;
  res.row(2) = Pri3.fn(s,svar,sparsity) ;

  PRIOR Pri4(mean,svar,4); // mode 1 = uniform moc
  // cout << "Prior mode = 4" << endl;
  // cout << Pri4.fn(s,svar,sparsity) << endl;
  res.row(3) = Pri4.fn(s,svar,sparsity);

  return(res);
}

////////////////////////////////////////////////////////////////////////////////
/// Fitness functions
////////////////////////////////////////////////////////////////////////////////


// /// Fitness class
// class FITNESS{
//   private:
//     int mode;
//     // int oldindex;
//     // int newindex;
//   public:
//     FITNESS(int m=1){
//       mode=m;
//     };
//     // fitness functions
//     double w(const double &s , const int &x, const double e=1){
//       switch(mode){
//         case 1:  // multiplicative
//           return 1 + (s * x) ;
//           break;
//         case 2:  // additive
//           return  (s * x) ;
//           break;
//         case 3:  // inverse multiplicative
//           return pow((1 + s),x);
//           break;
//         default:  // multiplicative
//           return 1 + (s * x) ;
//           break;
//       }
//     }
//     // void operator +*=(double w, double s,int x,int mode)  // operator not needed unless class
//     void wupdate(double &prevw, const double &s,const int &x) {
//       switch(mode){
//         case 2:  // additive
//           prevw= prevw + w(s,x);
//           break;
//         default:  // additive
//           prevw= prevw * w(s,x);
//           break;
//       }
//     }
// };
//
//
// /// Expectation functions
// // [[Rcpp::export]]
// arma::colvec Ey_go(const arma::Mat<double> & X, // careful arma::mat default is double
//                    const arma::colvec & s,
//                    const int & mode,
//                    double epi=1,
//                    double ref=1
//                    ){
//   int i,j;
//
//   FITNESS fit(mode);  // Initialize class and set fitness model
//
//   if(w.size()==1 ){ // Initialize vector of distribution means if it has not been provided
//     arma::colvec w(X.n_rows);
//     w.fill(ref); // IMPORTANT
//     for (i = 0; i < X.n_cols; i ++) {
//         for( j=0; j < X.n_rows ; j++){
//           fit.wupdate(w(j),s(i),X(j,i)); // works because these expressions generate a reference
//         }
//     }
//   }else{
//         for( j=0; j < X.n_rows ; j++){
//           fit.wupdate(w(j),s(i),X(j,i));
//         }
//   }
//
//   //checks
//   for( j=0; j<X.n_rows; j++) if(w(j)<0) w(j)=MIN_NUM; // **WARNING** this is necessary for non-NaN likelihood
//
//   if(epi!=1) w=pow(w,epi); // probably not very efficient
//   return(w);
// }

///// Find which selection coefficient is different in vector
// [[Rcpp::export]]
int whichS(arma::vec s1, arma::vec s2){
  int index=-9;
  int count=0;
  for(int i=0;i<s1.size();i++){
    if(s1(i) != s2(i)){ index=i; count++;}
  }
  if(count>1){
    cout<<"Counter of changes = "<<count<<endl;
    stop("More than one selection coefficient updated at a time!");
  }
  return(index);
}


/// Fitness class
class FITNESSW{
  private:
    int mode;
    int index;
    // arma::vec w = arma::randu<arma::vec>(1);
    arma::vec w = arma::vec(1).fill(1);
    arma::vec s;
    arma::Mat<double> X;
    double epi;
    double ref=1;
  public:
    FITNESSW(){
    };
    FITNESSW(int mode_,arma::Mat<double> X_,arma::vec s_,double epi_){
      mode=mode_;X=X_,s=s_,epi=epi_;
    };
    // fitness mode
    double wmode(const double &s , const int &x, const double e=1){
      switch(mode){
        case 1:  // multiplicative
          return 1 + (s * x) ;
          break;
        case 2:  // additive
          return  (s * x) ;
          break;
        case 3:  // inverse multiplicative
          return pow((1 + s),x);
          break;
        default:  // multiplicative
          return 1 + (s * x) ;
          break;
      }
    }
    void setw(arma::vec w_){w=w_;}
    // add fitness contribution of one SNP to overall fitness
    void wupdate(double &prevw, const double &s,const int &x) {     // void operator +*=(double w, double s,int x,int mode)  // operator not needed unless class
      switch(mode){
        case 2:  // additive
          prevw= prevw + wmode(s,x);
          break;
        default:  // non-additive
          prevw= prevw * wmode(s,x);
          break;
      }
    }
    // remove fitness contribution of one SNP to overall fitness (necessary to avoid repeating computations)
    void unwupdate(double &prevw, const double &s,const int &x) {
      switch(mode){
        case 2:  // additive
          prevw= prevw - wmode(s,x);
          break;
        default:  // non-additive
          prevw= prevw / wmode(s,x);
          break;
      }
    }
    // compute vector of fitness or update a previous one
    arma::vec run(
             const arma::colvec s2 ,
             double epi
            ){
              int i,j;
              if(w.size()==1 ){ // compute from scratch
                  cout << "start computing from scratch" <<endl;
                  arma::vec wres(X.n_rows);
                  // w= new arma::vec(X.n_rows)
                  wres.fill(ref);
                  cout << "start loops" <<endl;
                  for (i = 0; i < X.n_cols; i ++) {
                      for( j=0; j < X.n_rows ; j++){
                        wupdate(wres(j),s(i),X(j,i));
                      }
                  }
              w=wres;

              }else{ // only update
                int indecito=whichS(s,s2);
                cout << "indecito "<< indecito << endl;
                if(indecito != -9){ // only update if s and s2 are not identical
                  cout << "only update fitness due to position"<< indecito << endl;
                  for( j=0; j < X.n_rows ; j++){
                    unwupdate(w(j),s(indecito),X(j,indecito));
                   }
                  for( j=0; j < X.n_rows ; j++){
                    wupdate(w(j),s2(indecito),X(j,indecito));
                   }
                }
                }
              // // update with global epistatic parameter
              // if(epi!=1) w=pow(w,epi);
              // //checks
              // for( j=0; j<X.n_rows; j++) if(w(j)<0) w(j)=MIN_NUM; // **WARNING** this is necessary for non-NaN likelihood
            return(w);
            }
            // return fitness vector
      arma::vec getw(){return(w);}
};

// [[Rcpp::export]]
arma::vec W_go(
               const arma::Mat<double> X,
               const arma::vec s,
               const int mode,
               double epi=1,
               bool verbose=false){

    if(verbose) cout << "Initialize fitness" <<endl;
    FITNESSW W(mode,X,s,epi);

    if(verbose) cout << "Compute fitness" <<endl;
    arma::vec w= W.run(s,epi);

    if(verbose) cout << "done." <<endl;
    return(w);
}

// [[Rcpp::export]]
arma::vec W_update(
               const arma::Mat<double> X,
               const arma::vec s,
               const arma::vec s2,
               const int mode,
               arma::vec w,
               double epi=1,
               bool verbose=false){

    if(verbose) cout << "Initialize fitness" <<endl;
    FITNESSW W(mode,X,s,epi);

    if(verbose) cout << "Run first fitness" <<endl;
    W.run(s,epi);
    // arma::vec w_= W.run(s,epi);
    // W.setw(w_);

    // arma::vec w2_= W.run(s2,epi);
    W.run(s2,epi);

    if(verbose) cout << "done." <<endl;
    // return(w2_);
    return(W.getw());
}


////////////////////////////////////////////////////////////////////////////////
///// LIKELIHOOD

///// Fitness likelihood
// [[Rcpp::export]]
double LLGaussMix(double y,double e,double v,double p){
  double LL;
  if(y==0){
    LL = p  + (1-p) *  R::pnorm(0,e,v,true,false) ;
  }else{
    LL = (1-p) * R::dnorm(y,e,v,false);
  }
return log(LL);
}

///// Utilities subset genotypes
// [[Rcpp::export]]
arma::vec hsub(const arma::vec & h){
  arma::vec hunique = unique(h);
  arma::vec hpos(h.n_elem);
  for(int i=0; i<h.n_elem;i++){
    for(int j=0; j< hunique.n_elem;j++){
      if(h(i) == hunique(j)) hpos(i) = j;
    }
  }
  return(hpos);
}


///// Likelihood class
class LIKELIHOOD{
  private:
    int mode;
    bool TEST;
    bool verbose;
    arma::vec y;
    arma::vec h;
    arma::Mat<double> X;
    arma::vec s;
    double epi;
    arma::vec w;
    FITNESSW W;
  public:
  //Constructor
     LIKELIHOOD(
           const arma::vec  y_,
           const arma::vec  h_,
           const arma::Mat<double>  X_, // careful the arma::mat by default is double
           int mode_=1,
           bool TEST_=false,
           bool verbose_=false){
          y=y_; h=h_; X=X_;
          mode=mode_;TEST=TEST_;verbose=verbose_;

          FITNESSW W(mode,X,s,epi); //setup the fitness class
          W.run(s,epi); // and pre-compute first round
          }
  void printatributes(){
    cout <<"verbose = " << verbose << endl;
    cout <<"TEST = " <<TEST << endl;
    cout <<"mode = " <<mode << endl;
  }
  // likelihood function
  double fn(const arma::vec & snew, double b,double a, double p,double mu=1,double epi=1){
    if(TEST){
      return 1.0;
    }else{
      W.run(snew,epi); // update if there is a new set of selection coefficients
      s=snew; // update selectioncoefficients

      arma::vec e= W.getw();
      arma::vec v= a+abs(W.getw()*b);
      arma::vec hs=hsub(h);
      // Sum likelihood over all genotypes
      if(verbose) cout<< "Calculating likelihood over all genotypes..."<<  endl;
      int i;
      double L=0;
      double LL;
        for(i=0; i< y.n_elem ; i ++){
          LL= LLGaussMix(y(i)/mu,e(hs(i)),v(hs(i)),p);
          if(verbose and std::isinf(LL)){
            cout << "---" << endl;
            cout << i << endl;
            cout << y(i) << " "<< e(hs(i)) << " "<< v(hs(i)) <<" "<< p << endl;
            cout << LL << endl;
          }
          L += LL;
        }
      return(L);
    }
    }

};





