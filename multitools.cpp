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

#define ARMA_64BIT_WORD 1


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
/// A few declarations
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////////////////////


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


// [[Rcpp::export]]
double vectorSum(NumericVector x) {
   return std::accumulate(x.begin(), x.end(), 0.0);
}

// [[Rcpp::export]]
NumericVector rbinomvec(int n, double size, NumericVector prob) {
    NumericVector v(n);
    for (int i=0; i<n; i++) {v[i] = R::rbinom(1, prob[i]);}
    return(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Fitness functions
////////////////////////////////////////////////////////////////////////////////


/// Fitness class
class FITNESS{
  private:
    int mode;

  public:

    FITNESS(int m=1){
      mode=m;
    };

    // fitness functions
    double w(const double &s , const int &x, const double e=2){
      switch(mode){
        case 1:  // multiplicative
          return 1 + (s * x) ;
          break;
        case 2:  // additive
          return  (s * x) ;
          break;
        case 3:  // inverse multiplicative
          // return 1 / (1 + s * x);
          return pow((1 + s),x);
          break;
        default:  // multiplicative
          return 1 - (s * x) ;
          break;
      }
    }

    // void operator +*=(double w, double s,int x,int mode)  // operator not needed unless class
    void wupdate(double &prevw, const double &s,const int &x) {
      switch(mode){
        case 2:  // additive
          prevw= prevw + w(s,x);
          break;
        default:  // additive
          prevw= prevw * w(s,x);
          break;
      }
    }
};


/// Fintess generation functions
// [[Rcpp::export]]
arma::colvec wC(const arma::Mat<double> & X, // careful arma::mat default is double
                   const arma::colvec & s,
                   const int & mode,
                   double epi=1,
                   double mu=1
                   ){
  /*
  * Gamma expectation per haplotype - MULTIPLICATIVE MODEL
  */

  // Initialize class and set fitness model
  FITNESS fit(mode);

  // Initialize vector of distribution means
  arma::colvec myprod(X.n_rows);
  // myprod.fill(1);
  myprod.fill(mu); // IMPORTANT

  int i,j;
  for (i = 0; i < X.n_cols; i ++) {
      for( j=0; j < X.n_rows ; j++){
        fit.wupdate(myprod(j),s(i),X(j,i)); // works because these expressions generate a reference
      }
  }
  for( j=0; j<X.n_rows; j++) if(myprod(j)<0) myprod(j)=MIN_NUM; // **WARNING** this is necessary for non-NaN likelihood

  if(epi!=1) myprod=pow(myprod,epi); // probably not very efficient
  return(myprod);
}

// [[Rcpp::export]]
NumericVector samplew(NumericVector Eys,
                      NumericVector No,
                      double a=0,
                      double b=0.2,
                      double d=0.5
                      ){

  NumericVector N1(No.size());

  for(int i=1; i<No.size(); i++){
    N1(i) = vectorSum( ceil(Rcpp::rnorm(No(i),Eys(i),a+Eys(i)*b)) * rbinom(No(i),1,d) );
  }

  return(N1);
}


// [[Rcpp::export]]
NumericVector controlcapacity(NumericVector Nt,
                             double Nmax=100000
                              ){

  // arma::colvec Nupdate(Nt.n_elem);

  double Ntot=sum(Nt);
  NumericVector Nnew(Nt.size());

  if(Ntot > Nmax){
    Nnew=ceil((Nt/Ntot) * Nmax);
  }else{
    Nnew=Nt;
  }

  return(Nnew);
}




////////////////////////////////////////////////////////////////////////////////
//// Individual based simulations

//[[Rcpp::export]]
Rcpp::NumericMatrix multigenpopsimC(
                                    NumericVector fitness,
                                    double No=10000,
                                    double Nmax=100000,
                                    int t=5,
                                    double a=0,
                                    double b=0.1,
                                    double d=0.5,
                                    bool verbose = false
                                    ){


    // Total number of ecotypes
    int Ngeno=fitness.size();

    // Calculate expected individuals of each (geno)type
    NumericVector individuals(Ngeno,ceil( No * (1.0/Ngeno )));

    //// This stores the proportion of the type in the population
    NumericMatrix pop(Ngeno,t);
    pop(_,0) = individuals;

    // cout << "started matrix "<< endl;

    // if(verbose != false){
    // cout << "Simulation initiated" << endl;
    // cout << "The populations is of "<< No <<" individuals from "<< Ngeno << " (geno)types" << endl;
    // cout << "It will be evolved for "<< t << " generations" << endl;
    // // cout << "this is the (geno)types vector "<< types << endl;
    // cout << "this is the fitness vector "<< fitness << endl;
    // cout << "this is the start individuals "<< individuals << endl;
    // }



    // Iterate over t generations
    // cout << "starting iterations "<< endl;

    int i, g;
    for(i=1; i<t; i++){
      pop(_,i)=samplew(fitness,pop(_,(i-1)),a,b,d);
      pop(_,i)=controlcapacity(pop(_,i),Nmax);

    }
    return(pop);
}



// //[[Rcpp::export]]
// Rcpp::NumericMatrix populationsimC(
//                                     NumericVector types,
//                                     NumericVector fitness,
//                                     double No=10000,
//                                     int t=10,
//                                     double d=0.05,
//                                     bool verbose = false
//                                     ){
//
//
//     // Total number of ecotypes
//     int Ngeno=types.size();
//
//     // Calculate expected individuals of each (geno)type
//     NumericVector individuals(Ngeno,ceil( No * (1.0/Ngeno )));
//
//
//
//     if(verbose != false){
//     cout << "Simulation initiated" << endl;
//     cout << "The populations is of "<< No <<" individuals from "<< Ngeno << " (geno)types" << endl;
//     cout << "It will be evolved for "<< t << " generations" << endl;
//     cout << "this is the (geno)types vector "<< types << endl;
//     cout << "this is the fitness vector "<< fitness << endl;
//     cout << "this is the start individuals "<< individuals << endl;
//     }
//
//     //// This stores the proportion of the type in the population
//     Rcpp::NumericMatrix inds(Ngeno,t);
//
//     inds(_,0) = individuals;
//     rownames(inds) = types;
//
//
//     // Iterate over t generations
//     int i, g;
//     for(i=1; i<t; i++){
//        for(g=0 ; g < Ngeno ; g++){
//
//
//         // The drift Binomial sampling
//           if(inds(g,i-1) >0){
//             inds(g,i) =  R::rbinom(inds(g,i-1),1-d); // binomial sampling
//             // If NA take the previous value
//             if(inds(g,i) != inds(g,i) ){ // Na values are produced when the integer limit of the machine is exceded
//               inds(g,i)=inds(g,i-1) ;
//             }
//           }
//         // Deterministic reproduction as Poisson
//         inds(g,i) = floor(inds(g,i) * R::rpois(fitness(g))) ; // deterministic increase based on fitness
//         // inds(g,i) = floor(inds(g,i-1) * fitness(g)) ; // deterministic increase based on fitness
//
//       }
//
//     }
//     return(inds);
// }



//
// ////////////////////////////////////////////////////////////////////////////////
// /// Likelihood, Probabilities, Proposals
// ////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////
//
//
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
      // newval= runif_reflect(minhere,maxhere,min,max); // bad idea, segfault
  }

  // Check if it is out of bounds
  if(newval < min || newval>max){
      newval=(max-min)/2;
  }

  return(newval);
}
//
// // // Zero point mass proposal
// // // [[Rcpp::export]]
// // double pProposal(double & p, double min=0, double max=1, double bw=0.1){
// //     if(p<min || p> max){
// //       cout << "Parameter p needs to be in between minimum and maximum!" << endl;
// //     }
// //     double minhere=p-bw;
// //     double maxhere=p+bw;
// //     double newval = runif_reflect(minhere,maxhere,min,max);
// //     return(newval);
// // }
// //
// // // [[Rcpp::export]]
// // double bProposal(double & b, double min=-0.01, double max=10, double bw=0.1){
// //     if(b<min || b> max){
// //       cout <<"Parameter b needs to be in between minimum and maximum"<< endl;
// //     }
// //     double minhere=b-bw;
// //     double maxhere=b+bw;
// //     double newval = runif_reflect(minhere,maxhere,min,max);
// //   return(newval);
// // }
// //
// // // [[Rcpp::export]]
// // double aProposal(double & a, double min=-10, double max=10, double bw=0.1){
// //     if(a<min || a> max){
// //       cout <<"Parameter a needs to be in between minimum and maximum"<< endl;
// //     }
// //     double minhere=a-bw;
// //     double maxhere=a+bw;
// //     double newval = runif_reflect(minhere,maxhere,min,max);
// //   return(newval);
// // }
// //
// // // [[Rcpp::export]]
// // double epiProposal(double & a, double min=-0.1, double max=5, double bw=0.1){
// //     if(a<min || a> max){
// //       cout <<"Parameter epi needs to be in between minimum and maximum"<< endl;
// //     }
// //
// //     double minhere=a-bw;
// //     double maxhere=a+bw;
// //     double newval = runif_reflect(minhere,maxhere,min,max);
// //   return(newval);
// // }
// //
// // double muProposal(double & a, double min=-0.1, double max=10, double bw=0.1){
// //     if(a<min || a> max){
// //       cout <<"Parameter mu needs to be in between minimum and maximum"<< endl;
// //     }
// //
// //     double minhere=a-bw;
// //     double maxhere=a+bw;
// //     double newval = runif_reflect(minhere,maxhere,min,max);
// //   return(newval);
// // }
// //
// //
// //  // [[Rcpp::export]]
// // void test_pProposal(double p=0.5){
// //   cout << pProposal(p) << endl;
// // }
// //  // [[Rcpp::export]]
// // void test_bProposal(double b=0.2){
// //   cout << bProposal  (b) << endl;
// // }
// //  // [[Rcpp::export]]
// // void test_aProposal(double a=0.1){
// //   cout << bProposal  (a) << endl;
// // }
// //
// //  // [[Rcpp::export]]
// // void test_muProposal(double a=0.2){
// //   cout << muProposal  (a) << endl;
// // }
// //  // [[Rcpp::export]]
// // void test_epiProposal(double a=0.1){
// //   cout << epiProposal(a) << endl;
// // }
//
// class GPROPOSAL{
//   private:
//     double b; double bmin; double bmax;
//     double a; double amin; double amax;
//     double p; double pmin=0; double pmax=1;
//     double mu; double mumin; double mumax;
//     double epi; double epimin; double epimax;
//     double svar; double svarmin; double svarmax;
//     double ss; double ssmin; double ssmax;
//     int nupdates=1;
//     bool verbose;
//     double bw;
//
//   public:
//     GPROPOSAL(
//             double b_=0.5,double bmin_=0,double bmax_=1,
//             double a_=0.1,double amin_=0,double amax_=1,
//             double p_=0.5,
//             double mu_=1,double mumin_=0, double mumax_=10,
//             double epi_=1,double epimin_=0, double epimax_=5,
//             double svar_=0.5,double svarmin_=0, double svarmax_=5,
//             double ss_=0.1,double ssmin_=0, double ssmax_=1,
//             double bw_=0.1,
//             bool verbose_=false
//             ){
//       b=b_;bmin=bmin_;bmax=bmax_;
//       a=a_;amin=amin_;amax=amax_;
//       p=p_;
//       mu=mu_;mumin=mumin_;mumax=mumax_;
//       epi=epi_;epimin=epimin_;epimax=epimax_;
//       svar=svar_;svarmin=svarmin_;svarmax=svarmax_;
//       ss=ss_;ssmin=ssmin_;ssmax=ssmax_;
//       bw=bw_;
//       verbose=verbose_;
//     }
//     void setupdatesnum(int ups){nupdates=ups;}
//     void setverbose(bool verbose_){verbose=verbose_;}
//     void printatributes(){
//       cout <<"bw = " << bw << endl;
//       cout <<"b = " << b << " " << bmin << " " << bmax << endl;
//       cout <<"a = " << a << " " << amin << " " << amax << endl;
//       cout <<"p = " << p << " " << pmin << " " << pmax << endl;
//       cout <<"mu = " << mu << " " << mumin << " " << mumax << endl;
//       cout <<"epi = " << epi << " " << epimin << " " << epimax << endl;
//       cout <<"svar = " << svar << " " << svarmin << " " << svarmin << endl;
//       cout <<"ss = " << ss << " " << ssmin << " " << ssmin << endl;
//     }
//     arma::vec fn(arma::vec g){
//
//       // New proposal
//       arma::vec news=g;
//
//       // Update one position
//       double minhere,maxhere;
//       double newval;
//
//       if(verbose) cout << "Loop to substitute position" << endl;
//       for(int j=0; j< nupdates; j++){
//         int randomIndex = rand() % g.size();
//         minhere=g(randomIndex)-bw;
//         maxhere=g(randomIndex)+bw;
//         switch(randomIndex){
//             if(verbose) cout << g(randomIndex) << endl;
//             case 0: newval= runif_reflect(minhere,maxhere,bmin,bmax); break;
//             case 1: newval= runif_reflect(minhere,maxhere,amin,amax); break;
//             case 2: newval= runif_reflect(minhere,maxhere,pmin,pmax); break;
//             case 3: newval= runif_reflect(minhere,maxhere,mumin,mumax); break;
//             case 4: newval= runif_reflect(minhere,maxhere,epimin,epimax); break;
//             case 5: newval= runif_reflect(minhere,maxhere,svarmin,svarmax); break;
//             case 6: newval= runif_reflect(minhere,maxhere,ssmin,ssmax); break;
//         }
//         if(verbose) cout << newval << endl;
//         news(randomIndex) = newval;
//       }
//       if(verbose) cout << "End loop" << endl;
//
//       return(news);
//   }
// };
//
//
//
// // [[Rcpp::export]]
// void test_GProposal(double b=1,double a=1,
//                     double p=1, double mu=1,
//                     double epi=1, double svar=1,double ss=0.1,
//                     int iter=3,bool verbose = true){
//
//   arma::vec g(6);
//   g(0)= b;
//   g(1)= a;
//   g(2)= p;
//   g(3)= mu;
//   g(4)= epi;
//   g(5)= svar;
//   g(6)= ss;
//
//   GPROPOSAL GProp; // mode 1 = uniform ; mode 2 = LD
//   GProp.setverbose(verbose);
//   GProp.printatributes();
//   cout << "Original " << endl;
//   cout << g << endl;
//   cout << "Testing proposals "<< endl;
//   for(int i=0; i<iter; i++){
//     g=GProp.fn(g);
//     cout <<  g<< endl;
//   }
//
// }
//           //   gproposal=GProp.fn(par_chain.col(i));
//
// ///////////////////////////////////////////////////////////////////////////////
// Selection proposal with LD
class PROPOSAL {
  private:
    double bw;
    int nupdates;
    double min;
    double max;
    int mode;
    bool verbose;
    arma::mat R; // for LD implementatioon


  public:
  PROPOSAL(
            int nupdates_,
            int mode_=1
            ){
    nupdates=nupdates_;
    mode=mode_;
    verbose=false;
    }
  PROPOSAL(
            double bw_,
            int nupdates_,
            double min_,
            double max_,
            int mode_=1,
            bool verbose_= false
            ){
    bw=bw_;nupdates=nupdates_;min=min_;max=max_;mode=mode_;verbose=verbose_;
    // initialize R2 for default cases
    arma::mat onemat(1,1);
    onemat.fill(1);
    R=onemat;
    }
  PROPOSAL(arma::mat R2,
            double bw_,
            int nupdates_,
            double min_,
            double max_,
            int mode_=1,
            bool verbose_= false
            ){
    bw=bw_;nupdates=nupdates_;min=min_;max=max_;mode=mode_;verbose=verbose_;
    // initialize R2 for default cases
    R=R2;
    }

  void printatributes(){
    cout <<"bw = " << bw << endl;
    cout <<"nupdates = " <<nupdates << endl;
    cout <<"min = " << min << endl;
    cout <<"max = " <<max << endl;
    cout <<"mode = " <<mode << endl;
    cout <<"verbose = " <<verbose << endl;
  }
 arma::vec fn(arma::vec s){
    switch(mode){
      case 1:
        return update(s);
        break;
      case 2:
        return updateLD(s);
        break;
      default:
        return update(s);
        break;
    }
  }
  arma::vec fn(arma::vec s, double svar){
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
      /*
      * Make proposal change of one or more selection coefficients
      * from a previous vector.
      * Do not allow to go further thana bandwidth of 0.1
      */
      // New proposal
      arma::colvec news=s;

      // Update one position
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
      // New proposal
      arma::colvec news= s;

      // Update one position
      double meanhere,newval;

      if(verbose) cout << "Loop to substitute position" << endl;
      for(int j=0; j< nupdates; j++){

        int randomIndex = rand() % s.size();
        meanhere=log(1+s(randomIndex));
        if(std::isinf(meanhere)){
            meanhere= (max-min)/2;
          }

        newval = Rcpp::rnorm(1,0,svar)(0);

        news(randomIndex) = exp(newval)-1;
        if(verbose) cout << newval << endl;
        if(verbose) cout << news(randomIndex) << endl;
      }
      if(verbose) cout << "End loop" << endl;

      return(news);
  }
  arma::vec updateLD(arma::vec s){
    // New proposal
      arma::colvec news=s;

      // Update one position
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
arma::vec PropoS(int nupdates,double svar = 0.5){

  arma::vec s(nupdates);
  s.fill(0);

  PROPOSAL ps(nupdates,3);

  return(ps.fn(s,svar));
}

// // [[Rcpp::export]]
// void test_ProposalsLD(
//                     arma::mat X,
//                     double min=0,
//                     double max=1,
//                     double bw=0.1,
//                     int nupdates=1,
//                     int mode=2,
//                     int iterations=1,
//                     bool verbose=true
//                     ){
//
//   arma::vec s = Rcpp::runif(X.n_cols,0,1);
//   PROPOSAL Prop(LDrelative(X),bw,nupdates,min,max,2,verbose); // mode 1 = uniform ; mode 2 = LD
//
//   cout << "Original " << endl;
//   // cout << s << endl;
//   cout << "Testing proposals under mode = "<< mode << endl;
//   for(int i=0; i<3; i++){
//     s=Prop.fn(s);
//     // cout << s << endl;
//   }
// }

//
// // [[Rcpp::export]]
// arma::vec call_Proopsals(
//                         arma::vec s,
//                         int m=10,
//                         double min=0,
//                         double max=1,
//                         double bw=0.1,
//                         int nupdates=1,
//                         int mode=1,
//                         bool verbose=true){
//
//
//   PROPOSAL Prop(bw,nupdates,min,max,mode,verbose); // mode 1 = uniform ; mode 2 = LD
//   Prop.printatributes();
//
//   return Prop.fn(s);
// }
//
// // [[Rcpp::export]]
// void test_Proposals(int m=10,
//                     double min=0,
//                     double max=1,
//                     double bw=0.1,
//                     int nupdates=1,
//                     int mode=1,
//                     int iterations=3,
//                     bool verbose=true
//                     ){
//
//   arma::vec s = Rcpp::runif(m,0,1);
//   PROPOSAL Prop(bw,nupdates,min,max,mode,verbose); // mode 1 = uniform ; mode 2 = LD
//
//   Prop.printatributes();
//
//   cout << "Testing proposals under mode = "<< mode << endl;
//   for(int i=0; i<3; i++){
//     s=Prop.fn(s);
//     cout << s << endl;
//   }
// }
//
//
// ////////////////////////////////////////////////////////////////////////////////
//
// // class PRIOR{
// //
// // public:
// //   double min;
// //   double max;
// //   double mean;
// //   double svar;
// //   double ss;
// //   int mode;
// //
// // // Constructors
// //
// //  // PRIOR(double min_=0,double max_=1,
// //  //        double mean_=0,double variance_=1,
// //  //        int mode_=1){
// //  //  min=min_; max=max_;mean=mean_; variance=variance_; mode=mode_; };
// //
// //   PRIOR(double par1=0,double par2=1, int mode_=1){
// //        mode=mode_;
// //       switch(mode){
// //         case 1: // moc mode, return 1
// //           break;
// //         case 2: // true uniform
// //           min=par1;
// //           max=par2;
// //           break;
// //         case 3: // log +1 normal
// //           mean=par1;
// //           svar=par2;
// //           break;
// //       case 4: // log +1 mixture normal with sparcity
// //           svar=par1;
// //           ss=par2;
// //           break;
// //         default:
// //           min=0;max=1;mean=0,svar=0.5;
// //           break;
// //     }
// //   }
// //
// //   void printatributes(){
// //     cout <<"min = " << min << endl;
// //     cout <<"max = " <<max << endl;
// //     cout <<"s variance = " <<svar << endl;
// //     cout <<"mode = " <<mode << endl;
// //   }
// // // Prior functions
// //   double uniform(const arma::colvec & s){
// //     double L= 0;
// //     int N=s.n_elem;
// //     for(int i=0;i<N ;i++){
// //       L+= R::dunif(s(i),min,max,true);
// //     }
// //     return L;
// //   }
// //   double loggaussian(const arma::colvec & s, double svar){
// //     // int n=s.n_elem;
// //     arma::vec x = log(1+s);
// //     // double L = -.5*n*log(2*PI) -.5*n*log(svar) -(1/(2*svar))*sum(arma::pow((x-mean),2));
// //     double L=0;
// //     for(int i = 0; i<s.n_elem; i++){
// //       L+= R::dnorm(x(i),0,svar,true);
// //     }
// //     return L;
// //   }
// // double logmixgaussian(const arma::colvec & s, double svar, double ss){
// //
// //     arma::vec x = log(1+s);
// //
// //     double L=0;
// //     for(int i = 0; i<s.n_elem; i++){
// //       // L+= R::dnorm(x(i),0,svar,true);
// //       if(x(i)==0){
// //         L += log(ss  + (1-ss) * R::dnorm(x(i),0,svar,false)) ;
// //       }else{
// //         L += (1-ss) * R::dnorm(x(i),0,svar,true);
// //       }
// //     }
// //     return L;
// //   }
// //   // Prior distributor
// //   double fn(const arma::colvec & s){
// //     switch(mode){
// //       case 1: // moc mode, return 1
// //         return 1.0; break;
// //       case 2: // true uniform
// //         return uniform(s); break;
// //       default:
// //         return 1.0; break;
// //     }
// //   }
// //   double fn(const arma::colvec & s,const double & svar,const double & ss){
// //         switch(mode){
// //       case 1: // moc mode, return 1
// //         return 1.0; break;
// //       case 2: // true uniform
// //         return uniform(s); break;
// //       case 3:
// //         return loggaussian(s,svar); break;
// //       case 4:
// //         return logmixgaussian(s,svar,ss); break;
// //       default:
// //         return 1.0; break;
// //     }
// //   }
// // };
// //
// //
// // // [[Rcpp::export]]
// // void test_Prior(int m=10,
// //                 double min=0,
// //                 double max=1,
// //                 double mean=0,
// //                 double variance=1,
// //                 double sparsity=0.1,
// //                 int mode=1
// //                     ){
// //
// //   arma::vec s;
// //
// //   if(mode==1){
// //     s = exp( Rcpp::rnorm(m,0,variance) ) - 1;
// //   }else{
// //     s = Rcpp::runif(m,0,1);
// //   }
// //   cout << s << endl;
// //
// //   cout << "Prior mode = 1" << endl;
// //   PRIOR Pri; // mode 1 = uniform moc
// //   cout << Pri.fn(s) << endl;
// //
// //   cout << "Prior mode = 2" << endl;
// //   PRIOR Pri2(min,max,2); // mode 1 = uniform moc
// //   cout << Pri2.fn(s) << endl;
// //
// //   cout << "Prior mode = 3" << endl;
// //   PRIOR Pri3(mean,variance,3); // mode 1 = uniform moc
// //   cout << Pri3.fn(s,variance,sparsity) << endl;
// //
// //   cout << "Prior mode = 4" << endl;
// //   PRIOR Pri4(mean,variance,4); // mode 1 = uniform moc
// //   cout << Pri3.fn(s,variance,sparsity) << endl;
// // }
//
//
// ////////////////////////////////////////////////////////////////////////////////
//
// // // [[Rcpp::export]]
// // arma::vec hsub(const arma::vec & h){
// //   arma::vec hunique = unique(h);
// //   arma::vec hpos(h.n_elem);
// //   for(int i=0; i<h.n_elem;i++){
// //     for(int j=0; j< hunique.n_elem;j++){
// //       if(h(i) == hunique(j)) hpos(i) = j;
// //     }
// //   }
// //   return(hpos);
// // }
// //
// // // [[Rcpp::export]]
// // double trialLL(double hs=10){
// //
// //   arma::vec e(20);
// //   e.fill(2);
// //
// //   return e(hs);
// //
// // }
// // // [[Rcpp::export]]
// // double LLGaussMix(double y,double e,double v,double p){
// //   double LL;
// //   if(y==0){
// //     LL = p  + (1-p) *  R::pnorm(0,e,v,true,false) ;
// //   }else{
// //     LL = (1-p) * R::dnorm(y,e,v,false);
// //   }
// // return log(LL);
// // }
//
// ////////////////////////////////////////////////////////////////////////////////
// /// Likelihood class
//
// // class LIKELIHOOD{
// //   private:
// //     int mode;
// //     bool TEST;
// //     bool verbose;
// //     arma::vec y;
// //     arma::vec h;
// //     arma::Mat<double> X;
// //
// //
// //   public:
// //   //Constructor
// //      LIKELIHOOD(
// //            const arma::vec  y_,
// //            const arma::vec  h_,
// //            const arma::Mat<double>  X_, // careful the arma::mat by default is double
// //            int mode_=1,
// //            bool TEST_=false,
// //            bool verbose_=false){
// //
// //           y=y_; h=h_; X=X_;
// //           mode=mode_;TEST=TEST_;verbose=verbose_;
// //      }
// //
// //   void printatributes(){
// //     cout <<"verbose = " << verbose << endl;
// //     cout <<"TEST = " <<TEST << endl;
// //     cout <<"mode = " <<mode << endl;
// //   }
// //
// //   // calling function
// //     double fn(const arma::vec & s, double b,double a, double p,double mu=1,double epi=1){
// //       if(TEST) return 1.0;
// //       else return LLikfn(s,b,a,p,mu,epi);
// //     }
// //
// //   // likelihood function
// //   double LLikfn(const arma::vec & s, double b,double a, double p,double mu=1,double epi=1){
// //       /*
// //       * Summed log likelihood of all genotypes following each a Gammma distribution
// //       * inferred from sampling variance, mean observed genotype and selection a
// //       * set of selection coefficients
// //       */
// //
// //
// //       // Precompute all expectations of mean fitness values given genotypes X and s.
// //       if(verbose) cout<< "Precompute expectations..."<<  endl;
// //       arma::vec e= Ey_go(X,s,mode,epi,mu);
// //       // cout << e<< endl; // for debugging
// //       arma::vec v= a+abs(e*b);
// //       // cout << v<< endl; // for debugging
// //
// //       // Utilities
// //       arma::vec hs=hsub(h);
// //
// //       // Sum likelihood over all genotypes
// //       if(verbose) cout<< "Calculating likelihood over all genotypes..."<<  endl;
// //       int i;
// //       double L=0;
// //       double LL;
// //         for(i=0; i< y.n_elem ; i ++){
// //           LL= LLGaussMix(y(i),e(hs(i)),v(hs(i)),p);
// //           if(verbose and std::isinf(LL)){
// //             cout << "---" << endl;
// //             cout << i << endl;
// //             cout << y(i) << " "<< e(hs(i)) << " "<< v(hs(i)) <<" "<< p << endl;
// //             cout << LL << endl;
// //           }
// //           L += LL;
// //         }
// //
// //
// //       return(L);
// //     }
// // };
// //
// //
// // // [[Rcpp::export]]
// // void test_Likelihood(
// //                 SEXP A,
// //                 arma::vec y,
// //                 arma::vec h,
// //                 arma::vec s,
// //                 double b,
// //                 double a,
// //                 double p,
// //                 arma::uvec m,
// //                 arma::uvec n,
// //                 int Fitnessmode=1,
// //                 bool TEST=false,
// //                 bool verbose=true
// //                     ){
// //
// //   arma::Mat<double> X=BMsubset(A,n,m);
// //   cout << "Selection coefficients" << endl;
// //   cout << s << endl;
// //
// //   cout << "Likelihood" << endl;
// //   LIKELIHOOD LL(y,h,X,Fitnessmode,TEST,verbose);
// //   cout << LL.fn(s,b,a,p) << endl;
// // }
// //
// //
// //
// // // [[Rcpp::export]]
// // void test_Likelihoodall(
// //                 SEXP A,
// //                 arma::vec y,
// //                 arma::vec h,
// //                 arma::vec s,
// //                 double b,
// //                 double a,
// //                 double p,
// //                 arma::uvec m,
// //                 arma::uvec n,
// //                 int mode=1,
// //                 bool verbose=true
// //                     ){
// //
// //   arma::Mat<double> X=BMsubset(A,n,m);
// //   cout << "Selection coefficients" << endl;
// //   cout << s << endl;
// //
// //   cout << "Likelihood mode = 1 | TEST = false" << endl;
// //   LIKELIHOOD LL1(y,h,X,mode,false,verbose);
// //   cout << LL1.fn(s,b,a,p) << endl;
// //
// //   cout << "Likelihood mode = 1 | TEST = true" << endl;
// //   LIKELIHOOD LL2(y,h,X,mode,true,verbose);
// //   cout << LL2.fn(s,b,a,p) << endl;
// //
// //   cout << "Likelihood mode = 2 " << endl;
// //   LIKELIHOOD LL3(y,h,X,2,false,verbose);
// //   cout << LL3.fn(s,b,a,p) << endl;
// //
// //   cout << "Likelihood mode = 3 " << endl;
// //   LIKELIHOOD LL4(y,h,X,3,false,verbose);
// //   cout << LL4.fn(s,b,a,p) << endl;
// //
// // }
//
//
// ////////////////////////////////////////////////////////////////////////////////
// /// MCMC
// ////////////////////////////////////////////////////////////////////////////////
//
//
// // // [[Rcpp::export]]
// // List gwsMCMC(
// //             const arma::vec & y,
// //             const arma::vec & h,
// //             SEXP A, // instead of arma::mat X,
// //             const arma::colvec & s,
// //             const arma::uvec m, // the positions of SNPs
// //             const arma::uvec n , // the positions of individuals
// //             double b=0.5, double bmin=0, double bmax=1.0, // the mean variance transformation
// //             double a=0.1, double amin=0.0,   double amax=1, // the intercept of variance
// //             double p=0.5, // the proportion of zero values
// //             double mu=1.0, double mumin=0,double mumax=10,
// //             double epi=1.0, double epimin=1.0,double epimax=1.0,
// //             double svar=0.1, double svarmin=0,double svarmax=1,
// //             double ss=0.1, double ssmin=0,double ssmax=1,
// //             double bw= 0.1, // the maximum size of jumps of global parameters
// //             int nupdates=1,
// //             double min=1e-6,
// //             double max=1-1e-6 ,
// //             double iterations = 1e4,
// //             bool TEST =false ,
// //             bool verbose=false,
// //             bool debug=false,
// //             int Fitnessmode=1,
// //             int Priormode=1,
// //             int Proposalmode=1,
// //             std::string file2sink= "output.log",
// //             bool sink2file = false
// //             ){
// //
// //       if(sink2file) std::freopen(file2sink.c_str(), "w", stdout);
// //
// //       cout<< "Arguments:"<<endl;
// //       cout<< "----------"<<endl;
// //       cout<< "Range of s coefficients = ["<< min << ", " << max << "]" <<endl;
// //       cout<< "Total number of individual's observations = "<<  y.n_elem <<endl;
// //       cout<< "Total number of SNPs = "<<  s.n_elem <<endl;
// //
// //       cout<< "# iterations = "<< iterations <<endl;
// //       cout<< "TEST run = "<< TEST <<endl;
// //       cout<< "Verbose = "<< verbose <<endl;
// //       cout<< "Debug = "<< debug <<endl;
// //
// //       cout<< "----------"<<endl;
// //       cout<< "Initializing ... "<<endl;
// //       std::chrono::time_point<std::chrono::system_clock> start, end; // start chronometer values
// //       start = std::chrono::system_clock::now();
// //
// //
// //      ///////////////////////////////////////////////////////////////////////////
// //      // Deal with the matrix //
// //      ///////////////////////////////////////////////////////////////////////////
// //
// // //      arma::Mat<double> X;
// // // //      if(n.n_elem != A){
// // // //        cout<< "Reading and subsetting genome matrix ... "<<endl;
// // // //        X=BMsubset(A,n,m);
// // // //      }else{
// // // //       cout<< "Reading and subsetting columns ... "<<endl;
// // // //        X=BMsubset(A,m);
// // // //      }
// // //      cout<< "Reading and subsetting genome matrix ... "<<endl;
// // //      X=BMsubset(A,n,m);
// //
// //       arma::mat X(n.n_elem,m.n_elem);
// //       if(TYPEOF(A) == EXTPTRSXP){
// //       cout<< "Reading external pointer and subsetting genome matrix ... "<<endl;
// //          X=BMsubset(A,n,m);
// //       }else if(TYPEOF(A) == REALSXP){
// //       cout<< "Matrix provided already subsetted "<<endl;
// //         NumericMatrix Xr(A);
// //         cout << " nrow= " << Xr.nrow() << " ncol= " << Xr.ncol() << endl;
// //         arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
// //       }
// //
// //
// //      ///////////////////////////////////////////////////////////////////////////
// //      // Setup MCMC //
// //      ///////////////////////////////////////////////////////////////////////////
// //
// //       // Counter for printing
// //       int counter = 0;
// //       double percentunit = iterations/100;
// //
// //       // Initialize chain object
// //       arma::mat s_chain(s.n_elem,iterations+1);
// //       arma::mat par_chain(7, iterations+1);      /// update when adding new par
// //
// //       arma::vec  prob(iterations+1);
// //       arma::vec  paccepts(iterations+1);
// //
// //
// //       int naccepted=0;
// //       double prob0;
// //       double prob1;
// //       double Paccept;
// //       bool accept;
// //
// //       double selectionpar= s.n_elem ;
// //       double globalpar= par_chain.n_rows;
// //       double totpar= selectionpar + globalpar; // CHANGE WHEN ADDING MORE PARAMETERS
// //       if(debug) cout << "Total selection parameters = " << (selectionpar) << endl;
// //       if(debug) cout << "Total global parameters = " << (globalpar) << endl;
// //       if(debug) cout << "Total # of parameters = " << (totpar) << endl;
// //
// //
// //        ///////////////////////////////////////////////////////////////////////////
// //      // Setup MCMC conditions //
// //      ///////////////////////////////////////////////////////////////////////////
// //
// //
// //       // Proposals for variables
// //       arma::vec  sproposal=s;
// //       arma::vec gproposal(7); // Need to keep right order
// //       gproposal(0)=b;
// //       gproposal(1)=a;
// //       gproposal(2)=p;
// //       gproposal(3)=mu;
// //       gproposal(4)=epi;
// //       gproposal(5)=svar;
// //       gproposal(6)=ss;
// //       cout << "First proposal of global paramenters: "<< endl;
// //       cout <<"b = "<< gproposal(0)<<endl;
// //       cout <<"a = " <<gproposal(1)<<endl;
// //       cout <<"p = " <<gproposal(2)<<endl;
// //       cout <<"mu = " <<gproposal(3)<<endl;
// //       cout <<"epi = " <<gproposal(4)<<endl;
// //       cout <<"svar = " <<gproposal(5)<<endl;
// //       cout <<"ss = " <<gproposal(6)<<endl;
// //       Rcpp::StringVector parnames(gproposal.n_elem);
// //       parnames(0)="b";
// //       parnames(1)="a";
// //       parnames(2)="p";
// //       parnames(3)="mu";
// //       parnames(4)="epi";
// //       parnames(5)="svar";
// //       parnames(6)="ss";
// //
// //
// //       // Compute prob at first step
// //       cout<< "Start chains ..."<<endl;
// //       s_chain.col(0) = s; // set start as starting value
// //       par_chain.col(0) = gproposal;
// //
// //        // Setup probability objects
// //       LIKELIHOOD LL(y,h,X,Fitnessmode,TEST,false);
// //       PRIOR Pri(0,svar, Priormode); // loggauss
// //       // PRIOR Pri(min,max, Priormode); // uniform
// //
// //             // Proposal Prop(R,bw,nupdates,min,max,Proposalmode); // to implement R2
// //       PROPOSAL Prop(bw,
// //                     nupdates,
// //                     min,
// //                     max,
// //                     Proposalmode);
// //       GPROPOSAL GProp(b,bmin,bmax,
// //                        a, amin,amax,
// //                        p,
// //                        mu,mumin,mumax,
// //                        epi,epimin,epimax,
// //                        svar,svarmin,svarmax);
// //
// //        paccepts(0)=1;
// //       cout<< "Calculate posterior of starting point ..." ;
// //       prob(0) = Pri.fn(sproposal, gproposal(5),gproposal(6)) +
// //                 LL.fn(sproposal,
// //                  gproposal(0),
// //                  gproposal(1),
// //                  gproposal(2),
// //                  gproposal(3),
// //                  gproposal(4)
// //                           );
// //       cout << prob(0) << endl;
// //
// //    ///////////////////////////////////////////////////////////////////////////
// //     /// Handle -inf probability starts
// //     ///////////////////////////////////////////////////////////////////////////
// //
// //     if(std::isinf(prob(0)) || std::isnan(prob(0)) ){
// //       cout << "Posterior is infinite!!!. Attempt changing starting values" << endl;
// //       // stop("Posterior is infinite!!!. Attempt changing stargin values");
// //       // return List::create(Named("chain") = s_chain.t(),
// //       // Named("parchain") = par_chain.t(),
// //       // Named("parnames") = parnames,
// //       // Named("posterior") = prob,
// //       // Named("accept") = paccepts);
// //       // }else{
// //
// //
// //       PROPOSAL Propattempt(bw,
// //                       s.n_elem,
// //                       min,
// //                       max,
// //                       Proposalmode);
// //     int attemptcounter=1;
// //     while((std::isinf(prob(0)) || std::isnan(prob(0)) ) && attemptcounter < 1000 ){
// //       // cout << "attempt # " << attemptcounter << endl;
// //
// //       sproposal = Propattempt.fn(s);
// //       gproposal=GProp.fn(par_chain.col(0));
// //
// //       prob(0) =    Pri.fn(sproposal, gproposal(5),gproposal(6)) +
// //                           LL.fn(sproposal,
// //                           gproposal(0),
// //                           gproposal(1),
// //                           gproposal(2),
// //                           gproposal(3),
// //                           gproposal(4)
// //                           );
// //       attemptcounter++;
// //       }
// //      if(std::isinf(prob(0)) || std::isnan(prob(0))  ){ stop("Posterior is infinite!!!. Attempt changing stargin values"); }
// //       cout << "Successful after " << attemptcounter << " attempts " << endl;
// //       cout<< "Posterior = " ;
// //       cout << prob(0) << endl;
// //
// //     }
// //
// //
// //       ///////////////////////////////////////////////////////////////////////////
// //      // run MCMC  //
// //      ///////////////////////////////////////////////////////////////////////////
// //
// //       cout<< "Starting "<< iterations<< " MCMC iterations ..."<<endl;
// //       int counterupdate=0;
// //       int i;
// //         for(i=0; i<iterations;i++){
// //
// //           // // Propose new s and new probability
// //           if(verbose) cout<< "Generating new proposal "<<endl;
// //
// //           if(counterupdate <= selectionpar){
// //             sproposal=Prop.fn(s_chain.col(i),gproposal(5));
// //             counterupdate++;
// //           }else if(counterupdate < totpar ){
// //             gproposal=GProp.fn(par_chain.col(i));
// //             counterupdate++;
// //           }else{
// //             counterupdate=0;
// //           }
// //
// //           // Probabilities
// //           if(verbose) cout<< "Calculating next posterior "<<endl;
// //           prob0= prob(i);
// //           prob1 = Pri.fn(sproposal, gproposal(5),gproposal(6))+
// //                   LL.fn(sproposal,
// //                          gproposal(0),
// //                          gproposal(1),
// //                          gproposal(2),
// //                          gproposal(3),
// //                          gproposal(4)
// //                         );
// //           // if(debug) cout << prob1<< endl;
// //
// //             // Ratio of provabilities
// //             Paccept= exp( prob1 - prob0);
// //             accept = Rcpp::runif(1)(0)<Paccept;
// //             if(verbose) cout<< "Accept proposal " << accept <<endl;
// //
// //             if(accept){
// //               s_chain.col(i+1) = sproposal;
// //               par_chain.col(i+1) = gproposal;
// //
// //               prob(i+1) = prob1;
// //               paccepts(i+1) = Paccept;
// //               naccepted++;
// //             }else{
// //               s_chain.col(i+1) = s_chain.col(i);
// //               par_chain.col(i+1) = par_chain.col(i);
// //
// //               prob(i+1) = prob(i);
// //               paccepts(i+1) = paccepts(i);
// //             }
// //
// //             // Print row
// //             if(counter > percentunit){
// //               if(!debug & !sink2file) printProgress(i/iterations);
// //               counter =0;
// //             }else{
// //               counter++;
// //             }
// //
// //         }
// //         printProgress(1);
// //
// //       /////////
// //       // End //
// //       /////////
// //       cout<< endl<< "Summary:"<<endl;
// //       cout<< "----------"<<endl;
// //       cout<< "Acceptance ratio = "<< naccepted / iterations << endl;
// //       cout<< "Final posterior = "<< prob(iterations) << endl;
// //       cout<< "----------"<<endl;
// //       cout<<endl<< "Done."<<endl;
// //         end = std::chrono::system_clock::now();
// //         std::chrono::duration<double> elapsed_seconds1 = end-start;
// //         std::cout << "elapsed time: " << elapsed_seconds1.count() << " seconds" <<endl;
// //
// //      return List::create(Named("chain") = s_chain.t(),
// //                          Named("parchain") = par_chain.t(),
// //                          Named("parnames") = parnames,
// //                          Named("posterior") = prob,
// //                          Named("accept") = paccepts);
// //   }
