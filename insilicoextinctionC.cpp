#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <map>
#include <string>


//################################################################################
//## POSIBE POpulation SImulation By Leading Extinction
//################################################################################

// when armadillo is loaded, remove this below
// #include <Rcpp.h>

// when armadillo is loaded, include this below but remove #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


////////////////////////////////////////////////////////////////////////////////
// Function declarations

double sum(NumericVector x);


////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector extinCt(NumericVector pop,
                           double extfrac=0.9,
                           NumericVector extprob=NumericVector::create(),
                           bool verbose = false) {

  // get the population size
  int len = pop.size();

  // sample the survivors
  int size= floor((1-extfrac)*len);
  NumericVector survivors = RcppArmadillo::sample(pop,size,false,extprob);

  if(verbose!=false){
  cout << "Extinction simulated" << len << endl;
  cout << "The original (geno)type number was " << len << endl;
  cout << "After extinction size is " << size<< endl;
  }
  return(survivors);
}


Rcpp::NumericVector defaultfq(Rcpp::NumericVector freq, int Ngeno){
    if( freq(0) == 0){
      double eqfreq= 1.0/Ngeno;

      Rcpp::NumericVector freqdef(Ngeno); //, eqfreq
      freqdef= freqdef + eqfreq;
      freq=freqdef;
    }
    return(freq);
}

//[[Rcpp::export]]
void printfitness(NumericVector fitness){
    cout << "this is the fitness vector "<< fitness << endl;
}



////////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
Rcpp::NumericMatrix populationsimC(
                                    NumericVector types,
                                    NumericVector fitness,
                                    double No=10000,
                                    int t=20,
                                    double d=0.05,
                                    bool verbose = false
                                    ){


    // Total number of ecotypes
    int Ngeno=types.size();

    // Calculate expected individuals of each (geno)type
    NumericVector individuals(Ngeno,ceil( No * (1.0/Ngeno )));



    if(verbose != false){
    cout << "Simulation initiated" << endl;
    cout << "The populations is of "<< No <<" individuals from "<< Ngeno << " (geno)types" << endl;
    cout << "It will be evolved for "<< t << " generations" << endl;
    cout << "this is the (geno)types vector "<< types << endl;
    cout << "this is the fitness vector "<< fitness << endl;
    cout << "this is the start individuals "<< individuals << endl;
    }

    //// This stores the proportion of the type in the population
    Rcpp::NumericMatrix inds(Ngeno,t);

    inds(_,0) = individuals;
    rownames(inds) = types;

    // Iterate over t generations
    int i, g;
    for(i=1; i<t; i++){
       for(g=0 ; g < Ngeno ; g++){


        // The drift Binomial sampling
          if(inds(g,i-1) >0){
            inds(g,i) =  R::rbinom(inds(g,i-1),1-d); // binomial sampling
            // If NA take the previous value
            if(inds(g,i) != inds(g,i) ){ // Na values are produced when the integer limit of the machine is exceded
              inds(g,i)=inds(g,i-1) ;
            }
          }
        // Deterministic reproduction as Poisson
        inds(g,i) = floor(inds(g,i) * R::rpois(fitness(g))) ; // deterministic increase based on fitness
        // inds(g,i) = floor(inds(g,i-1) * fitness(g)) ; // deterministic increase based on fitness

      }

    }
    return(inds);
}


////////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
std::multimap<double, Rcpp::NumericMatrix>
    populationsimCgrid(Rcpp::NumericVector types,
                      Rcpp::NumericVector fitness,
                      double No=10000,
                      int t=20,
                      double dmin=0.800,
                      double dmax=0.999,
                      double dstep=0.001,
                      int reps=1,
                      bool verbose = false){

    if(verbose != false){
    cout << "Simulation grid initiated" << endl;
    cout << "The drift simulations will run from " << dmin << " to " << dmax << endl;
    }

    // 1 // Generate the data structure of the result object
    std::multimap<double, Rcpp::NumericMatrix> evopop;

    // 2 // Iterate over drift parameter grid
    for(float d = dmin; d <= dmax; d += dstep){

      // 2.1 // run
      Rcpp::NumericMatrix popi = populationsimC(types,
                                                fitness,
                                                No,
                                                t,
                                                d,
                                                verbose);

      // 2.2 // generate replicates and average
          for(int r=1; r<reps;r++ ){
            Rcpp::NumericMatrix popi2 = populationsimC(types,
                                                    fitness,
                                                    No,
                                                    t,
                                                    d,
                                                    verbose);

            for(int i = 0; i < popi.nrow(); ++i){
            for(int j = 0; j < popi.ncol(); ++j){
                popi(i,j) = popi(i,j) + popi2(i,j);
            }}
          }
          for(int j = 0; j < popi.ncol(); ++j){
              popi(_,j) = floor( popi(_,j) / reps ) ;
          }

      // 2.3 //  Subset simulations in which evolutionary rescue exist
        // if(sum(popi(_,0)) > sum(popi(_,1))  &
        //   sum(popi(_,0)) < sum(popi(_,t-1))  ){
        // if(sum(popi(_,t-1)) > 0 ){ // originally only evo rescue, but I lose all those that survive w/ evo rescue. Better to see limits, after easy to subset
        //   if(verbose != false){
        //         cout <<"With d="<<d << " Nt1="<< sum(popi(_,0)) << " and Nt2=" << sum(popi(_,1)) << endl;
        //         }
        //    evopop[d] =popi ;
        // }

      // 2.3 alternative //  No subset
      // evopop[d] =popi ;
      evopop.insert(pair<double, Rcpp::NumericMatrix >(d, popi));

     }// end of d loop

    return(evopop);
}



////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
std::multimap<double, Rcpp::NumericMatrix> popXsim(Rcpp::NumericVector types,
                                              Rcpp::NumericVector fitness,
                                              std::vector<int> N,
                                              std::vector<double> extfrac,
                                              int t=20,
                                              double dmin=0.800,
                                              double dmax=0.999,
                                              double dstep=0.001,
                                              int reps=1,
                                              Rcpp::NumericVector extprob = Rcpp::NumericVector::create(),
                                              bool verbose = false){

    // Generate the data structure of the result object
    std::multimap<double, Rcpp::NumericMatrix> evopop;

    // Update list after extinction of a fraction
    for(int i=0; i< extfrac.size();i++){
    for(int j=0; j< N.size();j++){

    cout << "Initial No: "<<N.at(j) << ". Extinction fraction: "<<extfrac.at(i) << endl;

      // Sample survivors
      Rcpp::NumericVector toextinct = extinCt(types,
                                              extfrac.at(i),
                                              extprob,
                                              verbose) ;

      // Update the (geno)types and other input vectors
      NumericVector survtypes ;
      NumericVector survfitness ;
        for(int i=0; i<types.size();i++){
          if(std::find(toextinct.begin(), toextinct.end(), types(i)) != toextinct.end()){ // != because want to select those not to extinct
             survtypes.push_back(types(i)) ;
             survfitness.push_back(fitness(i)) ;
          }
        }



    // Define the result data structure and run
    std::multimap<double, Rcpp::NumericMatrix> newevopop = populationsimCgrid(survtypes,
                                                                         survfitness,
                                                                         N.at(j),
                                                                         t,
                                                                         dmin,
                                                                         dmax,
                                                                         dstep,
                                                                         reps,
                                                                         verbose);
    evopop.insert(newevopop.begin(), newevopop.end()) ;

    }}

    return(evopop);
}

