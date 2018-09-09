// arma::colvec Ey_go(const arma::Mat<double> & X, // careful arma::mat default is double
//                    const arma::colvec & s,
//                    const int & mode,
//                    arma::vec w = arma::randu<arma::vec>(1),
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

////////////////////////////////////////////////////////////////////////////////
///// LIKELIHOOD

// ///// Fitness likelihood
// // [[Rcpp::export]]
// double LLGaussMix(double y,double e,double v,double p){
//   double LL;
//   if(y==0){
//     LL = p  + (1-p) *  R::pnorm(0,e,v,true,false) ;
//   }else{
//     LL = (1-p) * R::dnorm(y,e,v,false);
//   }
// return log(LL);
// }
//
// ///// Utilities subset genotypes
// // [[Rcpp::export]]
// arma::vec hsub(const arma::vec & h){
//   arma::vec hunique = unique(h);
//   arma::vec hpos(h.n_elem);
//   for(int i=0; i<h.n_elem;i++){
//     for(int j=0; j< hunique.n_elem;j++){
//       if(h(i) == hunique(j)) hpos(i) = j;
//     }
//   }
//   return(hpos);
// }
//
//
// ///// Likelihood class
// class LIKELIHOOD{
//   private:
//     int mode;
//     bool TEST;
//     bool verbose;
//     arma::vec y;
//     arma::vec h;
//     arma::Mat<double> X;
//     arma::vec s;
//     double epi;
//     FITNESSW W();
//   public:
//   //Constructor
//      LIKELIHOOD(
//            const arma::vec  y_,
//            const arma::vec  h_,
//            const arma::Mat<double>  X_, // careful the arma::mat by default is double
//            int mode_=1,
//            bool TEST_=false,
//            bool verbose_=false){
//           y=y_; h=h_; X=X_;
//           mode=mode_;TEST=TEST_;verbose=verbose_;
//           FITNESSW W(mode,X,s,epi);
//      }
//   void printatributes(){
//     cout <<"verbose = " << verbose << endl;
//     cout <<"TEST = " <<TEST << endl;
//     cout <<"mode = " <<mode << endl;
//   }
//   // likelihood function
//   double fn(const arma::vec & s, double b,double a, double p,double mu=1,double epi=1){
//     if(TEST){
//       return 1.0;
//     }else{
//       // Precompute all expectations of mean fitness values given genotypes X and s.
//       if(verbose) cout<< "Precompute expectations..."<<  endl;
//       W.run(s,epi);
//       // arma::vec
//       // arma::vec e= Ey_go(X,s,mode,epi);
//       // arma::vec v= a+abs(e*b);
//       arma::vec v= a+abs(W.getw()*b);
//       arma::vec hs=hsub(h);
//       // Sum likelihood over all genotypes
//       if(verbose) cout<< "Calculating likelihood over all genotypes..."<<  endl;
//       int i;
//       double L=0;
//       double LL;
//         for(i=0; i< y.n_elem ; i ++){
//           LL= LLGaussMix(y(i)/mu,e(hs(i)),v(hs(i)),p);
//           if(verbose and std::isinf(LL)){
//             cout << "---" << endl;
//             cout << i << endl;
//             cout << y(i) << " "<< e(hs(i)) << " "<< v(hs(i)) <<" "<< p << endl;
//             cout << LL << endl;
//           }
//           L += LL;
//         }
//       return(L);
//     }
// };







// // [[Rcpp::export]]
// void test_Likelihood(
//                 SEXP A,
//                 arma::vec y,
//                 arma::vec h,
//                 arma::vec s,
//                 double b,
//                 double a,
//                 double p,
//                 double mu,
//                 arma::uvec m,
//                 arma::uvec n,
//                 int Fitnessmode=1,
//                 bool TEST=false,
//                 bool verbose=true
//                     ){
//
//   arma::Mat<double> X=BMsubset(A,n,m);
//   cout << "Selection coefficients" << endl;
//   cout << s << endl;
//
//   cout << "Likelihood" << endl;
//   LIKELIHOOD LL(y,h,X,Fitnessmode,TEST,verbose);
//   cout << LL.fn(s,b,a,p,mu) << endl;
// }
//
//
//
// // [[Rcpp::export]]
// void test_Likelihoodall(
//                 SEXP A,
//                 arma::vec y,
//                 arma::vec h,
//                 arma::vec s,
//                 double b,
//                 double a,
//                 double p,
//                 arma::uvec m,
//                 arma::uvec n,
//                 int mode=1,
//                 bool verbose=true
//                     ){
//
//   arma::Mat<double> X=BMsubset(A,n,m);
//   cout << "Selection coefficients" << endl;
//   cout << s << endl;
//
//   cout << "Likelihood mode = 1 | TEST = false" << endl;
//   LIKELIHOOD LL1(y,h,X,mode,false,verbose);
//   cout << LL1.fn(s,b,a,p) << endl;
//
//   cout << "Likelihood mode = 1 | TEST = true" << endl;
//   LIKELIHOOD LL2(y,h,X,mode,true,verbose);
//   cout << LL2.fn(s,b,a,p) << endl;
//
//   cout << "Likelihood mode = 2 " << endl;
//   LIKELIHOOD LL3(y,h,X,2,false,verbose);
//   cout << LL3.fn(s,b,a,p) << endl;
//
//   cout << "Likelihood mode = 3 " << endl;
//   LIKELIHOOD LL4(y,h,X,3,false,verbose);
//   cout << LL4.fn(s,b,a,p) << endl;
//
// }






































               // arma::vec w = arma::randu<arma::vec>(1),
