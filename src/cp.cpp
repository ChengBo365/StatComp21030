#include <Rcpp.h>
using namespace Rcpp;

//' @title A markov sampler using Rcpp
//' @descriptionImplement a random walk Metropolis sampler for generating the standard Laplace distribution
//' @param n Number of random numbers
//' @param a Number of random numbers
//' @param b Number of random numbers
//' @return Returns random sequence and rejection probability
//' @examples
//' \dontrun{
//' C_Gibbs <-(25,2,2)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix C_Gibbs(int n,double a, double b){
  int N=10000;
  int x;
  double y;
  NumericMatrix M(N,2);
  M(0,0)=1;
  M(0,1)=0.5;
  for(int i=1;i<N;i++){
    y=M(i-1,1);
    M(i,0)=rbinom(1,n,y)[0];
    x=M(i,0);
    M(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  return(M) ;
}