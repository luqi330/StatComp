#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sample using Rcpp
//' @description A Gibbs sample using Rcpp
//' @param n Binomial distribution coefficient
//' @param a the Beta distribution with parameters shape1
//' @param b the Beta distribution with parameters shape2
//' @param M the Sample size
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' n <- 20
//' a <- 2
//' b <- 4
//' M <- 3500
//' randomnum.C <- gibbsC(n,a,b,M)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC_98(int n, int a, int b, int M){
  NumericMatrix mat(M, 2);
  double x1=0;
  double x2=0.5;
  double alpha=0;
  double beta=0;
  for(int i=0 ; i < M; i++){
    x1 = rbinom(1, n, x2)[0];
    alpha = x1 + a;
    beta = n - x1 + b;
    x2 = rbeta(1, alpha, beta)[0];
    mat(i,0) = x1;
    mat(i,1) = x2;
  }
  return(mat);
}


