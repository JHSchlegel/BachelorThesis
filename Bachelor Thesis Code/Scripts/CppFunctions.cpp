#include <Rcpp.h>
using namespace Rcpp;


//'@param constant Numeric Vector with N_boot replicates of riskfree rate
//'@param rets Numeric Matrix with simulated Factor returns
//'@param coefs Numeric Matrix with OLS coefficients
//'@param sim_errors Numeric Matrix with bootstrapped errors
//'@param N_Sim Number of Simulations done
// [[Rcpp::export]]
NumericMatrix columnwise_sum_cpp(NumericVector constant, NumericMatrix rets, NumericMatrix coefs, NumericMatrix sim_errors, int N_Sim){
  NumericMatrix sim_rets(N_Sim, 10);
  int n = 10;
  NumericVector unit(N_Sim, 1.0);
  for (int i=0; i<n; i++){
    sim_rets(_,i) = constant + coefs(i,0)*unit+ coefs(i,1)*rets(_,0)+coefs(i,2)*rets(_,1)+coefs(i,3)*rets(_,2)+coefs(i,4)*rets(_,3)+sim_errors(_,i);
  }
  return sim_rets;
}

