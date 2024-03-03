// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]
arma::mat columnwise_sum_cpp(arma::vec constant, arma::mat rets, arma::mat coefs, arma::mat sim_errors, int N_Sim){
  arma::mat sim_rets(N_Sim, 10);
  for (int i=0; i<10; i++){
    sim_rets.col(i) = constant + coefs(i,0) + coefs(i,1)*rets.col(0) + coefs(i,2)*rets.col(1) + coefs(i,3)*rets.col(2) + coefs(i,4)*rets.col(3) + sim_errors.col(i);
  }
  return sim_rets;
}
