// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp ;

// [[Rcpp::export]]
arma::uvec csample(IntegerVector x,
                      int size,
                      bool replace, 
                      NumericVector prob = NumericVector::create()
) {
  arma::uvec ret = as<arma::uvec>(RcppArmadillo::sample(x, size, replace, prob));
  return ret ;
}

// [[Rcpp::export]]
arma::vec permhelper(arma::mat m, int select, int nperm){
  
  // size of input matrix
  int n = int(m.n_rows);
  
  // initialize vector of permutation statistics (drawn under H0)
  arma::vec ret = arma::zeros(nperm);
  
  // loop over permutations
  for (int jj = 1; jj < (nperm + 1); jj++){
  
    // sample permutation indexes
    arma::vec auxu = arma::randu<arma::vec>(n);
    arma::uvec perms = find(auxu > 0.5);
    // switch signs for drawn indexes
    m.rows(perms) = -m.rows(perms);
    // take column means of permuted matrix
    arma::vec out = trans(mean(m, 0));
    // pick out positive column means
    arma::uvec inds = find(out > 0);
    arma::vec out2 = out.rows(inds);
    // depending on functional: take either sum of squares or sum
    arma::vec out3 = arma::vec();
    if (select == 1){
      out3 = square(out2);
    } else if (select == 2){
      out3 = out2;
    }
    // enter result
    ret(jj - 1) = sum(out3); 
    
  }
    
  return ret;
}