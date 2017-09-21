#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix matrix_manip(NumericVector heights, IntegerMatrix m, double k){
  NumericMatrix ml(m.nrow(),m.ncol());
  for(int i = 0; i < m.nrow(); ++i){
    for(int j = 0; j < m.ncol(); ++j){
      if(i != j && (heights[m(i,j)-1] - heights[i]) < k && (heights[m(j,i)-1] - heights[j]) < k){
	ml(i,j) = 1;
      } else {
	ml(i,j) = 0;
      }
    }
  }
  return ml;
} 
