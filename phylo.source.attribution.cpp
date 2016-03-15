/* lineages through time and hazard of co
 */
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]


using namespace std;
using namespace Rcpp; 

// [[Rcpp::export]]
NumericMatrix updateWCpp2(NumericMatrix W,  NumericVector psi_a
  ,  IntegerVector utips
  ,  IntegerVector vtips
)
{
	//~ NumericMatrix WW  = clone( W);
	int ut,vt,utW,vtW;
	for (int iu = 0; iu < utips.size(); iu++){
		for (int iv = 0; iv < vtips.size(); iv++){
			ut = utips(iu)-1;
			vt = vtips(iv)-1;
			W(ut,vt) = W(ut,vt) + psi_a(ut) * psi_a(vt) / 2.; 
			W(vt,ut) = W(ut,vt); 
		}
	}
	return W; 
}
