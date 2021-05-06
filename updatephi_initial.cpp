// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]

List updatephi2_initial(NumericMatrix proposal_phi,
                NumericMatrix ST_cur_phi,
                NumericVector cur_beta,
                NumericVector cur_tau2,
                NumericVector cur_theta,
                double rho,
                int n,int tp, NumericMatrix x, NumericMatrix offst,
                NumericMatrix y,
                List tripletf,
                NumericVector neighbour_size
                  
){
  double lik1;  // prior phi
  double lik2;  // data likelihood
  double logaccept, accepted=0;
  double sumphi;
  double acceptance; 
  double meanbelow,linearterm;
  double meanvalue,varvalue;
  for (int t=0;t<tp; t++){
    for (int s=0;s<n;s++){
      IntegerVector neighbourphi=tripletf[s]; //int len=neighbourphi.size();
      sumphi=0;
      for (int h=0; h<neighbour_size[s];h++){
        sumphi=sumphi+ST_cur_phi(neighbourphi[h]-1,t);
      };
      
      meanbelow=rho*neighbour_size[s]+1-rho;
      meanvalue=rho*sumphi/meanbelow;
      varvalue=(-2)*cur_tau2[t]/meanbelow; 
      linearterm=offst(s,t)+sum(x(s,_)*cur_beta);
      
      lik1=pow(proposal_phi(s,t)-meanvalue,2)/varvalue- (pow(ST_cur_phi(s,t)-meanvalue,2)/varvalue);
      
      lik2=(linearterm+cur_theta[t]+proposal_phi(s,t))*y(s,t)-exp(linearterm+cur_theta[t]+proposal_phi(s,t))- 
        ((linearterm+cur_theta[t]+ST_cur_phi(s,t))*y(s,t)-exp(linearterm+cur_theta[t]+ST_cur_phi(s,t)));
      
      logaccept = lik1 + lik2;
      acceptance = exp(logaccept);
      if(runif(1)[0] <= acceptance) 
      { ST_cur_phi(s,t)=proposal_phi(s,t);
        accepted=accepted+1;
      }
      else
      { 
      }
    }  
  }
  List out(2);
  out[0]=ST_cur_phi;
  out[1]=accepted;
  return out;
}
