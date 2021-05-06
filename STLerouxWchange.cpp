// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List updatebeta(NumericVector proposal_beta,
                NumericVector cur_beta,
                double beta_mu, double beta_v,
                int n,int tp,NumericMatrix x,NumericMatrix offst,
                NumericMatrix y, NumericMatrix ST_cur_phi, NumericVector cur_theta
){
  double lik1=0;
  double lik2=0;
  double logaccept=0, accepted;
  double acceptance; 
  for (int i=0; i<proposal_beta.size() ;i++){
    lik1=lik1+(pow((proposal_beta[i]-beta_mu),2)/(-2*beta_v)-
      pow((cur_beta[i]-beta_mu),2)/(-2*beta_v));
  }
  
  for (int i=0; i<n; i++){
    for (int t=0; t<tp; t++){
      lik2=lik2+(
        (offst(i,t)+sum(x(i,_)*proposal_beta)+cur_theta[t]+ST_cur_phi(i,t))*y(i,t)-exp(offst(i,t)+sum(x(i,_)*proposal_beta)+cur_theta[t]+ST_cur_phi(i,t))-
          ((offst(i,t)+sum(x(i,_)*cur_beta)+cur_theta[t]+ST_cur_phi(i,t))*y(i,t)-exp(offst(i,t)+sum(x(i,_)*cur_beta)+cur_theta[t]+ST_cur_phi(i,t)))
      );
    } 
  }
  
  logaccept = lik1 + lik2;
  acceptance = exp(logaccept);
  if(runif(1)[0] <= acceptance) 
  {  cur_beta=proposal_beta;
    accepted=1;
  }
  else
  { 
    accepted=0;
  }
  List out(2);
  out[0]=cur_beta;
  out[1]=accepted;
  return out;
}
// [[Rcpp::export]]
List updatetheta1(double proposal_theta1,
                  double cur_theta1,
                  double cur_theta1next,
                  double cur_alpha,
                  NumericVector cur_beta,
                  double cur_sigma2,
                  int n,int tp,NumericMatrix x,NumericMatrix offst,
                  NumericMatrix y, NumericMatrix ST_cur_phi
){
  double lik1;
  double lik2=0;
  double logaccept=0, accepted;
  double acceptance; 
  
  lik1= (proposal_theta1-cur_alpha/(1+cur_alpha*cur_alpha)*cur_theta1next)*(proposal_theta1-cur_alpha/(1+cur_alpha*cur_alpha)*cur_theta1next)/(-2*cur_sigma2/(1+cur_alpha*cur_alpha))-
    (cur_theta1-cur_alpha/(1+cur_alpha*cur_alpha)*cur_theta1next)*(cur_theta1-cur_alpha/(1+cur_alpha*cur_alpha)*cur_theta1next)/(-2*cur_sigma2/(1+cur_alpha*cur_alpha));
     
     //pow(proposal_theta1,2)/(-2*cur_sigma2)- pow(cur_theta1,2)/(-2*cur_sigma2);
  for (int i=0; i<n; i++){
    lik2=lik2+(
      (offst(i,0)+sum(x(i,_)*cur_beta)+proposal_theta1+ST_cur_phi(i,0))*y(i,0)-exp(offst(i,0)+sum(x(i,_)*cur_beta)+proposal_theta1+ST_cur_phi(i,0))-
        ((offst(i,0)+sum(x(i,_)*cur_beta)+cur_theta1+ST_cur_phi(i,0))*y(i,0)-exp(offst(i,0)+sum(x(i,_)*cur_beta)+cur_theta1+ST_cur_phi(i,0)))
    );
  }
  
  logaccept = lik1 + lik2;
  acceptance = exp(logaccept);
  if(runif(1)[0] <= acceptance) 
  {  cur_theta1=proposal_theta1;
    accepted=1;
  }
  else
  { 
    accepted=0;
  }
  List out(2);
  out[0]=cur_theta1;
  out[1]=accepted;
  return out;
}

// [[Rcpp::export]]
List updatetheta( NumericVector proposal_theta,
                  NumericVector cur_theta,
                  double cur_alpha,
                  NumericVector cur_beta,
                  double cur_sigma2,
                  int n,int tp,NumericMatrix x,NumericMatrix offst,
                  NumericMatrix y, NumericMatrix ST_cur_phi
){
  double lik1;
  double lik2=0;
  double logaccept=0, accepted=0;
  double acceptance; 
  
  for (int t=1;t<tp-1; t++){
  
    lik1=(proposal_theta[t]-cur_alpha/(1+cur_alpha*cur_alpha)*(cur_theta[t-1]+cur_theta[t+1]))*(proposal_theta[t]-cur_alpha/(1+cur_alpha*cur_alpha)*(cur_theta[t-1]+cur_theta[t+1]))/(-2*cur_sigma2/(1+cur_alpha*cur_alpha))-
         (cur_theta[t]-cur_alpha/(1+cur_alpha*cur_alpha)*(cur_theta[t-1]+cur_theta[t+1]))*(cur_theta[t]-cur_alpha/(1+cur_alpha*cur_alpha)*(cur_theta[t-1]+cur_theta[t+1]))/(-2*cur_sigma2/(1+cur_alpha*cur_alpha));
      
   
    lik2=0;// important
    for (int i=0; i<n; i++){
      lik2=lik2+(
        (offst(i,t)+sum(x(i,_)*cur_beta)+proposal_theta[t]+ST_cur_phi(i,t))*y(i,t)-exp(offst(i,t)+sum(x(i,_)*cur_beta)+proposal_theta[t]+ST_cur_phi(i,t))-
          ((offst(i,t)+sum(x(i,_)*cur_beta)+cur_theta[t]+ST_cur_phi(i,t))*y(i,t)-exp(offst(i,t)+sum(x(i,_)*cur_beta)+cur_theta[t]+ST_cur_phi(i,t)))
      );
    }
    logaccept = lik1 + lik2;
    acceptance = exp(logaccept);
    if(runif(1)[0] <= acceptance) 
    {  cur_theta[t]=proposal_theta[t];
      accepted=accepted+1;
    }
    else
    { 
    }
  }
  // when t= last year
  lik1=(proposal_theta[tp-1]-cur_alpha*cur_theta[tp-1-1])*(proposal_theta[tp-1]-cur_alpha*cur_theta[tp-1-1])/(-2*cur_sigma2)-
    (cur_theta[tp-1]-cur_alpha*cur_theta[tp-1-1])*(cur_theta[tp-1]-cur_alpha*cur_theta[tp-1-1])/(-2*cur_sigma2);
  
  lik2=0;// important
  for (int i=0; i<n; i++){
    lik2=lik2+(
      (offst(i,tp-1)+sum(x(i,_)*cur_beta)+proposal_theta[tp-1]+ST_cur_phi(i,tp-1))*y(i,tp-1)-exp(offst(i,tp-1)+sum(x(i,_)*cur_beta)+proposal_theta[tp-1]+ST_cur_phi(i,tp-1))-
        ((offst(i,tp-1)+sum(x(i,_)*cur_beta)+cur_theta[tp-1]+ST_cur_phi(i,tp-1))*y(i,tp-1)-exp(offst(i,tp-1)+sum(x(i,_)*cur_beta)+cur_theta[tp-1]+ST_cur_phi(i,tp-1)))
    );
  }
  logaccept = lik1 + lik2;
  acceptance = exp(logaccept);
  if(runif(1)[0] <= acceptance) 
  {  cur_theta[tp-1]=proposal_theta[tp-1];
    accepted=accepted+1;
  }
  else
  { 
  }
  
  List out(2);
  out[0]=cur_theta;
  out[1]=accepted;
  return out;
}



// [[Rcpp::export]]
NumericMatrix updateW(int n,
                      NumericMatrix w,
                      NumericVector Cluster
){
  NumericMatrix wclone=clone(w);
  for (int nro=0;nro<n;nro++){
    for (int nco=0; nco<n;nco++){
      double ele=wclone(nro,nco);
      double ele2=Cluster[nro];
      double ele3=Cluster[nco];
      LogicalVector v1 =(ele2==ele3);
      LogicalVector v2 =(ele2!=ele3);
      LogicalVector v3 =(ele==1);
      if(v3[0] & v1[0] )
      {
        wclone(nro,nco)=1;
      }
      else if(v3[0] & v2[0])
      {
        wclone(nro,nco)=0;
      }
    }
  }
  return wclone;
}



// [[Rcpp::export]]
List updatealpha(NumericVector cur_theta,
                 int n,
                 int tp ){
  double meantop=0;
  double meanbottom=0;
  for (int timepoint=1;timepoint<tp;timepoint++){
    meantop= meantop+cur_theta[timepoint]*cur_theta[timepoint-1];
    meanbottom= meanbottom+cur_theta[timepoint-1]*cur_theta[timepoint-1];
  }
  
  double meanvalue=meantop/meanbottom;
  List out(2);
  out[0]=meanvalue;
  out[1]=meanbottom;
  return out;
}


// [[Rcpp::export]]
double quadraticform(NumericMatrix Q_w, 
                     NumericVector temp_cur_phi,
                     int n
                       
){
  double sum1=0;
  
  for (int j=0;j<n;j++){
    for (int i=0;i<n;i++){
      sum1= sum1+temp_cur_phi[i]*Q_w(i,j)*temp_cur_phi[j];
    }
  }
  return sum1;
}  


// [[Rcpp::export]]
NumericMatrix model_fit(int n,int tp, NumericMatrix x, NumericVector offst,
                        NumericMatrix beta_mc, NumericMatrix theta_mc,
                        NumericMatrix phi_mc
){
  int m=beta_mc.nrow();
  NumericMatrix y_fitted(m,n*tp);
  for (int it=0;it<m;it++){
    for (int i=0;i<n*tp;i++){
      y_fitted(it,i)=exp(offst[i]+sum(x(i,_)*beta_mc(it,_))+phi_mc(it,i)+theta_mc(it,i));
    } 
  }
  return y_fitted;
}


// [[Rcpp::export]]
NumericMatrix model_fit_W_change(int n,int tp, NumericMatrix x, NumericVector offst,
                        NumericMatrix beta_mc, NumericMatrix theta_mc,
                        NumericMatrix phi_mc
){
  int m=beta_mc.nrow();
  NumericMatrix y_fitted(m,n);
  for (int it=0;it<m;it++){
    for (int i=0;i<n;i++){
      y_fitted(it,i)=exp(offst[i]+sum(x(i,_)*beta_mc(it,_))+phi_mc(it,i)+theta_mc(it,0));
    } 
  }
  return y_fitted;
}


// [[Rcpp::export]]
List DIC(int n, int tp, NumericVector y, 
         NumericMatrix y_fitted, NumericVector offst,
         NumericMatrix beta_mc,NumericMatrix phi_mc,
         NumericVector beta_mean, 
         NumericVector phi_mean,NumericVector theta_mean
){
  Environment stats("package:stats");
  Function dpois = stats["dpois"];
  double term1=0;
  int m=beta_mc.nrow();
  NumericMatrix deviance_all(m,n*tp);
  NumericVector deviance_mean(n*tp);
  NumericVector y_fitted_mean(n*tp);
  
  for (int i=0;i<n*tp;i++){
    y_fitted_mean[i]=exp(offst[i]+sum(beta_mean)+phi_mean[i]+theta_mean[i]);
  }
  //y_fitted m, n*tp
  for (int i =0;i<n*tp;i++){
    deviance_all(_,i)= as<NumericVector>(dpois(y[i], y_fitted(_,i)));
  }
  
  deviance_mean =log(as<NumericVector>(dpois(y, y_fitted_mean)));
  term1=sum(deviance_mean);
  
  List out(2);
  out[0]=term1;
  out[1]=deviance_all;
  return out;
}



// [[Rcpp::export]]
List DIC_changeW(int n, NumericVector y, 
         NumericMatrix y_fitted, NumericVector offst,
         NumericMatrix beta_mc,NumericMatrix phi_mc,
         NumericVector beta_mean, 
         NumericVector phi_mean,NumericVector theta_mean
){
  Environment stats("package:stats");
  Function dpois = stats["dpois"];
  double term1=0;
  int m=beta_mc.nrow();
  NumericMatrix deviance_all(m,n*1);
  NumericVector deviance_mean(n*1);
  NumericVector y_fitted_mean(n*1);
  
  for (int i=0;i<n;i++){
    y_fitted_mean[i]=exp(offst[i]+sum(beta_mean)+phi_mean[i]+theta_mean[0]);
  }
  //y_fitted m, n
  for (int i =0;i<n;i++){
    deviance_all(_,i)= as<NumericVector>(dpois(y[i], y_fitted(_,i)));
  }
  
  deviance_mean =log(as<NumericVector>(dpois(y, y_fitted_mean))); // actually here is not quite stable, better calcualte in r using dpois(y, y_fitted_mean,log=true)
  term1=sum(deviance_mean);
  
  List out(2);
  out[0]=term1; // 
  out[1]=deviance_all;
  return out;
}



// [[Rcpp::export]]
List updatephi2(NumericMatrix proposal_phi,
                NumericMatrix ST_cur_phi,
                NumericVector cur_beta,
                NumericVector cur_tau2,
                NumericVector cur_theta,
                double rho,
                int n,int tp, NumericMatrix x, NumericMatrix offst,
                NumericMatrix y,
                List tripletf,
                NumericVector neighbour_size, int t
                  
){
  double lik1;  // prior phi
  double lik2;  // data likelihood
  double logaccept, accepted=0;
  double sumphi;
  double acceptance; 
  double meanbelow,linearterm;
  double meanvalue,varvalue;
  //for (int t=0;t<tp; t++){
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
 // }
  List out(2);
  out[0]=ST_cur_phi(_,t);
  out[1]=accepted;
  return out;
}


// [[Rcpp::export]]
double Qwdet( NumericVector WstarValue,
            double rho,
            int tp
){
  double det_Q_w; 
  det_Q_w=0.5*tp*sum(log((rho *WstarValue + (1-rho))));
  return det_Q_w;
}


