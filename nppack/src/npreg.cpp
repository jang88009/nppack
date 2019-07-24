#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//







// [[Rcpp::export]]
double box_kernel(double x){
  double res=0;
  if(abs(x)<=1) res=1.0/2.0;
  return res;
}

// [[Rcpp::export]]
double normal_kernel(double x){
  double res;
  res=R::dnorm(x, 0, 1, FALSE);
  return res;
}

// [[Rcpp::export]]
List imput_np(arma::vec& x_vec, arma::vec& y_vec, arma::vec& d_vec, double h, int type) {
  double n = y_vec.size();
  arma::vec m_vec(n); m_vec.zeros();
  arma::vec s2_vec(n); s2_vec.zeros();
  arma::vec p_vec(n); p_vec.zeros();
  arma::vec u_vec(n); u_vec.zeros();
  arma::vec g_vec(n); g_vec.zeros();
  
  double (*my_kernel)(double);
  switch(type) {
  case 1: my_kernel=box_kernel; break;
  case 0: my_kernel=normal_kernel; break;
  default: my_kernel=normal_kernel; break;
  
  }
  
  double nu1;
  double dn1;
  double nu2;
  double dn3;
  double dn4;
  double temp;
  double temp2;
  double theta_hat=0;
  double theta_tilde=0;
  double s2_hat=0;
  for(int i=0; i<n; ++i){
    nu1=0; dn1=0; nu2=0;dn3=0; temp2=0;
    for(int j=0; j<n; ++j){
      temp=my_kernel((x_vec(i)-x_vec(j))/h);
      dn1+=temp*d_vec(j);
      nu1+=temp*d_vec(j)*y_vec(j);
      nu2+=temp*d_vec(j)*pow(y_vec(j), 2.0);
      dn3+=temp;
      dn4=0;
      for(int l=0; l<n; ++l){
        dn4+=my_kernel((x_vec(l)-x_vec(j))/h)*d_vec(l);
      }
      temp2+=temp/dn4;
    }
    m_vec(i)=nu1/dn1;
    theta_hat+=1/n*(d_vec(i)*y_vec(i)+(1-d_vec(i))*m_vec(i));
    theta_tilde+=1/n*m_vec(i);
    s2_vec(i)=nu2/dn1-pow(m_vec(i), 2.0);
    p_vec(i)=dn1/dn3;
    s2_hat+=1/n*(s2_vec(i)/p_vec(i)+pow(m_vec(i),2));
    g_vec(i)=temp2;
    u_vec(i)=m_vec(i)+d_vec(i)*(y_vec(i)-m_vec(i))*g_vec(i);
    
  }
  s2_hat=s2_hat-pow(theta_tilde, 2);
  double s2_tilde=var(u_vec);
  
  List L=List::create(theta_hat, theta_tilde, s2_hat, s2_tilde);
  return L;
}

// [[Rcpp::export]]
arma::vec np_cv(arma::vec& trax_vec, arma::vec& tray_vec, arma::vec& trad_vec, arma::vec& tesx_vec, arma::vec& tesy_vec, arma::vec& tesd_vec, arma::vec& h_vec, int type){
  double ntra = tray_vec.size();
  double ntes = tesy_vec.size();
  double mtes = sum(tesd_vec);
  double m = h_vec.size();
  arma::vec m_vec(ntes);
  
  double (*my_kernel)(double);
  switch(type) {
  case 1: my_kernel=box_kernel; break;
  case 0: my_kernel=normal_kernel; break;
  default: my_kernel=normal_kernel; break;
  
  }
  
  double nu1;
  double dn1;
  double h;
  double temp;
  double cv_err;
  double m_i;
  arma::vec cv_err_vec(m); cv_err_vec.zeros();
  for(int k=0; k<m; ++k){
    h=h_vec(k);
    cv_err=0;
    for(int i=0; i<ntes; ++i){
      nu1=0; dn1=0;
      for(int j=0; j<ntra; ++j){
        temp=my_kernel((tesx_vec(i)-trax_vec(j))/h);
        dn1+=temp*trad_vec(j);
        nu1+=temp*trad_vec(j)*tray_vec(j);
      }
      m_i=nu1/dn1;
      cv_err+=1/mtes*pow(tesy_vec(i)-m_i, 2)*tesd_vec(i);
    }
    cv_err_vec(k)=cv_err;
  }
  
  return cv_err_vec;
}

// [[Rcpp::export]]
List imput_npt(arma::vec& x_vec, arma::vec& y_vec, arma::vec& d_vec, double h, int type, double c=0.2) {
  double n = y_vec.size();
  arma::vec m_vec(n); m_vec.zeros();
  arma::vec s2_vec(n); s2_vec.zeros();
  arma::vec p_vec(n); p_vec.zeros();
  arma::vec u_vec(n); u_vec.zeros();
  arma::vec g_vec(n); g_vec.zeros();
  double (*my_kernel)(double);
  switch(type) {
  case 1: my_kernel=box_kernel; break;
  case 0: my_kernel=normal_kernel; break;
  default: my_kernel=normal_kernel; break;
  
  }
  
  double nu1;
  double dn1;
  double nu2;
  double dn3;
  double dn4;
  double temp;
  double temp2;
  double theta_hat=0;
  double theta_tilde=0;
  double s2_hat=0;
  
  for(int i=0; i<n; ++i){
    nu1=0; dn1=0; nu2=0;dn3=0; temp2=0;
    for(int j=0; j<n; ++j){
      temp=my_kernel((x_vec(i)-x_vec(j))/h);
      dn1+=temp*d_vec(j);
      nu1+=temp*d_vec(j)*y_vec(j);
      nu2+=temp*d_vec(j)*pow(y_vec(j), 2.0);
      dn3+=temp;
      dn4=0;
      for(int l=0; l<n; ++l){
        dn4+=my_kernel((x_vec(l)-x_vec(j))/h)*d_vec(l);
      }
      temp2+=temp/dn4;
    }
    if(dn1>=c/h*log(n)) {m_vec(i)=nu1/dn1; s2_vec(i)=nu2/dn1-pow(m_vec(i), 2.0);}
    theta_hat+=1/n*(d_vec(i)*y_vec(i)+(1-d_vec(i))*m_vec(i));
    theta_tilde+=1/n*m_vec(i);
    p_vec(i)=dn1/dn3;
    s2_hat+=1/n*(s2_vec(i)/p_vec(i)+pow(m_vec(i),2));
    g_vec(i)=temp2;
    u_vec(i)=m_vec(i)+d_vec(i)*(y_vec(i)-m_vec(i))*g_vec(i);
    
  }
  s2_hat=s2_hat-pow(theta_tilde, 2);
  double s2_tilde=var(u_vec);
  
  List L=List::create(theta_hat, theta_tilde, s2_hat, s2_tilde);
  return L;
}

// [[Rcpp::export]]
arma::vec npt_cv(arma::vec& trax_vec, arma::vec& tray_vec, arma::vec& trad_vec, arma::vec& tesx_vec, arma::vec& tesy_vec, arma::vec& tesd_vec, arma::vec& h_vec, int type, double c){
  double ntra = tray_vec.size();
  double ntes = tesy_vec.size();
  double mtes = sum(tesd_vec);
  double m = h_vec.size();
  arma::vec m_vec(ntes);
  
  double (*my_kernel)(double);
  switch(type) {
  case 1: my_kernel=box_kernel; break;
  case 0: my_kernel=normal_kernel; break;
  default: my_kernel=normal_kernel; break;
  
  }
  
  double nu1;
  double dn1;
  double h;
  double temp;
  double cv_err;
  double m_i;
  arma::vec cv_err_vec(m); cv_err_vec.zeros();
  for(int k=0; k<m; ++k){
    h=h_vec(k);
    cv_err=0;
    for(int i=0; i<ntes; ++i){
      nu1=0; dn1=0; m_i=0;
      for(int j=0; j<ntra; ++j){
        temp=my_kernel((tesx_vec(i)-trax_vec(j))/h);
        dn1+=temp*trad_vec(j);
        nu1+=temp*trad_vec(j)*tray_vec(j);
      }
      if(dn1>=c/h*log(ntra)) m_i=nu1/dn1;
      cv_err+=1/mtes*pow(tesy_vec(i)-m_i, 2)*tesd_vec(i);
    }
    cv_err_vec(k)=cv_err;
  }
  
  return cv_err_vec;
}



// [[Rcpp::export]]
List pse_npt(arma::vec& x_vec, arma::vec& y_vec, arma::vec& d_vec, double h, int type, double c=0.2) {
  double n = y_vec.size();
  //  arma::vec m_vec(n); m_vec.zeros();
  //  arma::vec s2_vec(n); s2_vec.zeros();
  arma::vec p_vec(n); p_vec.zeros();
  //  arma::vec u_vec(n); u_vec.zeros();
  arma::vec g_vec(n); g_vec.zeros();
  double (*my_kernel)(double);
  switch(type) {
  case 1: my_kernel=box_kernel; break;
  case 0: my_kernel=normal_kernel; break;
  default: my_kernel=normal_kernel; break;
  
  }
  
  //  double nu1;
  double dn1;
  //  double nu2;
  double dn3;
  double dn4;
  double temp;
  double temp2;
  double theta_hat_n=0;
  double theta_hat_d=0;
  double theta_hat=0;
  double theta_tilde_n=0;
  double theta_tilde_d=0;
  double theta_tilde=0;
  //  double s2_hat=0;
  
  for(int i=0; i<n; ++i){
    /*nu1=0; nu2=0;*/ dn1=0; dn3=0; temp2=0;
    for(int j=0; j<n; ++j){
      temp=my_kernel((x_vec(i)-x_vec(j))/h);
      dn1+=temp*d_vec(j);
      // nu1+=temp*d_vec(j)*y_vec(j);
      // nu2+=temp*d_vec(j)*pow(y_vec(j), 2.0);
      dn3+=temp;
      dn4=0;
      for(int l=0; l<n; ++l){
        dn4+=my_kernel((x_vec(l)-x_vec(j))/h)*d_vec(l);
      }
      temp2+=temp/dn4;
    }
    // if(dn1>=c/h*log(n)) {m_vec(i)=nu1/dn1; s2_vec(i)=nu2/dn1-pow(m_vec(i), 2.0);}
    p_vec(i)=dn1/dn3;
    g_vec(i)=temp2;
    theta_hat_n+=(d_vec(i)*y_vec(i)/p_vec(i));
    theta_hat_d+=d_vec(i)/p_vec(i);
    theta_tilde_n+=(d_vec(i)*y_vec(i)*g_vec(i));
    theta_tilde_d+=d_vec(i)*g_vec(i);
    
    //   s2_hat+=1/n*(s2_vec(i)/p_vec(i)+pow(m_vec(i),2));
    
    //   u_vec(i)=m_vec(i)+d_vec(i)*(y_vec(i)-m_vec(i))*g_vec(i);
    
  }
  theta_hat=theta_hat_n/theta_hat_d;
  theta_tilde=theta_tilde_n/theta_tilde_d;
  //  s2_hat=s2_hat-pow(theta_tilde, 2);
  //  double s2_tilde=var(u_vec);
  
  List L=List::create(theta_hat, theta_tilde);
  return L;
}

// [[Rcpp::export]]
List gpi_fun(arma::vec& x_vec, arma::vec& d_vec, arma::vec& x_new, double h, int type) {
  double n = d_vec.size();
  double m = x_new.size();
  //  arma::vec m_vec(n); m_vec.zeros();
  //  arma::vec s2_vec(n); s2_vec.zeros();
  arma::vec p_vec(m); p_vec.zeros();
  //  arma::vec u_vec(n); u_vec.zeros();
  arma::vec g_vec(m); g_vec.zeros();
  double (*my_kernel)(double);
  switch(type) {
  case 1: my_kernel=box_kernel; break;
  case 0: my_kernel=normal_kernel; break;
  default: my_kernel=normal_kernel; break;
  
  }
  
  //  double nu1;
  double dn1;
  //  double nu2;
  double dn3;
  double dn4;
  double temp;
  double temp2;

  
  for(int i=0; i<m; ++i){
    /*nu1=0; nu2=0;*/ dn1=0; dn3=0; temp2=0;
    for(int j=0; j<n; ++j){
      temp=my_kernel((x_new(i)-x_vec(j))/h);
      dn1+=temp*d_vec(j);
      // nu1+=temp*d_vec(j)*y_vec(j);
      // nu2+=temp*d_vec(j)*pow(y_vec(j), 2.0);
      dn3+=temp;
      dn4=0;
      for(int l=0; l<n; ++l){
        dn4+=my_kernel((x_vec(l)-x_vec(j))/h)*d_vec(l);
      }
      temp2+=temp/dn4;
    }
    // if(dn1>=c/h*log(n)) {m_vec(i)=nu1/dn1; s2_vec(i)=nu2/dn1-pow(m_vec(i), 2.0);}
    p_vec(i)=dn1/dn3;
    g_vec(i)=temp2;

    
    //   s2_hat+=1/n*(s2_vec(i)/p_vec(i)+pow(m_vec(i),2));
    
    //   u_vec(i)=m_vec(i)+d_vec(i)*(y_vec(i)-m_vec(i))*g_vec(i);
    
  }

  
  List L=List::create(p_vec, g_vec);
  return L;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


