#ifndef INTEGRAND_27042017
#define INTEGRAND_27042017

#include <integrate_new.h>
#include <approxxpp.h>
#include <typeinfo>


using namespace std;

/*use this to compute the integral over the integrand g(x)*f(x), where g(x) is an extern static function */

template<class T, class K> class Integrand_static_function{
 public:
 K (*gfunction) (double);
 linear_ipol_bin<T> function;
 matrix<double> select(T &M);
 Integrand_static_function(linear_ipol_bin<T> function_in, K (*gfunction_in) (double));
 T operator()(double x);
};

template <class T, class K> Integrand_static_function<T,K>::Integrand_static_function(linear_ipol_bin<T> function_in, K (*gfunction_in) (double)): function(function_in), gfunction(gfunction_in){};
template <class T, class K> T Integrand_static_function<T,K>::operator()(double x){
 return function(x)*gfunction(x); 
}

matrix<double> Integrand_static_function<matrix<complex<double> >, double >::select(matrix<complex<double> >  &M){
 matrix<double> ret(2*M.dim_r*M.dim_c);
 
 for(int i=0, z=0; i<M.dim_r; ++i){
  for(int j=0; j<M.dim_c; ++j){
   ret(2*z)=   real(M(i,j));
   ret(2*z+1)= imag(M(i,j));
  }
 }
 return ret;
}

matrix<double> Integrand_static_function<matrix<complex<double> >, complex<double> >::select(matrix<complex<double> >  &M){
 matrix<double> ret(2*M.dim_r*M.dim_c);
 
 for(int i=0, z=0; i<M.dim_r; ++i){
  for(int j=0; j<M.dim_c; ++j){
   ret(2*z)=   real(M(i,j));
   ret(2*z+1)= imag(M(i,j));
  }
 }
 return ret;
}

matrix<double> Integrand_static_function<syma<complex<double> >, double >::select(syma<complex<double> >  &M){
 matrix<double> ret(M.dim*M.dim + 2*M.dim);
 
 for(int i=0, z=0; i<M.dim; ++i){
  for(int j=0; j<=i; ++j){
   ret(2*z)=   real(M(i,j));
   ret(2*z+1)= imag(M(i,j));
  }
 }
 return ret;
}

matrix<double> Integrand_static_function<syma<complex<double> >, complex<double> >::select(syma<complex<double> >  &M){
 matrix<double> ret(M.dim*M.dim + 2*M.dim);
 
 for(int i=0, z=0; i<M.dim; ++i){
  for(int j=0; j<=i; ++j){
   ret(2*z)=   real(M(i,j));
   ret(2*z+1)= imag(M(i,j));
  }
 }
 return ret;
}
    
matrix<double> Integrand_static_function<matrix<double>, double >::select(matrix<double>  &M){
 matrix<double> ret(M.dim_r*M.dim_c);
 
 for(int i=0, z=0; i<M.dim_r; ++i){
  for(int j=0; j<M.dim_c; ++j){
   ret(z)=   M(i,j);
  }
 }
 return ret;
}

matrix<double> Integrand_static_function<syma<double>, double >::select(syma<double>  &M){
 matrix<double> ret(M.dim*M.dim/2 + M.dim);
 
 for(int i=0, z=0; i<M.dim; ++i){
  for(int j=0; j<=i; ++j){
   ret(z)=   M(i,j);
  }
 }
 return ret;
}

/*use this to compute the integral over the integrand g(x)*f(x), where g(x) is a member function */

template<class T, class G> class Integrand_member_function{
 public:
 linear_ipol_bin<T> function;
 matrix<double> select(T &M);
 G gfunction;
 Integrand_member_function(linear_ipol_bin<T> function_in);
 T operator()(double x);
};

template <class T, class G> Integrand_member_function<T,G>::Integrand_member_function(linear_ipol_bin<T> function_in): function(function_in), gfunction(){};

template <class T, class G> T Integrand_member_function<T,G>::operator()(double x){
 return function(x)*gfunction.function(x); 
}

template<class T, class G> matrix<double> Integrand_member_function<T,G>::select(T &M){
 matrix<double> ret;
   ret.resize(M.dim*(M.dim +2) );
   ret=0.0;
   for(int i=0, z=0; i<M.dim; ++i){
    for(int j=0; j<=i; ++j){
     ret(2*z)=   real(M(i,j));
     ret(2*z+1)= imag(M(i,j));
    }
   }
 return ret;
}



#endif
