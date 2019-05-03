#include <iostream>
#include <string.h> 

#include "matrix.h" 
#include "odesolverpp.h" 

class Ode_object{
 	public:
		void operator()(double x, matrix<double> &y, matrix<double> &dydx);  
};

void Ode_object::operator()(double x, matrix<double> &y, matrix<double> &dydx){
 	dydx.resize(y.dim_r, y.dim_c);
	for(int i=0; i<y.dim_r; ++i){
		for(int j=0; j<y.dim_c; ++j){
 			dydx(i,j) = x*x;
		}
	}
}
	

int main(){
 	Ode_object f;
	double x_start = 0.0;
	double x_end = 1.0;
	double tolerance_abs = 1e-4;
	double tolerance_rel = 1e-4;
	double first_step =1e-3; //Try for first step
	double step_min = 1e-6; //minimal step size
	long nok=0;
	long nbad=0;
	{ //For scalar valued integrand:
		matrix<double> y(1); 
		y(0) = 0.0;
		odeint3(y, x_start, x_end, 1e-27, tolerance_abs, tolerance_rel, first_step,step_min,nok,nbad,f);	
		cout<<"y(0)="<<endl;
		cout<<y(0)<<endl;
	}
	{ //For matrix valued integrand:
		int N_rows=5, N_cols=3;
		matrix<double> y(N_rows,N_cols); 
		y(0) = 0.0;
		odeint3(y, x_start, x_end, 1e-27, tolerance_abs, tolerance_rel, first_step,step_min,nok,nbad,f);	
		cout<<"y(0,0)="<<endl;
		cout<<y(0,0)<<endl;
	}
	return 0;
}

