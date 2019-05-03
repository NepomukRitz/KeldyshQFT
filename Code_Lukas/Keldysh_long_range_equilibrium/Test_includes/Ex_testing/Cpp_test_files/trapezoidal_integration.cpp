#include <iostream>
#include <string.h> 
#include <time.h> 



#include "Ex_testing.h"

class test_function{
	public:
		double operator()(double x){ return x*x;} 
};

int main(int argc, char *argv[]){
	test_function f;
	double lower = -1.;
	double upper = +1.;
	int N=1000;
	double eps=1e-5;
	double X = trapezoidal_integration<double, test_function>(lower, upper, N, f, eps); 
	cout<<"X="<<X<<endl;
	return 0;
}


