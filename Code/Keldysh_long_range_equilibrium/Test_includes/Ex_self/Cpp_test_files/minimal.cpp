#include <iostream>
#include <string.h> 
#include <time.h> 

#include "matrix.h"


using namespace std;

int main(){
	#ifdef DEBUG
		cout<<"bla"<<endl;
	#endif
	syma<double> A(2); 
	A = 0.0;
	A(0,1) = 2;
	cout<<"A="<<endl;
	cout<<A(0,1)<<endl;
	return 0;
}

