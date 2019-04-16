#ifndef BAREVERTEX_10032017
#define BAREVERTEX_10032017

#include <iostream> 
#include <stdio.h>
#include <string.h>

#include "matrix.h" 


using namespace std;

class Barevertex{
	public:
		int N;
		int Lu;
		double U1;
		double U2;
		double Xi;
		matrix<double> U;
		Barevertex(int N_in, int Lu_in, matrix<double>  U_in);
		Barevertex(int N_in,int Lu_in, double U1_in, double U2_in, double Xi_in);
		double operator()(int j1, bool s1, int j2, bool s2, int j3, bool s3, int j4, bool s4);
		void save(char *filename);
};

 
Barevertex::Barevertex(int N_in, int Lu_in, matrix<double>  U_in): N(N_in), Lu(Lu_in), U1(99.), U2(99.), Xi(99.), U(U_in){
	if(U_in.dim_c!=U_in.dim_r){
		cout<<"Error Barevertex: Interaction is not a squarematrix"<<endl;
		throw 1;
	}
	if((U_in.dim_c-1)%2!=0){
		cout<<"Error Barevertex: Interaction has even number of sites"<<endl;
		throw 1;
	}
	if(N_in!=(U_in.dim_c-1)/2){
		cout<<"Error Barevertex: Interaction does not fit chain length"<<endl;
		throw 1;
	}
}
 
 
Barevertex::Barevertex(int N_in,int Lu_in, double U1_in, double U2_in, double Xi_in): N(N_in), Lu(Lu_in), U1(U1_in), U2(U2_in), Xi(Xi_in), U(2*N+1,2*N+1){
	for(int i=-N;i<=N;++i){
		for(int j=-N;j<i;j++){
			double x=((double) i)/N;
			double y=((double) j)/N;
			double z=max(abs(x),abs(y)); //Dampening
			U(i+N,j+N)=U(j+N,i+N)=(U2/abs(i-j))*exp(-abs(i-j)/Xi)*exp(-pow(z,6)/(1-pow(z,2)));
		}
	}
	for(int i=-N;i<=N;++i){
		double z=((double) i)/N;
		U(i+N,i+N)=U1*exp(-pow(z,6)/(1-pow(z,2)));
	}
	#if BAREVERTEX_RPA_INV==1 //Caveat: only for testing, to make the interaction invertible. Do not use this for data production!
		for(int i=0; i<U.dim_r; ++i){ 
		    	U(i,i) += 1e-9;
		}
		//for(int i=0; i<U.dim_r-1; ++i){ 
		//    	U(i,i+1) += 1e-1;
		//    	U(i+1,i) += 1e-1;
		//}
	#endif
} 
 
 
double Barevertex::operator()(int j1, bool s1, int j2, bool s2, int j3, bool s3, int j4, bool s4){
	double value=0.;
	if(j1==j2 && j2==j3 && j3==j4 && s1!=s2 && s3!=s4){
		if(s1==s3){
			value=U(j1,j1);
		}
		else{
			value=-U(j1,j1);
		}
	}
	if(j1!=j2 && (abs(j1-j2)<=Lu)){
		if(j1==j3 && j2==j4 && s1==s3 && s2==s4){
			value=U(j1,j2);
		}
		if(j1==j4 && j2==j3 && s1==s4 && s2==s3){
			value=-U(j1,j2);
		}
	}
	#if BAREVERTEX_RPA_INV==1 //Caveat: only for testing, to make the interaction invertible. Do not use this for data production!
		//if( j1==j3 && j2==j4 && (s1==s2 && s2==s3 & s3==s4) ){
		if( j1==j3 && j2==j4){
			value+=1e-9;
		}
	#endif
	return value;
}
 
 
void Barevertex::save(char *filename){
	matrix<double> sv(1);
	//sv(0)=N;
	//sv.save(filename,"N");
	sv(0)=Lu;
	sv.save(filename,"Lu");
	sv(0)=U1;
	sv.save(filename,"U1");
	sv(0)=U2;
	sv.save(filename,"U2");
	sv(0)=Xi;
	sv.save(filename,"Xi");
	U.save(filename,"U");
}


#endif
