#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Ex_Generalmatrix.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int L=1;
	int N=1;
	int Nff=4;
	int NfbP=2;
	int NfbX=2;
	int num_freq_pre=30000;
	int D=100;
	double Vg=0.25;
	double h=0; 
	double mu=-1.475;
	double T=0.0;
	Physics phy(N,Vg,h,mu,T);
	Numerics num(L,N,Nff,NfbP,NfbX,num_freq_pre,phy);
	
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;

	num.Lp_structure = L_structure;
	init_random(L_structure,L);
	num.Lx_structure = L_structure;

	cout<<"num.Lp_structure="<<endl;
	cout<<num.Lp_structure<<endl;
	cout<<"num.Lx_structure="<<endl;
	cout<<num.Lx_structure<<endl;
	cout<<"num.Nff="<<endl;
	cout<<num.Nff<<endl;

	//Test:
	
	Ex_Generalmatrix data(num);	
	data.initialize(0.3);
	
	matrix<syma<complex<double> > > sigma;
	resize_str(sigma,num.Nff,num.N);
	init(sigma,(complex<double>) 0.3);
	cout<<"abs(sigma - data.self_str(0))="<<abs(sigma - data.self_str(0))<<endl;
	cout<<"abs(sigma - data.self_str(1))="<<abs(sigma - data.self_str(1))<<endl;
	
	matrix<matrix<matrix<complex<double> > > > dyn;
	resize_str(dyn,num.Lp_structure,num.N);
	init(dyn,(complex<double>) 0.3);
	cout<<"abs(dyn - data.dynamic_str(0))="<<abs(dyn - data.dynamic_str(0))<<endl;
	cout<<"abs(dyn - data.dynamic_str(1))="<<abs(dyn - data.dynamic_str(1))<<endl;
	cout<<"abs(dyn - data.dynamic_str(2))="<<abs(dyn - data.dynamic_str(2))<<endl;
	resize_str(dyn,num.Lx_structure,num.N);
	init(dyn,(complex<double>) 0.3);
	cout<<"abs(dyn - data.dynamic_str(3))="<<abs(dyn - data.dynamic_str(3))<<endl;
	cout<<"abs(dyn - data.dynamic_str(4))="<<abs(dyn - data.dynamic_str(4))<<endl;
	cout<<"abs(dyn - data.dynamic_str(5))="<<abs(dyn - data.dynamic_str(5))<<endl;
	cout<<"abs(dyn - data.dynamic_str(6))="<<abs(dyn - data.dynamic_str(6))<<endl;
	
	data.initialize_random(D);

	Ex_Generalmatrix data2(num);	
	data2.initialize_random(D);
	Ex_Generalmatrix data3 = data + data2;
	cout<<"abs(data3.self_str - (data.self_str + data2.self_str))="<<abs(data3.self_str - (data.self_str + data2.self_str))<<endl; 
	cout<<"abs(data3.dynamic_str - (data.dynamic_str + data2.dynamic_str))="<<abs(data3.dynamic_str - (data.dynamic_str + data2.dynamic_str))<<endl; 
	cout<<"abs(data3.static_str - (data.static_str + data2.static_str))="<<abs(data3.static_str - (data.static_str + data2.static_str))<<endl; 

	double a = 0.7;	
	data3 = a*data;
	cout<<"abs(data3.self_str - a*data.self_str )="<<abs(data3.self_str - a*data.self_str)<<endl; 
	cout<<"abs(data3.dynamic_str -a*data.dynamic_str)="<<abs(data3.dynamic_str - a*data.dynamic_str)<<endl; 
	cout<<"abs(data3.static_str - a*data.static_str)="<<abs(data3.static_str - a*data.static_str)<<endl; 
	
	data.save("Ex_Generalmatrix.mat","data");
	Ex_Generalmatrix data_load;
	data_load.load("Ex_Generalmatrix.mat","data");
	cout<<"abs(data.self_str - data_load.self_str)="<<abs(data.self_str - data_load.self_str)<<endl;
	cout<<"abs(data.dynamic_str - data_load.dynamic_str)="<<abs(data.dynamic_str - data_load.dynamic_str)<<endl;
	cout<<"abs(data.static_str - data_load.static_str)="<<abs(data.static_str - data_load.static_str)<<endl;

	
	return 0;
}

