#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Blockmatrix.h"
#include "Ex_testing.h"
#include "Ex_freq_str.h" 


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int Nff=2;
	int L=2;
	int N=1;
	int D=100;
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;

	matrix<double> wb(Nff);
	for(int i=0; i<Nff; ++i){
		wb(i) = (double) i;
	}
	matrix<matrix<matrix<complex<double> > > > A;
	resize_str(A,L_structure,N);
	init_random(A,D);
	matrix<matrix<double> > B;
	resize_str(B,L,N);
	init_random(B,D);

	cout.precision(2);
	cout.setf( std::ios::fixed, std:: ios::floatfield );
	
	//Test:
	Ex_freq_str Str(L,N,L_structure,wb,A,B);
	
	cout<<"N - Str.N="<<N-Str.N<<endl;
	cout<<"L - Str.L="<<L-Str.L<<endl;
	cout<<"abs(L_structure - Str.L_structure)="<<abs(L_structure - Str.L_structure)<<endl;
	cout<<"abs(wb - Str.wb)="<<abs(wb - Str.wb)<<endl;
	cout<<"abs(A - Str.dynamic_str)="<<abs(A - Str.dynamic_str)<<endl;
	cout<<"abs(B - Str.static_str)="<<abs(B - Str.static_str)<<endl;
	cout<<"2*N+1 - Str.Nges="<<2*N+1 - Str.Nges<<endl;
	cout<<"Nff - Str.N_freq="<<Nff - Str.N_freq<<endl;
	cout<<"Str.dynamic_str="<<endl;
	cout<<Str.dynamic_str<<endl;
	cout<<"Str.static_str="<<endl;
	cout<<Str.static_str<<endl;

	matrix<int> L_structure2(Nff);
	init_random(L_structure2,L);
	cout<<"L_structure2="<<endl;
	cout<<L_structure2<<endl;
	
	Str.L_structure = L_structure2;
	Str.resize();
	Str.initialize(0.99);

	cout<<"Str.dynamic_str="<<endl;
	cout<<Str.dynamic_str<<endl;
	cout<<"Str.static_str="<<endl;
	cout<<Str.static_str<<endl;
	
	Str.initialize_random(D);

	cout<<"Str.dynamic_str="<<endl;
	cout<<Str.dynamic_str<<endl;
	cout<<"Str.static_str="<<endl;
	cout<<Str.static_str<<endl;


	matrix<matrix<complex<double> > > tmp_ipol = Str.dyn_ipol(0.3);
	cout<<"tmp_ipol="<<endl;
	cout<<tmp_ipol<<endl;

	Str.save("Ex_freq_str.mat","Str");

	int N_load;
	int L_load;
	matrix<int> L_structure_load;
	matrix<double> wb_load;
	matrix<matrix<matrix<complex<double> > > > dyn_str_load;
	matrix<matrix<double> > stat_str_load;
	
	Ex_freq_str Str_load(L_load, N_load, L_structure_load, wb_load, dyn_str_load, stat_str_load);
	Str_load.load("Ex_freq_str.mat","Str");
	cout<<"Str_load.N - Str.N="<<Str_load.N - Str.N<<endl;
	cout<<"Str_load.L - Str.L="<<Str_load.L-Str.L<<endl;
	cout<<"abs(Str_load.L_structure - Str.L_structure)="<<abs(Str_load.L_structure - Str.L_structure)<<endl;
	cout<<"abs(Str_load.wb - Str.wb)="<<abs(Str_load.wb - Str.wb)<<endl;
	cout<<"abs(Str_load.dynamic_str - Str.dynamic_str)="<<abs(Str_load.dynamic_str - Str.dynamic_str)<<endl;
	cout<<"abs(Str_load.static_str - Str.static_str)="<<abs(Str_load.static_str - Str.static_str)<<endl;
	cout<<"Str_load.Nges - Str.Nges="<<Str_load.Nges - Str.Nges<<endl;
	cout<<"Str_load.N_freq - Str.N_freq="<<Str_load.N_freq - Str.N_freq<<endl;


	return 0;
}

