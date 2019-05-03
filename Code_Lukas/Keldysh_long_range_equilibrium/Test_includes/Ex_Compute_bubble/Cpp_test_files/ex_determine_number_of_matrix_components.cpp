#include <iostream>
#include <string.h> 
#include <time.h> 

#include "Physics.h"
#include "Numerics.h"
#include "Substitution.h"
#include "Ex_Compute_bubble.h"
#include "Ex_testing.h"


int main(){
	unsigned int seed=time(NULL);
	srand(seed);
	int const mode=0;
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
	
	matrix<int> L_structure(Nff);
	init_random(L_structure,L);
	cout<<"L_structure="<<endl;
	cout<<L_structure<<endl;

	//Test:
	cout<<"ex_determine_number_of_matrix_components(L_structure)="<<ex_determine_number_of_matrix_components(L_structure)<<endl;
	
	
	return 0;
}

