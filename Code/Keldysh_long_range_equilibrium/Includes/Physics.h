#ifndef PHYSICS_10032017
#define PHYSICS_10032017

#include <iostream> 
#include <stdio.h>
#include <string>
#include <cstring>

#include "matrix.h" 
#include "Hamiltonian.h" 
#include "Norm.h"

using namespace std;

class Physics{
	public:
		int N;
		double Vg;
		double h;
		double mu;
		double T;
		syma<complex<double> > hamiltonian;
		Physics(){};
		Physics(int N_in, double Vg_in, double h_in, double mu_in, double T_in);
		Physics(int N_in, double Vg_in, double h_in, double mu_in, double T_in, double V_sg, int pot_width);
		Physics(int N_in, double Vg_in, double h_in, double mu_in, double T_in, syma<complex<double> > hamiltonian_in);
		void save(char *filename);
		void save(string filename);
		double compare(Physics phy);
};

 
 
Physics::Physics(int N_in, double Vg_in, double h_in, double mu_in, double T_in): N(N_in), Vg(Vg_in), h(h_in), mu(mu_in), T(T_in), hamiltonian(2*N+1) {
	Hamiltonian ham(N,Vg);
	hamiltonian=ham.H();
}

Physics::Physics(int N_in, double Vg_in, double h_in, double mu_in, double T_in, double V_sg, int pot_width): N(N_in), Vg(Vg_in), h(h_in), mu(mu_in), T(T_in), hamiltonian(2*N+1) {
	Hamiltonian ham(N,Vg);
	hamiltonian=ham.H_dot(V_sg, pot_width);
}

Physics::Physics(int N_in, double Vg_in, double h_in, double mu_in, double T_in, syma<complex<double> > hamiltonian_in): N(N_in), Vg(Vg_in), h(h_in), mu(mu_in), T(T_in), hamiltonian(hamiltonian_in) {}

void Physics::save(char *filename){
	matrix<double> tmp(1);
	tmp(0)=Vg;
	tmp.save(filename,"Vg");
	tmp(0)=h;
	tmp.save(filename,"h");
	tmp(0)=mu;
	tmp.save(filename,"mu");
	tmp(0)=T;
	tmp.save(filename,"T");
	hamiltonian.save(filename,"hamiltonian");
}

void Physics::save(string filename){
	char * cstr = new char [filename.length()+1];
  	std::strcpy (cstr, filename.c_str());
	save(cstr);
}

double Physics::compare(Physics phy){
	double diff=0.0;
	diff=max(diff, abs( (double)(N - phy.N)) );
	diff=max(diff, abs(Vg - phy.Vg));
	diff=max(diff, abs(h - phy.h));
	diff=max(diff, abs(mu - phy.mu));
	diff=max(diff, abs(T - phy.T));
	if(hamiltonian.dim != phy.hamiltonian.dim){
		diff=1e16;
	}
	else{
		diff=max(diff, abs( maximumsnorm(hamiltonian - phy.hamiltonian) ) );
	}
	
	return diff;
}

#endif
