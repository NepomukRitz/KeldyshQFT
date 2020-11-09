//
//  main.cpp
//  Multiloop fRG for X-ray-edge singularity, paper: multiloop functional renormalization group that sums up all parquet diagrams
//
//  Created by Fabian Kugler on 10.01.17.
//  Copyright © 2016 Fabian Kugler. All rights reserved.
//

#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
//#include<iomanip> // for ofstream precision, cout << setprecision(2);
#include <cmath>
#include <functional>
#include <omp.h>
using namespace std;

// parameters are specified in compiling process, default values are given here
// parameters are: inverse temperature, bandgap, chemical potential, interaction value, number of frequencies in parametrization (that accounts for vertex asymptotics) for four-point vertex (NUMFREQ1-3), number of frequencies in bosonic self-energies (NFREQ1) and three-point vertices (NFREQ2), number of frequencies for the precomputation of propagators, mode, regulator, multiloop order, for mode 4: inclusion of bosonic self-energies, weights in exchange and pairing channel
// physical parameters: inverse temperature, bandgap, chemical potential, interaction value
#ifndef BETA 
	#define BETA 100. 
#endif
#ifndef XID 
	#define XID -0.2 
#endif
#ifndef MU 
	#define MU 5. 
#endif
#ifndef UVAL 
	#define UVAL 0.28 
#endif 
// sizes for four-point vertices (NUMFREQ1-3) in parametrization that accoutns for vertex asymptotics (note: default values are chosen small for testing)
#ifndef NUMFREQ1
	#define NUMFREQ1 100//8000//2N^3
#endif
#ifndef NUMFREQ2
	#define NUMFREQ2 20//80//N^(3/2)
#endif
#ifndef NUMFREQ3
	#define NUMFREQ3 10//20//N
#endif
// sizes for bosonic two- and three-point vertices (NFREQ1,2)
#ifndef NFREQ1
	#define NFREQ1 1
#endif
#ifndef NFREQ2
	#define NFREQ2 1
#endif
// size for precomputed arrays (NAUX)
#ifndef NAUX
	#define NAUX 2000//30000// >= BETA*MU + NUMFREQ1
#endif

//test mode, MODE=2
// #define NUMFREQ1 20
// #define NUMFREQ2 4
// #define NUMFREQ3 3
// #define NFREQ1 20
// #define NFREQ2 5

#ifndef MODE
	#define MODE 1//0: parquet equations; 1: fRG I_{(a)p}; 2: fRG I_{(a)p} & Pi^{gamma}; 3: fRG I_{(a)p} & Pi^{gamma} with trivial (constant) 4-point vertex; 4: HS; only 0 and 1 thoroughly tested!
#endif
#ifndef REGULATOR
	#define REGULATOR 1//0: delta; 1: Litim; 2: smooth; 3: oscillating
#endif
#ifndef MULTILOOP
	#define MULTILOOP 2//1: standard fRG; >=2 multiloop fRG
#endif
// specific parameters for each regulator
#if REGULATOR==2
#ifndef REGEXP
	#define REGEXP 2.
#endif
#endif
#if REGULATOR==3
#ifndef OSCFACTOR
	#define OSCFACTOR 1.
#endif
#endif
// inclusion of bosonic self-energies in mixed fermionic-bosonic flow
#ifndef WITH_BOSONIC_SELFENERGIES
	#define WITH_BOSONIC_SELFENERGIES 1
#endif
// weights in exchange and pairing channel and mixed fermionic-bosonic flow
#ifndef UXVAL
	#define UXVAL 0.
#endif
#ifndef UYVAL
	#define UYVAL 0.
#endif

//#define NOMOVE 1//if NOMOVE is defined to 1, all move semantics are ignored

#include "myheader.h"
#include "myfuncs.h"


State derivative(const double Lambda, const State& state) {
	// return derivative of all vertex functions given current state and current Lambda
	State dState;
	
	// MODE>1 for mixed fermionic-bosonic flow	
#if MODE>1
	bubble_Pi(state, dState.PiX, Lambda, 1); // \partial_{\Lambda} \Pi^{\chi}
#endif
#if MODE==2
	bubble_Gamma3(state, dState.GammaCDX, Lambda, 3); // \partial_{\Lambda} \Gamma^{c,d,\chi}
#elif MODE==3
	bubble_Gamma3(state, dState.GammaCDX, Lambda, 5); // \partial_{\Lambda} \Gamma^{c,d,\chi}
#elif MODE==4
	bubble_Gamma3(state, dState.GammaCDX, Lambda, 2); // \partial_{\Lambda} \Gamma^{c,d,\chi}
	bubble_Pi(state, dState.PiY, Lambda, 2);  // \partial_{\Lambda} \Pi^{\psi}
	bubble_Gamma3(state, dState.GammaCDY, Lambda, 2); // \partial_{\Lambda} \Gamma^{c,d,\psi}
#endif

// vertex flow: \partial_{\Lambda} \Gamma^{(4)}
#if MODE<3
	// 1-loop fRG
	bubble(state, dState.Ip, Lambda, 1, 1); // \dot{gamma}_a \equiv \dot{I}_p: bubble in channel a (1) with single-scale propagatar and full vertices (1)
	bubble(state, dState.Iap, Lambda, 2, 1); // \dot{gamma}_p \equiv \dot{I}_ap: bubble in channel p (2)
	
	// 2-loop fRG
	if(MULTILOOP > 1) {
		// allocate memory for left and right vertex; these are nonsymmetric (ns) in the fermionic frequencies, right one is obtained from left one by symmetry
		CVertex dIpL; dIpL.K2ns = initCMat(NUMFREQ2, 2*NUMFREQ2); dIpL.K3ns = initCTen(NUMFREQ3, 2*NUMFREQ3, 2*NUMFREQ3);
		CVertex dIapL; dIapL.K2ns = initCMat(NUMFREQ2, 2*NUMFREQ2); dIapL.K3ns = initCTen(NUMFREQ3, 2*NUMFREQ3, 2*NUMFREQ3);
		
		bubble(state, dState.Iap, dIpL, Lambda, 1, 2); // bubble in a channel (1) with full propagators and differentiated vertex left (2)
		bubble(state, dState.Ip, dIapL, Lambda, 2, 2); // bubble in p channel (2)
		
		// numerically: dIpR(w, v, w_bar) = dIpL(v, w, w_bar) very well fulfilled. can be checked using below lines
// 		CVertex dIpR; dIpR.K2ns = initCMat(NUMFREQ2, 2*NUMFREQ2); dIpR.K3ns = initCTen(NUMFREQ3, 2*NUMFREQ3, 2*NUMFREQ3);
// 		CVertex dIapR; dIapR.K2ns = initCMat(NUMFREQ2, 2*NUMFREQ2); dIapR.K3ns = initCTen(NUMFREQ3, 2*NUMFREQ3, 2*NUMFREQ3);
// 		bubble(state, dState.Iap, dIpR, Lambda, 1, 3); 
// 		bubble(state, dState.Ip, dIapR, Lambda, 2, 3); 
// 		combineCVertex(dIpL, dIpR, dIpT);
// 		combineCVertex(dIapL, dIapR, dIapT);

		CVertex dIpT, dIapT; // allocate memory for total vertex, symmetric in fermionic frequencies
		
		symmetrizeCVertex(dIpL, dIpT); // get total vertex by summing left and right vertex \equiv symmetrizing left vertex
		symmetrizeCVertex(dIapL, dIapT);

		dState.Ip = dState.Ip + dIpT; // update r.h.s.
		dState.Iap = dState.Iap + dIapT;
		
		// higher-loop fRG
		if(MULTILOOP > 2) {
			CVertex dIpC, dIapC; // allocate memory for center vertex
//			CVertex dIpC2, dIapC2;
			
			for (int loop = 3; loop <= MULTILOOP; ++loop) {
				bubble(state, dIpL, dIpC, Lambda, 1, 4); // bubble in a channel (1) with full propagators and nonsymmetric vertex on the left
				bubble(state, dIapL, dIapC, Lambda, 2, 4); // bubble in p channel (2)

				// compute mean between left and right realization of center (unnecessary)
// 				bubble(state, dIpL, dIpC2, Lambda, 1, 7);
// 				bubble(state, dIapL, dIapC2, Lambda, 2, 7);
				
				bubble(state, dIapT, dIpL, Lambda, 1, 2); // see above
				bubble(state, dIpT, dIapL, Lambda, 2, 2); 
				
				//symmetrizeCVertex(dIpL, dIpT);
				//symmetrizeCVertex(dIapL, dIapT);
				symmetrizeCVertex_withC(dIpL, dIpC, dIpT); // get total vertex by summing left and right vertex (\equiv symmetrizing left vertex) and center vertex
				symmetrizeCVertex_withC(dIapL, dIapC, dIapT);
				
				dState.Ip = dState.Ip + dIpT; // update r.h.s.
				dState.Iap = dState.Iap + dIapT;
			}
		}
	}
#endif
	return dState;
}


// substitution functions in ODE solver, Lambda(l) and d Lambda(l) / dl
// simple choice:
// double LambdaL(const double l, const double Lambda_0) {
// 	return Lambda_0 * exp(-l);
// }
// double dLambdaL(const double l, const double Lambda_0) {
// 	return - Lambda_0 * exp(-l);
// }
// more complicated choice decreases extremely fast for large Lambda
double LambdaL(const double l, const double Lambda_0) {
	const static double a = 1.2;
	const double temp = exp(-l);
	return Lambda_0 * (a - 1.) * temp /( a - temp );
}
double dLambdaL(const double l, const double Lambda_0) {
	const static double a = 1.2;
	const double temp = exp(-l);
	return - Lambda_0 * (a - 1.) * a * temp / ( (a-temp)*(a-temp) );
}

// if fixed step Runge-Kutta solver is used, use most efficient RK4 algorithm by defining fastRK4:
#define fastRK4

// propagatate state by solving the flow equation 
void propagate(State& state, const double Lambda_0, const double Lambda_f_fraction, int steps, int ODE_solver, bool print, string suffix) {
	//solve ODE, i.e., evolve state from \Lambda = \Lambda_0 to \Lambda = 0. In some cases, Lambda_f < minimal Matsubara frequency sufficient.
	double Lambda_f;
#if REGULATOR==0
	Lambda_f = 0.002 * PiOverBeta;
#elif REGULATOR==1
	Lambda_f = 0.9 * PiOverBeta; // derivative vanishes for Lambda < M_PI/state.par.beta
#elif REGULATOR==2
	Lambda_f = 0.2 * PiOverBeta; // derivative vanishes for Lambda << M_PI/state.par.beta
#elif REGULATOR==3
	Lambda_f = 0.01 * PiOverBeta;
#else
	Lambda_f = Lambda_f_fraction * PiOverBeta;
#endif
	double h = 0.025; // initial step size
	if (steps > 1000) h = 0.015;
	
	if (ODE_solver==0) { // adaptive RungeKutta
		DVec c_coeff(6); DMat a_coeff(6, DVec(6)); DVec b_coeff(6); DVec d_coeff(6); cashKarp(c_coeff, a_coeff, b_coeff, d_coeff); vector<State> auxStates(6);
		State preState, refState, tempState;
		int ngood=0; int nretry=0; int nbad=0; double adaptiveError; double epsilon = 1.e-4; 
		double hmax = 1.; double hmin = 0.00001; double l = 0.;
		
		for (int m = 0; m < steps; ++m) {
			if (print and m%1==0) cout << "iteration at step " << m << " at Lambda*beta/Pi " << LambdaL(l, Lambda_0)/PiOverBeta <<  " with previous error " << adaptiveError << " \t\t\t\t\t\t\twith h = " << h <<  " \t\t and l = " << l << endl;
			if (LambdaL(l, Lambda_0) < Lambda_f) {
				cout << "steps " << ngood << " " << nretry << " " << nbad << "\t" << "Lambda/Lambda_f " << LambdaL(l, Lambda_0)/Lambda_f << endl;
				break;
			}
			
			for (int i = 0; i < 6; ++i) {
				tempState = state;
				for (int j = 0; j < i; ++j) tempState = tempState + auxStates[i-1] * a_coeff[i][j];
				auxStates[i] = derivative(LambdaL(l+h*c_coeff[i], Lambda_0), tempState) * ( dLambdaL(l+h*c_coeff[i], Lambda_0) * h );
			}
			
			preState = state; refState = state; 
			for (int i = 0; i < 6; ++i) preState = preState + auxStates[i] * b_coeff[i];
			for (int i = 0; i < 6; ++i) refState = refState + move(auxStates[i]) * d_coeff[i];
			
			adaptiveError = getError(preState, refState)/epsilon;
			
			if (adaptiveError > 1.) {
				if (h <= hmin) {
					nbad += 1;
					if (print) cout << "hmin reached at step " << m << " with error " << adaptiveError << " at Lambda*beta/Pi " << LambdaL(l, Lambda_0)/PiOverBeta <<  " with h = " << h << endl;
					//state = move(refState);
					state = move(preState);
					l += h; h = hmin;
				}
				else {
					h = h * max( 0.8*pow(adaptiveError, -0.2), 0.01 );//6
					nretry += 1; 
					if (print) cout << "\tretry at step " << m << " with error " << adaptiveError << " at Lambda*beta/Pi " << LambdaL(l, Lambda_0)/PiOverBeta <<  " with h = " << h << endl;
					continue;}
			}
			else {
				l += h;
				h = h * min( 0.6*pow(adaptiveError, -0.2), 100. );//7
				ngood += 1;

				if (h > hmax) {
					h = hmax;
				//if (print) cout << "hmax reached at step " << m << " with error " << adaptiveError << " at Lambda*beta/Pi " << Lambda/PiOverBeta <<  " with h = " << h << endl;
				}
				
				//state = move(refState);
				state = move(preState); 
			}
		}
		// last step
		if (print) preState = state;
		State deriv = derivative(LambdaL(l, Lambda_0), state);
		state = state + deriv * (0.-LambdaL(l, Lambda_0));
		if (print) cout << "change in last step " << getError(preState, state) << endl;
	}
	else if (ODE_solver == 1) {//RungeKutta
		int actsteps = 0; double l = 0.;
#ifndef fastRK4
		const int order = 4; const int variation = 0;
		DVec c_coeff(order); DVec b_coeff(order); DMat a_coeff(order, DVec(order));
		rungeKutta(c_coeff, b_coeff, a_coeff, order, variation);
		vector<State> auxStates(order); State tempState;
#endif		
		for (int m = 0; m < steps; ++m) {
			if (LambdaL(l, Lambda_0) < Lambda_f) break;
			actsteps += 1;
#ifdef fastRK4
			State k1 = derivative(LambdaL(l, Lambda_0), state) * ( dLambdaL(l, Lambda_0) * h );
			State k2 = derivative(LambdaL(l+h/2., Lambda_0), state + (k1*0.5) ) * ( dLambdaL(l+h/2., Lambda_0) * h );
			State k3 = derivative(LambdaL(l+h/2., Lambda_0), state + (k2*0.5) ) * ( dLambdaL(l+h/2., Lambda_0) * h );
			State k4 = derivative(LambdaL(l+h, Lambda_0), state + k3) * ( dLambdaL(l+h, Lambda_0) * h );
			state = state + ( k1 + (k2*2.) + (k3*2.) + k4 ) * (1./6.);
#else
			for (int i = 0; i < order; ++i) {
				tempState = state;
				for (int j = 0; j < i; ++j) tempState = tempState + auxStates[i-1] * a_coeff[i][j];
				auxStates[i] = derivative(LambdaL(l+h*c_coeff[i], Lambda_0), tempState) * ( dLambdaL(l+h*c_coeff[i], Lambda_0) * h );
			}
			for (int i = 0; i < order; ++i) state = state + move(auxStates[i]) * b_coeff[i];
#endif
			l += h;
		}
		// last step
		State memoryState; if (print or true) memoryState = state;
		State deriv = derivative(LambdaL(l, Lambda_0), state);
		state = state + deriv * (0.-LambdaL(l, Lambda_0));
		if (print or true) cout << "steps " << actsteps << " change in last step " << getError(memoryState, state) << endl;
	}
	else if (ODE_solver == 2 or ODE_solver == 3) {//AdamsBashforth, AdamsMoulton
		State memoryState; State compareState; 
		const int orderB = 6; const int orderMmB = 0; //orderMmB can only be 0 or 1, 0 is recommended
		DVec ab_coeff(orderB); DVec am_coeff(orderB+orderMmB); DVec abCompare_coeff(orderB-1);
		vector<State> auxStates( orderB+1 ); 
		adamsB(ab_coeff, orderB);
	   	if (ODE_solver == 3) adamsM(am_coeff, orderB+orderMmB); 
		if (print) adamsB(abCompare_coeff, orderB-1); 
		const int rkdiv = 2;
		double l = 0.;
		
		DVec L0(steps); DVec l0(steps); CVec C0(steps); CVec C1(steps); CVec C2(steps); CVec C3(steps); CVec C4(steps); CVec Pi(1);
		
		h = h/rkdiv; l = l-h*( 1.-1./((double)rkdiv) ); // l = l-h*(rkdiv-1);
		if (print) cout << "start at Lambda = " << LambdaL(l, Lambda_0) << endl;
		for (int m = 0; m < orderB*rkdiv; ++m) {//initializer, RK4
			auxStates[m/rkdiv] = derivative(LambdaL(l, Lambda_0), state) * ( dLambdaL(l, Lambda_0) * h );
			State k2 = derivative(LambdaL(l+h/2., Lambda_0), state + (auxStates[m/rkdiv]*0.5) ) * ( dLambdaL(l+h/2., Lambda_0) * h );
			State k3 = derivative(LambdaL(l+h/2., Lambda_0), state + (k2*0.5) ) * ( dLambdaL(l+h/2., Lambda_0) * h );
			State k4 = derivative(LambdaL(l+h, Lambda_0), state + k3) * ( dLambdaL(l+h, Lambda_0) * h );
			if (m < orderB*rkdiv-1) {
				state = state + ( auxStates[m/rkdiv] + (k2*2.) + (k3*2.) + k4 ) * (1./6.);
				l += h;
			}
			if (print and m==orderB*rkdiv-1) compareState = state + ( auxStates[m/rkdiv] + (k2*2.) + (k3*2.) + k4 ) * ( ((double)rkdiv) *1./6.);
		}
		
		h = h*rkdiv; l = h*(orderB-1); 
		for (int m = 0; m < steps; ++m) {
			if (ODE_solver == 3) memoryState = state;
			for (int i = 0; i < orderB; ++i) state = state + auxStates[i] * ab_coeff[i];
			if (print and m==0) cout << "difference between RK and first Adams " << getError(compareState, state) << endl;
			if (print and m%10==0) {
				compareState = state;
				for (int i = 0; i < orderB-1; ++i) compareState = compareState + auxStates[i+1] * abCompare_coeff[i];
			}
				
			l += h;
			
			l0[m] = l; L0[m] = LambdaL(l, Lambda_0); state.getCorr(Pi); C0[m] = Pi[0]; state.getCorr(Pi, LambdaL(l, Lambda_0)); C1[m] = Pi[0], C2[m] = state.Ip.K1[0]; C3[m] = state.Ip.K2[0][NUMFREQ2]; C4[m] = state.Ip.K3[0][triangularIndex(NUMFREQ3, NUMFREQ3, 2*NUMFREQ3)];
			
			if (ODE_solver == 3) {
				auxStates[orderB] = derivative(LambdaL(l, Lambda_0), state) * ( dLambdaL(l, Lambda_0) * h );
				state = move(memoryState);
				for (int i = 0; i < orderB+orderMmB; ++i) state = state + auxStates[i+1-orderMmB] * am_coeff[i];
			}
						
			if (print and m%10==0) cout << "iteration at step " << m << " at Lambda " << LambdaL(l, Lambda_0) << " at Lambda*beta/Pi " << LambdaL(l, Lambda_0)/PiOverBeta << "\tdifference between Adams and Adams one order lower " << getError(compareState, state) << endl;
			
			if (LambdaL(l, Lambda_0) < Lambda_f) break;
			
			for (int i = 0; i < orderB-1; ++i) auxStates[i] = auxStates[i+1];//move( auxStates[i+1] );
			 
			auxStates[orderB-1] = derivative(LambdaL(l, Lambda_0), state) * ( dLambdaL(l, Lambda_0) * h );
		}		
		// last step
		if (print) memoryState = state;
		State deriv = derivative(LambdaL(l, Lambda_0), state);
		state = state + deriv * (0.-LambdaL(l, Lambda_0));
		if (print or true) cout << "change in last step " << getError(memoryState, state) << endl;
		writeFlow(suffix, l0, L0, C0, C1, C2, C3, C4);
	}
	
}

#include "mytests.h"
int main(int argc, const char * argv[]) {
	//bool test_only = false; if(test_only) {bool print = true; test_run(print); return 0;}
	double Time1, Time2;
	Time1 = omp_get_wtime();
	
	string suffix; int steps, ODE_solver; double Lambda_0, Lambda_f_fraction;
    argc > 1 ? suffix = argv[1] : suffix="test_U28"; // suffix for file name
    argc > 2 ? steps = (int)strtod(argv[2], NULL) : steps=600;
    argc > 3 ? ODE_solver = (int)strtod(argv[3], NULL) : ODE_solver=2;
    argc > 4 ? Lambda_0 = strtod(argv[4], NULL) : Lambda_0=5000.;
    argc > 5 ? Lambda_f_fraction = strtod(argv[5], NULL) : Lambda_f_fraction=0.002;

	State mystate;
	
#if MODE==0
	parquetSolver(mystate, steps);
#else 
	propagate(mystate, Lambda_0, Lambda_f_fraction, steps, ODE_solver, true, suffix);
#endif
	
	mystate.write(suffix);
	
	Time2 = omp_get_wtime();
	printf("Program took %.5f seconds.\n", Time2-Time1);
	return 0;
}
