#include <cstdlib>
#include<bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <complex>
#include "parameters.h"
#include "vertex.h"
#include "state.h"
#include "loop.h"
#include "a_bubble.h"
#include "p_bubble.h"
#include "t_bubble.h"
#include "propagator.h"
#include "selfenergy.h"

using namespace std;

typedef complex<double> comp;

void writeOutFile(rvec& freqs, double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy);
template <typename Q> State<Q> derivative(double Lambda, State<Q>& state);



int main() {

    for(int i=0; i<nSE; ++i)
    {
        bfreqs[i] = w_lower_b + i*dw;
        ffreqs[i] = w_lower_f + i*dv;

        freqs_a[i] = w_lower_b + i*dw_a;
        freqs_p[i] = w_lower_b + i*dw_p;
        freqs_t[i] = w_lower_b + i*dw_t;

        simpson_weights[i] = (double)(2+2*(i%2));
    }
    simpson_weights[0] = 1.;
    simpson_weights[nSE-1] = 1.;

    for(int i=0; i<nEVO; ++i)
    {
        flow_grid[i] = Lambda_ini + i*dL;
    }



    double  t0 = get_time();
    State<comp> state(Lambda_ini);
    cout << "State created. ";
    get_time(t0);

    //Initial conditions
    for (int i = 0; i < nSE; ++i) {
        state.selfenergy.setself(0, i, 0.5 * U);
        state.selfenergy.setself(1, i, 0.);
        state.diffselfenergy.setself(0, i, 0.);
        state.diffselfenergy.setself(1, i, 0.);
    }
    cout << "self energy and diff self energy assigned" << endl;

    for (auto i:odd_Keldysh) {
        state.vertex.densvertex.irred.setvert(i,  0.5*U);
        state.vertex.spinvertex.irred.setvert(i, -0.5*U);
    }
    cout << "vertex assigned" << endl;


    cout << "Start of flow" << endl;
    for(auto Lambda:flow_grid) {
        State<comp> dPsi = derivative(Lambda, state);

        state += dPsi;

        Propagator control = propag(Lambda, state.selfenergy, state.diffselfenergy, 'g');

        writeOutFile(ffreqs, Lambda, control, state.selfenergy);
    }


    return 0;
}


template <typename Q>
State<Q> derivative(double Lambda, State<Q>& state)
{
    State<Q> resp(Lambda);

    /*Here I begin implementing Fabian's pseudocode*/
    //Line 1
    Propagator S = propag(Lambda, state.selfenergy, state.diffselfenergy, 's');
    cout << "S calculated" << endl;

    //Line 2
    Propagator G = propag(Lambda, state.selfenergy, state.diffselfenergy, 'g');
    cout << "G calculated" << endl;

    //Line 3
    SelfEnergy<comp> Sigma_std = loop(state.vertex, S);
    cout << "loop calculated" << endl;

    //Line 4
    resp.selfenergy = Sigma_std;

    //Line 6
    Propagator dG = propag(Lambda, resp.selfenergy, resp.diffselfenergy, 'k');

    cout << "diff bubble started" << endl;
    double t2 = get_time();
    //Lines 7, 8 & 9
    Vertex<avert<comp> > dgammaa = diff_a_bubble_function(state.vertex, state.vertex, G, dG);
    Vertex<pvert<comp> > dgammap = diff_p_bubble_function(state.vertex, state.vertex, G, dG);
    Vertex<tvert<comp> > dgammat = diff_t_bubble_function(state.vertex, state.vertex, G, dG);
    cout << "diff bubble finished. ";
    get_time(t2);


//    resp.vertex.densvertex.irred = state.vertex.densvertex.irred;
    resp.vertex.densvertex.avertex = dgammaa.densvertex;
    resp.vertex.densvertex.pvertex = dgammap.densvertex;
    resp.vertex.densvertex.tvertex = dgammat.densvertex;
    return resp;
}


//Writes of .dat files of the propagator, the self energy and selected values of the vertex for different values of Lambda during the fRG flow
void writeOutFile(rvec& freqs, double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy)
{
    int i = fconv_Lambda(Lambda);

    ostringstream self_energyR, self_energyA, self_energyK, propR, propA, propK;
    self_energyR << "self_energyR"<<i<<".dat";
    self_energyA << "self_energyA"<<i<<".dat";
    self_energyK << "self_energyK"<<i<<".dat";

    propR << "propagatorR"<<i<<".dat";
    propA << "propagatorA"<<i<<".dat";
    propK << "propagatorK"<<i<<".dat";

    ofstream my_file_sigmaR, my_file_sigmaA, my_file_sigmaK, my_file_propR, my_file_propA, my_file_propK;
    my_file_sigmaR.open(self_energyR.str());
    my_file_sigmaA.open(self_energyA.str());
    my_file_sigmaK.open(self_energyK.str());

    my_file_propR.open(propR.str());
    my_file_propA.open(propA.str());
    my_file_propK.open(propK.str());


    for (int j = 0; j < ffreqs.size(); j++) {
        my_file_sigmaR << ffreqs[j] << " " << selfEnergy.svalsmooth(0, ffreqs[j]).real() << " " <<  selfEnergy.svalsmooth(0,freqs[j]).imag() << "\n";
        my_file_sigmaA << ffreqs[j] << " " << selfEnergy.svalsmooth(0, ffreqs[j]).real() << " " << -selfEnergy.svalsmooth(0,freqs[j]).imag() << "\n";
        my_file_sigmaK << ffreqs[j] << " " << selfEnergy.svalsmooth(1, ffreqs[j]).real() << " " <<  selfEnergy.svalsmooth(1,freqs[j]).imag() << "\n";

        my_file_propR << ffreqs[j] << " " << propagator.pvalsmooth(0, ffreqs[j]).real() << " " <<  propagator.pvalsmooth(0, ffreqs[j]).imag() << "\n";
        my_file_propA << ffreqs[j] << " " << propagator.pvalsmooth(0, ffreqs[j]).real() << " " << -propagator.pvalsmooth(0, ffreqs[j]).imag() << "\n";
        my_file_propK << ffreqs[j] << " " << propagator.pvalsmooth(1, ffreqs[j]).real() << " " <<  propagator.pvalsmooth(1, ffreqs[j]).imag() << "\n";
    }
    my_file_sigmaR.close();
    my_file_sigmaA.close();
    my_file_sigmaK.close();

    my_file_propR.close();
    my_file_propA.close();
    my_file_propK.close();
}
