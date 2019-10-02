#include <cstdlib>
#include<bits/stdc++.h> // TODO: what is this needed for?
#include <iostream>
#include <fstream>
#include <complex>
#include "parameters.h"
#include "vertex.h"
#include "state.h"
#include "loop.h"
#include "a_bubble.h"
#include "propagator.h"
#include "selfenergy.h"

using namespace std;

typedef complex<double> comp;

void writeOutFile(rvec& freqs, SelfEnergy<comp>& selfEnergy, Propagator& propagator);




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


    double Lambda = 1.0;

    State state(Lambda);


    //Initial conditions
    SelfEnergy<comp> self_ini;
    SelfEnergy<comp> diff_self_ini;
    for(int i=0; i<nSE; ++i)
    {
        self_ini.setself(0, i, 0.5*U);
        self_ini.setself(1, i, 0.);
        diff_self_ini.setself(0, i, 0.);
        diff_self_ini.setself(1, i, 0.);
    }

    double t0 = get_time();

    Vertex<fullvert<comp> > vertex = Vertex<fullvert<comp> >();
    get_time(t0);

    cout << "vertex calculated" << endl;
    for(auto i:odd_Keldysh)
    {
        vertex.densvertex.irred.setvert(i, 0.5*U);
        vertex.spinvertex.irred.setvert(i,-0.5*U);
    }
    state.vertex = vertex;
    cout << "vertex assigned" << endl;

    /*Here I begin implementing Fabian's pseudocode*/
    //Line 1
    Propagator S = propag(Lambda, self_ini, diff_self_ini, 's');
    cout << "S calculated" << endl;

    //Line 2
    Propagator G = propag(Lambda, self_ini, diff_self_ini, 'g');
    cout << "G calculated" << endl;


    //Line 3
    SelfEnergy<comp> Sigma_std = loop(vertex, S);
    cout << "loop calculated" << endl;


//    //Line 4
//    state.selfenergy = Sigma_std;
//
//    //Line 6
//    Propagator dG = propag(Lambda, state.selfenergy, diff_self_ini, 'k');
//
//    cout << "diff bubble started" << endl;
//    //Line 8
//    avert<comp> dgammaa = diff_a_bubble(state.vertex.densvertex.avertex, state.vertex.densvertex.avertex, G, dG);
//    cout << "diff bubble finished" << endl;
//
//
//
    writeOutFile(ffreqs, Sigma_std, G);


    return 0;
}




//Writes of .dat files of the propagator, the self energy and selected values of the vertex for different values of Lambda during the fRG flow
void writeOutFile(rvec& freqs, SelfEnergy<comp>& selfEnergy, Propagator& propagator)
{
        ostringstream self_energyR, self_energyA, self_energyK, propR, propA, propK;
        self_energyR << "self_energyR.dat";
        self_energyA << "self_energyA.dat";
        self_energyK << "self_energyK.dat";

        propR << "propagatorR.dat";
        propA << "propagatorA.dat";
        propK << "propagatorK.dat";

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
