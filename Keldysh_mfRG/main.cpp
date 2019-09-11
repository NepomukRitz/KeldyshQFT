#include <cstdlib>
#include <iostream>
#include <fstream>
#include <complex>

#include "vertex.h"
//#include "bubbles.h"
#include "state.h"
#include "a_bubble.h"
#include "propagator.h"
#include "selfenergy.h"

using namespace std;

typedef complex<double> comp;

vector<double> ffreqs1 (nSE);                                                                                                      // NOLINT(cert-err58-cpp)

//TODO figure out why ffreqs is not exactly how we want it

void writeOutFile(rvec& freqs, SelfEnergy<comp>& selfEnergy, Propagator& propagator);
bool isAntisymmetric(Propagator& propagator);
bool isAntisymmetric();


int main() {

    auto dw = (w_upper_b-w_lower_b)/((double)(nSE-1));
    auto dv = (w_upper_f-w_lower_f)/((double)(nSE-1));

    for(int i=0; i<nSE; ++i)
    {
        bfreqs[i] = w_lower_b + i*dw;
        ffreqs[i] = w_lower_f + i*dv;
        ffreqs1[i] = w_lower_f + i*dv;


        freqs_a[i] = w_lower_b + i*dw_a;
        freqs_p[i] = w_lower_b + i*dw_p;
        freqs_t[i] = w_lower_b + i*dw_t;

        simpson_weights[i] = (double)(2+2*(i%2));

    }
    simpson_weights[0] = 1.;
    simpson_weights[nSE-1] = 1.;


    double Lambda = 1.0;

//    state State(Lambda);


    //Initial conditions
    SelfEnergy<comp> self_ini;
    SelfEnergy<comp> diff_self_ini;
    for(int i=0; i<nSE; ++i)
    {
        self_ini.setself(0, i, U/2.);
        self_ini.setself(1, 0, 0.);
        diff_self_ini.setself(0, i, 0.);
        diff_self_ini.setself(1, i, 0.);
    }
//    Vertex<fullvert<comp> > vertex = Vertex<fullvert<comp> >();
//
//    for(auto i:odd_Keldysh)
//    {
//        vertex.densvertex.irred.setvert(i, 0.5*U);
//        vertex.spinvertex.irred.setvert(i,-0.5*U);
//    }

    /*Here I begin implementing Fabian's pseudocode*/
    //Line 1
    Propagator S = propag(Lambda, self_ini, diff_self_ini, 's');
    //Line 2
    Propagator G = propag(Lambda, self_ini, diff_self_ini, 'g');
    //Line 3
//    SelfEnergy<comp> Sigma_std = State.loop(vertex, S);
    //Line4
//    State.selfenergy = Sigma_std;


//    writeOutFile(ffreqs, Sigma_std, S);


//    bool aid1 = isAntisymmetric(S);
//    bool aid2 = isAntisymmetric(G);
//
//    cout << "aid1: " << aid1 << " and aid2: " << aid2 << endl;

    bool aidusmaximus = isAntisymmetric();

//    test.spinvertex.avertex.K1_setvert(0, 1, 2, 3);
//
//    Vertex<avert<comp> > test1;
//    Vertex<avert<comp> > test2;
//    Vertex<avert<comp> > test3;

//    test3 = test1 + test2;

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

bool isAntisymmetric(Propagator& propa)
{
    bool resp = true;
    for(int i=0; i<nPROP; ++i)
    {
        if(propa.pval(0,i).real() + propa.pval(0, nPROP-1-i).real()>10e-3)
            cout << propa.pval(0,i).real() + propa.pval(0, nPROP-1-i).real()<<endl;
            resp = false;
    }
    return resp;
}

bool isAntisymmetric()
{
    bool resp = true;
    for(int i=0; i<ffreqs1.size(); ++i)
    {
//        cout << ffreqs[i] << " " << -ffreqs[ffreqs.size()-1-i] << endl;
        if(ffreqs[i] != -ffreqs[ffreqs.size()-1-i]) {
            resp = false;
            cout << ffreqs [i]+ ffreqs [ffreqs.size() - 1 - i] << endl;
            cout << ffreqs1[i]+ ffreqs1[ffreqs1.size() - 1 - i] << endl;
        }
    }
    return resp;
}