#include <cstdlib>
#include <bits/stdc++.h> // TODO: what is this needed for?
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
void setUpBosGrid();
void setUpFerGrid();
void setUpFlowGrid();

template <typename Q> void SOPT(double Lambda, State<Q>& state);
void writeOutSOPT(Propagator& propagator, SelfEnergy<comp>& selfEnergy);


int main() {

    setUpBosGrid();
    setUpFerGrid();
    setUpFlowGrid();


    double  t0 = get_time();
    State<comp> state(Lambda_ini);
    cout << "State created. ";
    get_time(t0);

    //Initial conditions
    for (int i = 0; i < nSE; ++i) {
        state.selfenergy.setself(0, i, 0.5*U);
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
        double tder = get_time();
        state.Lambda = Lambda;
        State<comp> dPsi = derivative(Lambda, state)*dL;

        Propagator control = propag(Lambda, state.selfenergy, state.diffselfenergy, 'g','f');

        double tadd = get_time();

        //TODO Check addition operator for Vertex
//        state = state + dPsi;
        state += dPsi;
        cout << "Added:";
        get_time(tadd);

        writeOutFile(ffreqs, Lambda, control, state.selfenergy);
        cout << "Wrote out" <<endl;
        cout << "One derivatie step: ";
        get_time(tder);
    }

//    SOPT(1., state);
//    Propagator control = propag(1., state.selfenergy, state.diffselfenergy, 'g', '.');
//    writeOutSOPT(control, state.selfenergy);


    return 0;
}



template <typename Q>
State<Q> derivative(double Lambda, State<Q>& state)
{
    State<Q> resp(Lambda);

    /*Here I begin implementing Fabian's pseudocode*/
    //Line 1
    Propagator S = propag(Lambda, state.selfenergy, state.diffselfenergy, 's', '.');
    cout << "S calculated" << endl;

    //Line 2
    Propagator G = propag(Lambda, state.selfenergy, state.diffselfenergy, 'g', '.');
    cout << "G calculated" << endl;

    //Line 3
    SelfEnergy<comp> Sigma_std = loop(state.vertex, S);
    cout << "loop calculated" << endl;

    //Line 4
    resp.selfenergy = Sigma_std;

    //Line 6
    Propagator dG = propag(Lambda, resp.selfenergy, resp.diffselfenergy, 'k', '.');

    cout << "diff bubble started" << endl;
    double t2 = get_time();
    //Lines 7-9
    double ta = get_time();
    Vertex<avert<comp> > dgammaa = diff_a_bubble_function(state.vertex, state.vertex, G, dG);
    cout<<  "a - Bubble:";
    get_time(ta);

    double tp = get_time();
    Vertex<pvert<comp> > dgammap = diff_p_bubble_function(state.vertex, state.vertex, G, dG);
    cout<<  "p - Bubble:";
    get_time(tp);

    double  tt=get_time();
    Vertex<tvert<comp> > dgammat = diff_t_bubble_function(state.vertex, state.vertex, G, dG);
    cout<<  "t - Bubble:";
    get_time(tt);

    cout << "diff bubble finished. ";
    get_time(t2);


//    double t3 = get_time();
////    Lines 10-13
//    Vertex<fullvert<Q> > dGamma = Vertex<fullvert<Q> >();
//    dGamma.densvertex.avertex = dgammaa.densvertex;
//    dGamma.densvertex.pvertex = dgammap.densvertex;
//    dGamma.densvertex.tvertex = dgammat.densvertex;
//    dGamma.spinvertex.avertex = dgammaa.spinvertex;
//    dGamma.spinvertex.pvertex = dgammap.spinvertex;
//    dGamma.spinvertex.tvertex = dgammat.spinvertex;
//    cout<< "dGamma assigned: " <<endl;
//    get_time(t3);
//    Vertex<fullvert<Q> > dgammaabar = dgammap + dgammat;
//    Vertex<fullvert<Q> > dgammapbar = dgammat + dgammaa;
//    Vertex<fullvert<Q> > dgammatbar = dgammaa + dgammap;

//    double t4 = get_time();
//    //The r_bubble_function pics the gamma_r bar contributions
//    Vertex<avert<Q> > dgammaLa = a_bubble_function(dGamma, dGamma, G, 'L');
//    Vertex<pvert<Q> > dgammaLp = p_bubble_function(dGamma, dGamma, G, 'L');
//    Vertex<tvert<Q> > dgammaLt = t_bubble_function(dGamma, dGamma, G, 'L');
//
//    Vertex<avert<Q> > dgammaRa = a_bubble_function(dGamma, dGamma, G, 'R');
//    Vertex<pvert<Q> > dgammaRp = p_bubble_function(dGamma, dGamma, G, 'R');
//    Vertex<tvert<Q> > dgammaRt = t_bubble_function(dGamma, dGamma, G, 'R');
//
//    cout<< "Bubbles calculated: " << endl;
//    get_time(t4);


    //Line 14-17
//    Vertex<fullvert<Q> > dGammaT = Vertex<fullvert<Q> >();
//    dGammaT.densvertex.avertex = dgammaLa.densvertex + dgammaRa.densvertex;
//    dGammaT.densvertex.pvertex = dgammaLp.densvertex + dgammaRp.densvertex;
//    dGammaT.densvertex.tvertex = dgammaLt.densvertex + dgammaRt.densvertex;
//    dGammaT.spinvertex.avertex = dgammaLa.spinvertex + dgammaRa.spinvertex;
//    dGammaT.spinvertex.pvertex = dgammaLp.spinvertex + dgammaRp.spinvertex;
//    dGammaT.spinvertex.tvertex = dgammaLt.spinvertex + dgammaRt.spinvertex;
//
//    dgammaa.densvertex += dGammaT.densvertex.avertex;
//    dgammap.densvertex += dGammaT.densvertex.pvertex;
//    dgammat.densvertex += dGammaT.densvertex.tvertex;
//    dgammaa.spinvertex += dGammaT.spinvertex.avertex;
//    dgammap.spinvertex += dGammaT.spinvertex.pvertex;
//    dgammat.spinvertex += dGammaT.spinvertex.tvertex;


    //Line 41
    resp.vertex.densvertex.avertex = dgammaa.densvertex;
    resp.vertex.densvertex.pvertex = dgammap.densvertex;
    resp.vertex.densvertex.tvertex = dgammat.densvertex;
    resp.vertex.spinvertex.avertex = dgammaa.spinvertex;
    resp.vertex.spinvertex.pvertex = dgammap.spinvertex;
    resp.vertex.spinvertex.tvertex = dgammat.spinvertex;
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


    for (int j = 0; j < freqs.size(); j++) {
        my_file_sigmaR << freqs[j] << " " << selfEnergy.sval(0, j).real() << " " <<  selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaA << freqs[j] << " " << selfEnergy.sval(0, j).real() << " " << -selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaK << freqs[j] << " " << selfEnergy.sval(1, j).real() << " " <<  selfEnergy.sval(1, j).imag() << "\n";

        my_file_propR << freqs[j] << " " << propagator.pval(0, j).real() << " " <<  propagator.pval(0, j).imag() << "\n";
        my_file_propA << freqs[j] << " " << propagator.pval(0, j).real() << " " << -propagator.pval(0, j).imag() << "\n";
        my_file_propK << freqs[j] << " " << propagator.pval(1, j).real() << " " <<  propagator.pval(1, j).imag() << "\n";
    }
    my_file_sigmaR.close();
    my_file_sigmaA.close();
    my_file_sigmaK.close();

    my_file_propR.close();
    my_file_propA.close();
    my_file_propK.close();
}

void setUpBosGrid()
{
#if GRID==1

#elif GRID==2
    for(int i=0; i<nBOS; ++i)
        bfreqs[i] = w_lower_b + i*dw;

#endif
}
void setUpFerGrid()
{
#if GRID==1

#elif GRID==2
    for(int i=0; i<nFER; ++i)
        ffreqs[i] = w_lower_f + i*dv;
#endif
}
void setUpFlowGrid()
{
    for(int i=0; i<nEVO; ++i)
        flow_grid[i] = Lambda_ini + i*dL;
}

template<typename Q>
void SOPT(double Lambda, State<Q> &state) {

    Propagator g = propag(Lambda, state.selfenergy, state.diffselfenergy, 'g', 'f');
    cout << "G calculated" << endl;


    cout << "bubble started" << endl;
    double t2 = get_time();
    //Lines 7-9
    double ta = get_time();
    Vertex<avert<comp> > dgammaa = a_bubble_function(state.vertex, state.vertex, g, '.');
    cout<<  "a - Bubble:";
    get_time(ta);

    double tp = get_time();
    Vertex<pvert<comp> > dgammap = p_bubble_function(state.vertex, state.vertex, g, '.');
    cout<<  "p - Bubble:";
    get_time(tp);

    double tt = get_time();
    Vertex<tvert<comp> > dgammat = t_bubble_function(state.vertex, state.vertex, g, '.');
    cout<<  "t - Bubble:";
    get_time(tt);

    cout << "bubble finished. ";
    get_time(t2);

    Propagator s = propag(Lambda, state.selfenergy, state.diffselfenergy, 's', 'f');

    //Line 41
    state.vertex.densvertex.avertex += dgammaa.densvertex;
    state.vertex.densvertex.pvertex += dgammap.densvertex;
    state.vertex.densvertex.tvertex += dgammat.densvertex;


    double tloop = get_time();
    SelfEnergy<comp> Sigma_std = loop(state.vertex, s);
    cout << "loop calculated ";
    get_time(tloop);

    state.selfenergy += Sigma_std;
}

void writeOutSOPT(Propagator& propagator, SelfEnergy<comp>& selfEnergy)
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

    for (int j = 0; j < nFER; j++) {
        my_file_sigmaR << ffreqs[j] << " " << selfEnergy.sval(0, j).real() << " " <<  selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaA << ffreqs[j] << " " << selfEnergy.sval(0, j).real() << " " << -selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaK << ffreqs[j] << " " << selfEnergy.sval(1, j).real() << " " <<  selfEnergy.sval(1, j).imag() << "\n";

        my_file_propR << ffreqs[j] << " " << propagator.pval(0, j).real() << " " <<  propagator.pval(0, j).imag() << "\n";
        my_file_propA << ffreqs[j] << " " << propagator.pval(0, j).real() << " " << -propagator.pval(0, j).imag() << "\n";
        my_file_propK << ffreqs[j] << " " << propagator.pval(1, j).real() << " " <<  propagator.pval(1, j).imag() << "\n";
    }
    my_file_sigmaR.close();
    my_file_sigmaA.close();
    my_file_sigmaK.close();

    my_file_propR.close();
    my_file_propA.close();
    my_file_propK.close();
}