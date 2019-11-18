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

void writeOutFile(rvec& freqs, double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex);
template <typename Q> void derivative(State<Q>& dPsi, double Lambda, State<Q>& state);
template <typename Q> void RungeKutta4thOrder(State<Q>& dPsi, double Lambda, State<Q>& state);

template <typename Q> void SOPT(State<Q>& bare, double Lambda, State<Q>& state);
void writeOutSOPT(double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex);


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
    }
    cout << "self energy and diff self energy assigned" << endl;

    for (auto i:odd_Keldysh) {
        state.vertex.densvertex.irred.setvert(i,  0.5*U);
        state.vertex.spinvertex.irred.setvert(i, -0.5*U);
    }
    cout << "vertex assigned" << endl;


    for(int i=0; i<nEVO; ++i) {
        double Lambda = flow_grid[i];

        State<comp> bare (Lambda);
        for (auto j:odd_Keldysh) {
            bare.vertex.densvertex.irred.setvert(j,  0.5*U);
        }

        SOPT(bare, Lambda, state);

        Propagator control = propag(Lambda, state.selfenergy, 'g', '.');
        writeOutSOPT(Lambda, control, state.selfenergy, bare.vertex);
    }



//    Propagator initial = propag(state.Lambda, state.selfenergy, 'g', '.');
//    Propagator compare = propag(state.Lambda, state.selfenergy, 'g', 'f');
//
//    writeOutFile(ffreqs, state.Lambda, initial, state.selfenergy, state.vertex);
//
//    cout << "Start of flow" << endl;
//    for(int i=1; i<nEVO; ++i) {
//        double Lambda = flow_grid[i-1];
//        double tder = get_time();
//        state.Lambda = Lambda;
//
//        State<comp> dPsi(Lambda);
//
//        derivative(dPsi, Lambda, state);
////        RungeKutta4thOrder(dPsi, Lambda, state);
//
//
//        double tadd = get_time();
//        state += dPsi;
//        cout << "Added:";
//        get_time(tadd);
//
//        double next_Lambda = flow_grid[i];
//
//        Propagator control = propag(next_Lambda, state.selfenergy, 'g','.');
//
//        writeOutFile(ffreqs, next_Lambda, control, state.selfenergy, state.vertex);
//        cout << "Wrote out" <<endl;
//        cout << "One RK-derivative step: ";
//        get_time(tder);
//    }
//
//    cout << "Total execution time: ";
//    get_time(t0);


    return 0;
}



template <typename Q>
void derivative(State<Q>& dPsi, double Lambda, State<Q>& state)
{
    /*Here I begin implementing Fabian's pseudocode*/
    //Line 1
    Propagator S = propag(Lambda, state.selfenergy, 's', 'f');
    cout << "S calculated" << endl;

    //Line 2
    Propagator G = propag(Lambda, state.selfenergy, 'g', 'f');
    cout << "G calculated" << endl;

    //Line 3
    SelfEnergy<comp> Sigma_std = loop(state.vertex, S);
    cout << "loop calculated" << endl;

    //Line 4
    dPsi.selfenergy=Sigma_std;

    //Line 6
    Propagator extension = propag(Lambda, dPsi.selfenergy, 'e', 'f');\
    Propagator dG = S + extension;

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

//    Lines 10-13
//    double t4 = get_time();
//    //The r_bubble_function pics the gamma_r bar contributions
//    Vertex<avert<Q> > dgammaLa = a_bubble_function(state.vertex, state.vertex, G, 'L');
//    Vertex<pvert<Q> > dgammaLp = p_bubble_function(state.vertex, state.vertex, G, 'L');
//    Vertex<tvert<Q> > dgammaLt = t_bubble_function(state.vertex, state.vertex, G, 'L');
//
//    Vertex<avert<Q> > dgammaRa = a_bubble_function(state.vertex, state.vertex, G, 'R');
//    Vertex<pvert<Q> > dgammaRp = p_bubble_function(state.vertex, state.vertex, G, 'R');
//    Vertex<tvert<Q> > dgammaRt = t_bubble_function(state.vertex, state.vertex, G, 'R');
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
    dPsi.vertex.densvertex.avertex = dgammaa.densvertex;
    dPsi.vertex.densvertex.pvertex = dgammap.densvertex;
    dPsi.vertex.densvertex.tvertex = dgammat.densvertex;
    dPsi.vertex.spinvertex.avertex = dgammaa.spinvertex;
    dPsi.vertex.spinvertex.pvertex = dgammap.spinvertex;
    dPsi.vertex.spinvertex.tvertex = dgammat.spinvertex;

    dPsi*=dL;
}

template <typename Q> void RungeKutta4thOrder(State<Q>& dPsi, double Lambda, State<Q>& state)
{
    State<Q> k1(Lambda);
    State<Q> k2(Lambda);
    State<Q> k3(Lambda);
    State<Q> k4(Lambda);

    derivative(k1, Lambda, state);

    k1 *= dL/2.;
    k1 += state;

    derivative(k2, Lambda + dL/2., k1);

    k2 *= dL/2.;
    k2 += state;

    derivative(k3, Lambda + dL/2., k2);

    k3 *= dL;
    k3 += state;

    derivative(k4, Lambda + dL, k3);

    k1 -= state;
    k1 *= 2./dL;
    k2 -= state;
    k2 *= 2./dL;
    k3 -= state;
    k3 *= 1./dL;

    dPsi += (k1 + k2*2. + k3*2. + k4)*(1./6.);

}

//Writes of .dat files of the propagator, the self energy and selected values of the vertex for different values of Lambda during the fRG flow
void writeOutFile(rvec& freqs, double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex)
{
    int i = fconv_Lambda(Lambda);

    ostringstream self_energyR, self_energyA, self_energyK, propR, propA, propK;
    ostringstream avert1, pvert1, tvert1, avert3, pvert5, tvert3;
    self_energyR << "self_energyR"<<i<<".dat";
    self_energyA << "self_energyA"<<i<<".dat";
    self_energyK << "self_energyK"<<i<<".dat";

    propR << "propagatorR"<<i<<".dat";
    propA << "propagatorA"<<i<<".dat";
    propK << "propagatorK"<<i<<".dat";

    avert1 << "avert1"<<i<<".dat";
    pvert1 << "pvert1"<<i<<".dat";
    tvert1 << "tvert1"<<i<<".dat";

    avert3 << "avert3"<<i<<".dat";
    pvert5 << "pvert5"<<i<<".dat";
    tvert3 << "tvert3"<<i<<".dat";

    ofstream my_file_sigmaR, my_file_sigmaA, my_file_sigmaK, my_file_propR, my_file_propA, my_file_propK;
    ofstream  my_file_avert1, my_file_pvert1, my_file_tvert1, my_file_avert3, my_file_pvert5, my_file_tvert3;
    my_file_sigmaR.open(self_energyR.str());
    my_file_sigmaA.open(self_energyA.str());
    my_file_sigmaK.open(self_energyK.str());

    my_file_propR.open(propR.str());
    my_file_propA.open(propA.str());
    my_file_propK.open(propK.str());

    my_file_avert1.open(avert1.str());
    my_file_pvert1.open(pvert1.str());
    my_file_tvert1.open(tvert1.str());

    my_file_avert3.open(avert3.str());
    my_file_pvert5.open(pvert5.str());
    my_file_tvert3.open(tvert3.str());

    for (int j = 0; j < freqs.size(); j++) {
        my_file_sigmaR << freqs[j] << " " << selfEnergy.sval(0, j).real() << " " <<  selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaA << freqs[j] << " " << selfEnergy.sval(0, j).real() << " " << -selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaK << freqs[j] << " " << selfEnergy.sval(1, j).real() << " " <<  selfEnergy.sval(1, j).imag() << "\n";

        my_file_propR << freqs[j] << " " << propagator.pval(0, j).real() << " " <<  propagator.pval(0, j).imag() << "\n";
        my_file_propA << freqs[j] << " " << propagator.pval(0, j).real() << " " << -propagator.pval(0, j).imag() << "\n";
        my_file_propK << freqs[j] << " " << propagator.pval(1, j).real() << " " <<  propagator.pval(1, j).imag() << "\n";
    }

    for (int j = 0; j<bfreqs.size(); j++){
        my_file_avert1 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_pvert1 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_tvert1 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).imag()<< "\n";

        my_file_avert3 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(3, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(3, j, 0).imag()<< "\n";
        my_file_pvert5 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(5, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(5, j, 0).imag()<< "\n";
        my_file_tvert3 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(3, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(3, j, 0).imag()<< "\n";

    }

    my_file_sigmaR.close();
    my_file_sigmaA.close();
    my_file_sigmaK.close();

    my_file_propR.close();
    my_file_propA.close();
    my_file_propK.close();

    my_file_avert1.close();
    my_file_pvert1.close();
    my_file_tvert1.close();

    my_file_avert3.close();
    my_file_pvert5.close();
    my_file_tvert3.close();
}


template<typename Q>
void SOPT(State<Q>& bare, double Lambda, State<Q> &state) {

    Propagator g = propag(Lambda, state.selfenergy, 'g', 'f');
    cout << "G calculated" << endl;


    cout << "bubble started" << endl;
    double t2 = get_time();
    //Lines 7-9
    double ta = get_time();
    Vertex<avert<comp> > dgammaa = a_bubble_function(bare.vertex, bare.vertex, g, '.');
    cout<<  "a - Bubble:";
    get_time(ta);

    double tp = get_time();
    Vertex<pvert<comp> > dgammap = p_bubble_function(bare.vertex, bare.vertex, g, '.');
    cout<<  "p - Bubble:";
    get_time(tp);

    double tt = get_time();
    Vertex<tvert<comp> > dgammat = t_bubble_function(bare.vertex, bare.vertex, g, '.');
    cout<<  "t - Bubble:";
    get_time(tt);

    cout << "bubble finished. ";
    get_time(t2);

    Propagator s = propag(Lambda, state.selfenergy, 's', 'f');

    //Line 41
    bare.vertex.densvertex.avertex = dgammaa.densvertex;
    bare.vertex.densvertex.pvertex = dgammap.densvertex;
    bare.vertex.densvertex.tvertex = dgammat.densvertex;

    double tloop = get_time();
    SelfEnergy<comp> Sigma_std = loop(bare.vertex, s);
    cout << "loop calculated ";
    get_time(tloop);

    Sigma_std*= dL;

    state.selfenergy += Sigma_std;

}

void writeOutSOPT(double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex)
{
    int i = fconv_Lambda(Lambda);

    ostringstream self_energyR, self_energyA, self_energyK, propR, propA, propK;
    ostringstream avert1, pvert1, tvert1, avert3, pvert5, tvert3;
    self_energyR << "self_energyR"<<i<<".dat";
    self_energyA << "self_energyA"<<i<<".dat";
    self_energyK << "self_energyK"<<i<<".dat";

    propR << "propagatorR"<<i<<".dat";
    propA << "propagatorA"<<i<<".dat";
    propK << "propagatorK"<<i<<".dat";

    avert1 << "avert1"<<i<<".dat";
    pvert1 << "pvert1"<<i<<".dat";
    tvert1 << "tvert1"<<i<<".dat";

    avert3 << "avert3"<<i<<".dat";
    pvert5 << "pvert5"<<i<<".dat";
    tvert3 << "tvert3"<<i<<".dat";

    ofstream my_file_sigmaR, my_file_sigmaA, my_file_sigmaK, my_file_propR, my_file_propA, my_file_propK;
    ofstream  my_file_avert1, my_file_pvert1, my_file_tvert1, my_file_avert3, my_file_pvert5, my_file_tvert3;
    my_file_sigmaR.open(self_energyR.str());
    my_file_sigmaA.open(self_energyA.str());
    my_file_sigmaK.open(self_energyK.str());

    my_file_propR.open(propR.str());
    my_file_propA.open(propA.str());
    my_file_propK.open(propK.str());

    my_file_avert1.open(avert1.str());
    my_file_pvert1.open(pvert1.str());
    my_file_tvert1.open(tvert1.str());

    my_file_avert3.open(avert3.str());
    my_file_pvert5.open(pvert5.str());
    my_file_tvert3.open(tvert3.str());

    for (int j = 0; j < ffreqs.size(); j++) {
        my_file_sigmaR << ffreqs[j] << " " << selfEnergy.sval(0, j).real() << " " <<  selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaA << ffreqs[j] << " " << selfEnergy.sval(0, j).real() << " " << -selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaK << ffreqs[j] << " " << selfEnergy.sval(1, j).real() << " " <<  selfEnergy.sval(1, j).imag() << "\n";

        my_file_propR << ffreqs[j] << " " << propagator.pval(0, j).real() << " " <<  propagator.pval(0, j).imag() << "\n";
        my_file_propA << ffreqs[j] << " " << propagator.pval(0, j).real() << " " << -propagator.pval(0, j).imag() << "\n";
        my_file_propK << ffreqs[j] << " " << propagator.pval(1, j).real() << " " <<  propagator.pval(1, j).imag() << "\n";
    }

    for (int j = 0; j<bfreqs.size(); j++){
        my_file_avert1 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_pvert1 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_tvert1 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).imag()<< "\n";

        my_file_avert3 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(3, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(3, j, 0).imag()<< "\n";
        my_file_pvert5 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(5, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(5, j, 0).imag()<< "\n";
        my_file_tvert3 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(3, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(3, j, 0).imag()<< "\n";

    }

    my_file_sigmaR.close();
    my_file_sigmaA.close();
    my_file_sigmaK.close();

    my_file_propR.close();
    my_file_propA.close();
    my_file_propK.close();

    my_file_avert1.close();
    my_file_pvert1.close();
    my_file_tvert1.close();

    my_file_avert3.close();
    my_file_pvert5.close();
    my_file_tvert3.close();
}