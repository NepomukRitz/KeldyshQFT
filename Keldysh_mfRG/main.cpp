#include <cstdlib>
#include <bits/stdc++.h> // TODO: what is this needed for?
#include <iostream>
#include <fstream>
#include <complex>

#include "parameters.h"
#include "vertex.h"
#include "state.h"
#include "loop.h"
//#include "a_bubble.h"
//#include "p_bubble.h"
//#include "t_bubble.h"
#include "bubbles.h"
#include "propagator.h"
#include "selfenergy.h"
#include "hdf5_routines.h"

#include "testFunctions.h"


using namespace std;

typedef complex<double> comp;  // TODO: complex double as two doubles??

void writeOutFile(double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex);
void flow();
template <typename Q> void setInitialConditions(State<Q>& state);
template <typename Q> void derivative(State<Q>& dPsi, double Lambda, State<Q>& state);
template <typename Q> void RungeKutta4thOrder(State<Q>& dPsi, double Lambda, State<Q>& state);

template <typename Q> void sopt(State<Q>& bare, double Lambda, State<Q>& state);
void writeOutSOPT(double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex);
void writePropagators(Propagator& free, Propagator& full);

auto main() -> int {

    setUpBosGrid();
    setUpFerGrid();
    setUpFlowGrid();

    MPI_Init(NULL, NULL);


#ifndef SOPT // standard fRG flow
    double t0 = get_time();
    print("Start of flow", true);
    flow();
    print("Total execution time: ");
    get_time(t0);

#else // SOPT Code here:
    SelfEnergy<comp> diffZero;

    for(int i=0; i<nEVO; ++i) {
        double Lambda = flow_grid[i];

        State<comp> bare (Lambda);
        for (auto j:odd_Keldysh) {
            bare.vertex.spinvertex.irred.setvert(j, 0, 0.5*U);
        }

        sopt(bare, Lambda, bare);

        Propagator control = propag(Lambda, bare.selfenergy, diffZero, 'g');
        writeOutSOPT(Lambda, control, bare.selfenergy, bare.vertex);
    }
#endif

//    SelfEnergy<comp> self;
//    for (int i = 0; i < nSE; ++i) {
//        self.setself(0, i, U/2.);
//        self.setself(1, i, 0.);
//    }
//    SelfEnergy<comp> zero;
//    Propagator bubbles_prop = propag(1.0, self, zero, 'g');
//    testBubbles(bubbles_prop, bubbles_prop);

    MPI_Finalize();

    return 0;
}



template <typename Q>
void derivative(State<Q>& dPsi, double Lambda, State<Q>& state) {
    /*Here I begin implementing Fabian's pseudocode*/
    //Line 1
    Propagator S = propag(Lambda, state.selfenergy, state.selfenergy, 's');
    print("S calculated", true);
    //Line 2
    Propagator G = propag(Lambda, state.selfenergy, state.selfenergy, 'g');
    print("G calculated", true);
    //Line 3
    SelfEnergy<comp> Sigma_std = loop(state.vertex, S);
    print("loop calculated", true);
    //Line 4
    dPsi.selfenergy = Sigma_std;
    //Line 6
    Propagator extension = propag(Lambda, state.selfenergy, dPsi.selfenergy, 'e');\
    Propagator dG = S + extension;


    print("diff bubble started", true);
    bool diff = true;
    double t2 = get_time();
    //Lines 7-9
    double ta = get_time();
    bubble_function(dPsi.vertex, state.vertex, state.vertex, G, dG, 'a', diff, '.');
    print("a - Bubble:");
    get_time(ta);

    double tp = get_time();
    bubble_function(dPsi.vertex, state.vertex, state.vertex, G, dG, 'p', diff, '.');
    print("p - Bubble:");
    get_time(tp);

    double tt = get_time();
    bubble_function(dPsi.vertex, state.vertex, state.vertex, G, dG, 't', diff, '.');
    print("t - Bubble:");
    get_time(tt);

    print("diff bubble finished. ");
    get_time(t2);

#ifdef SOPT
    Propagator s = propag(Lambda, state.selfenergy, state.selfenergy, 's');
    dPsi.selfenergy = loop(state.vertex, s);
#endif
#ifdef NLOOPS
    #if NLOOPS > 1
//    Lines 10-13   => Multi-loop
    double t4 = get_time();
    /*Create two new vertices to accommodate the contributions on each side */
    Vertex<fullvert<Q> > dGammaL;
    Vertex<fullvert<Q> > dGammaR;
    //Change from differentiated to regular bubbles
    diff = false;

    bubble_function(dGammaL, dPsi.vertex, state.vertex, G, G, 'a', diff, 'L');
    bubble_function(dGammaL, dPsi.vertex, state.vertex, G, G, 'p', diff, 'L');
    bubble_function(dGammaL, dPsi.vertex, state.vertex, G, G, 't', diff, 'L');

    bubble_function(dGammaR, state.vertex, dPsi.vertex, G, G, 'a', diff, 'R');
    bubble_function(dGammaR, state.vertex, dPsi.vertex, G, G, 'p', diff, 'R');
    bubble_function(dGammaR, state.vertex, dPsi.vertex, G, G, 't', diff, 'R');

    print("Bubbles calculated: ", true);
    get_time(t4);


    //Lines 14-17
    Vertex<fullvert<Q> > dGammaT = dGammaL + dGammaR;
    dPsi.vertex += dGammaT;
    print("2-loops done. \n");

    //Lines 18-33
    #if NLOOPS >=3
    Vertex<fullvert<Q> > dGammaC;
    Vertex<fullvert<Q> > dGammaCtb;
    for(int i=3; i<=NLOOPS; i++){
        bubble_function(dGammaC, state.vertex, dGammaL, G, G, 'a', diff, 'C');
        bubble_function(dGammaC, state.vertex, dGammaL, G, G, 'p', diff, 'C');

        //This corresponds to Line 29 in the pseudo-code and is important for self-energy corrections.
        dGammaCtb +=dGammaC;

        bubble_function(dGammaC, state.vertex, dGammaL, G, G, 't', diff, 'C');

        bubble_function(dGammaL, dGammaT, state.vertex, G, G, 'a', diff, 'L');
        bubble_function(dGammaL, dGammaT, state.vertex, G, G, 'p', diff, 'L');
        bubble_function(dGammaL, dGammaT, state.vertex, G, G, 't', diff, 'L');

        bubble_function(dGammaR, state.vertex, dGammaT, G, G, 'a', diff, 'R');
        bubble_function(dGammaR, state.vertex, dGammaT, G, G, 'p', diff, 'R');
        bubble_function(dGammaR, state.vertex, dGammaT, G, G, 't', diff, 'R');

        dGammaT = dGammaL + dGammaC + dGammaR;
        dPsi.vertex += dGammaT;

//        if(max_r(norm(dGammaT)/norm(dPsi.vertex)) < tol_vertex){ //TODO define a sensible norm for the vertices and a good way to implement this condition
//            break;
//        }
        printf("%i-loops done. \n", i);
    }
    #endif
    #endif
#endif

//#if PROP_TYPE==2
//    //Lines 33-41
//    SelfEnergy<comp> dSigma_tbar = loop(dGammaCtb, G);
//    Propagator corrected = propag(Lambda, dPsi.selfenergy, dSigma_tbar, 'e', false);
//    SelfEnergy<comp> dSigma_t = loop(dPsi.vertex, corrected);
//    dPsi.selfenergy += (dSigma_t + dSigma_tbar);
//#endif


    double t_multiply = get_time();
    dPsi *= dL;
    print("dPsi multiplied. ");
    get_time(t_multiply);
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

template <typename Q> void setInitialConditions (State<Q>& state){
    //Initial conditions
    for (int i = 0; i < nSE; ++i) {
        state.selfenergy.setself(0, i, U/2.);
        state.selfenergy.setself(1, i, 0.);
    }
    print("self energy and diff self energy assigned", true);

    for (auto i:odd_Keldysh) {
        state.vertex.densvertex.irred.setvert(i, 0, 0.);
        state.vertex.spinvertex.irred.setvert(i, 0, U/2.);
    }
    print("vertex assigned", true);
}

void flow(){


    int world_rank = mpi_world_rank();

    double  t0 = get_time();
    State<comp> state(Lambda_ini);
    print("State created. ");
    get_time(t0);

//    vector<double> Lambdas(10);
    const H5std_string FILE_NAME("testfile.h5");
//    write_hdf(FILE_NAME,10.,10,state);
//    add_hdf(FILE_NAME,1,10,state,Lambdas);
//    State<comp> out= read_hdf<comp>(FILE_NAME,1,10,Lambdas);

    //Set initial conditions
    setInitialConditions(state);

    //Create objects to print at zero-th evolution step
    SelfEnergy<comp> diffZero = SelfEnergy<comp>();
    Propagator initial = propag(state.Lambda, state.selfenergy, diffZero, 's');

    if (world_rank == 0){
        write_hdf(FILE_NAME, 0, nEVO, state);
        writeOutFile(state.Lambda, initial, state.selfenergy, state.vertex);
    }

    for(int i=1; i<nEVO; ++i) {
        double Lambda = flow_grid[i-1];
        double tder = get_time();
        state.Lambda = Lambda;

        State<comp> dPsi(Lambda);

        derivative(dPsi, Lambda, state);
//        RungeKutta4thOrder(dPsi, Lambda, state);


        double tadd = get_time();
        state += dPsi;
        print("Added. ");
        get_time(tadd);

        double next_Lambda = flow_grid[i];
        //state.Lambda = next_Lambda;

        Propagator control = propag(state.Lambda, state.selfenergy, diffZero, 's');


        if (world_rank == 0){
            add_hdf(FILE_NAME, i, nEVO, state, flow_grid);
            writeOutFile(next_Lambda, control, state.selfenergy, state.vertex);
//            test_hdf5(FILE_NAME, i, state);
        }

        print("One RK-derivative step. ");
        get_time(tder);
    }



}

//Writes of .dat files of the propagator, the self energy and selected values of the vertex for different values of Lambda during the fRG flow
void writeOutFile(double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex)
{

    double t_write = get_time();
    int i = fconv_Lambda(Lambda);

    ostringstream self_energyR, self_energyA, self_energyK, propR, propA, propK;
    ofstream my_file_sigmaR, my_file_sigmaA, my_file_sigmaK, my_file_propR, my_file_propA, my_file_propK;

    self_energyR << "Output/self_energyR"<<i<<".dat";
    self_energyA << "Output/self_energyA"<<i<<".dat";
    self_energyK << "Output/self_energyK"<<i<<".dat";
    propR << "Output/propagatorR"<<i<<".dat";
    propA << "Output/propagatorA"<<i<<".dat";
    propK << "Output/propagatorK"<<i<<".dat";

    my_file_sigmaR.open(self_energyR.str());
    my_file_sigmaA.open(self_energyA.str());
    my_file_sigmaK.open(self_energyK.str());
    my_file_propR.open(propR.str());
    my_file_propA.open(propA.str());
    my_file_propK.open(propK.str());

    for (int j = 0; j < ffreqs.size(); j++) {
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

#if DIAG_CLASS >=1
    ostringstream avert1K1, pvert1K1, tvert1K1, avert3K1, pvert5K1, tvert3K1;
    avert1K1 << "Output/avert1K1"<<i<<".dat";
    pvert1K1 << "Output/pvert1K1"<<i<<".dat";
    tvert1K1 << "Output/tvert1K1"<<i<<".dat";
    avert3K1 << "Output/avert3K1"<<i<<".dat";
    pvert5K1 << "Output/pvert5K1"<<i<<".dat";
    tvert3K1 << "Output/tvert3K1"<<i<<".dat";

    ofstream  my_file_avert1K1, my_file_pvert1K1, my_file_tvert1K1, my_file_avert3K1, my_file_pvert5K1, my_file_tvert3K1;

    my_file_avert1K1.open(avert1K1.str());
    my_file_pvert1K1.open(pvert1K1.str());
    my_file_tvert1K1.open(tvert1K1.str());
    my_file_avert3K1.open(avert3K1.str());
    my_file_pvert5K1.open(pvert5K1.str());
    my_file_tvert3K1.open(tvert3K1.str());

    for (int j = 0; j<bfreqs.size(); j++){
        my_file_avert1K1 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(0, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(0, j, 0).imag()<< "\n";
        my_file_pvert1K1 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(0, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(0, j, 0).imag()<< "\n";
        my_file_tvert1K1 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(0, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(0, j, 0).imag()<< "\n";
        my_file_avert3K1 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_pvert5K1 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_tvert3K1 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).imag()<< "\n";
    }

    my_file_avert1K1.close();
    my_file_pvert1K1.close();
    my_file_tvert1K1.close();
    my_file_avert3K1.close();
    my_file_pvert5K1.close();
    my_file_tvert3K1.close();
#endif

#if DIAG_CLASS>=2
    ostringstream avert0K2, avert1K2, avert2K2, avert3K2, avert11K2;
    ostringstream pvert0K2, pvert1K2, pvert4K2, pvert5K2, pvert13K2;
    ostringstream tvert0K2, tvert1K2, tvert2K2, tvert3K2, tvert11K2;


    avert0K2 << "Output/avert0K2"<<i<<".dat";
    avert1K2 << "Output/avert1K2"<<i<<".dat";
    avert2K2 << "Output/avert2K2"<<i<<".dat";
    avert3K2 << "Output/avert3K2"<<i<<".dat";
    avert11K2<< "Output/avert11K2"<<i<<".dat";

    pvert0K2 << "Output/pvert0K2"<<i<<".dat";
    pvert1K2 << "Output/pvert1K2"<<i<<".dat";
    pvert4K2 << "Output/pvert4K2"<<i<<".dat";
    pvert5K2 << "Output/pvert5K2"<<i<<".dat";
    pvert13K2<< "Output/pvert13K2"<<i<<".dat";

    tvert0K2 << "Output/tvert0K2"<<i<<".dat";
    tvert1K2 << "Output/tvert1K2"<<i<<".dat";
    tvert2K2 << "Output/tvert2K2"<<i<<".dat";
    tvert3K2 << "Output/tvert3K2"<<i<<".dat";
    tvert11K2<< "Output/tvert11K2"<<i<<".dat";

    ofstream  my_file_avert0K2, my_file_avert1K2, my_file_avert2K2, my_file_avert3K2, my_file_avert11K2;
    ofstream  my_file_pvert0K2, my_file_pvert1K2, my_file_pvert4K2, my_file_pvert5K2, my_file_pvert13K2;
    ofstream  my_file_tvert0K2, my_file_tvert1K2, my_file_tvert2K2, my_file_tvert3K2, my_file_tvert11K2;

    my_file_avert0K2.open(avert0K2.str());
    my_file_avert1K2.open(avert1K2.str());
    my_file_avert2K2.open(avert2K2.str());
    my_file_avert3K2.open(avert3K2.str());
    my_file_avert11K2.open(avert11K2.str());

    my_file_pvert0K2.open(pvert0K2.str());
    my_file_pvert1K2.open(pvert1K2.str());
    my_file_pvert4K2.open(pvert4K2.str());
    my_file_pvert5K2.open(pvert5K2.str());
    my_file_pvert13K2.open(pvert13K2.str());

    my_file_tvert0K2.open(tvert0K2.str());
    my_file_tvert1K2.open(tvert1K2.str());
    my_file_tvert2K2.open(tvert2K2.str());
    my_file_tvert3K2.open(tvert3K2.str());
    my_file_tvert11K2.open(tvert11K2.str());


    for(int j=0; i<bfreqs.size(); ++i){
        for(int k=0; j<ffreqs.size(); j++){

            my_file_avert0K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K2_vval(0, j, k, 0).real() << " " << vertex.densvertex.avertex.K2_vval(0, j, k, 0).imag()<< "\n";
            my_file_avert1K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K2_vval(1, j, k, 0).real() << " " << vertex.densvertex.avertex.K2_vval(1, j, k, 0).imag()<< "\n";
            my_file_avert2K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K2_vval(2, j, k, 0).real() << " " << vertex.densvertex.avertex.K2_vval(2, j, k, 0).imag()<< "\n";
            my_file_avert3K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K2_vval(3, j, k, 0).real() << " " << vertex.densvertex.avertex.K2_vval(3, j, k, 0).imag()<< "\n";
            my_file_avert11K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K2_vval(4, j, k, 0).real() << " " << vertex.densvertex.avertex.K2_vval(4, j, k, 0).imag()<< "\n";

            my_file_pvert0K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K2_vval(0, j, k, 0).real() << " " << vertex.densvertex.pvertex.K2_vval(0, j, k, 0).imag()<< "\n";
            my_file_pvert1K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K2_vval(1, j, k, 0).real() << " " << vertex.densvertex.pvertex.K2_vval(1, j, k, 0).imag()<< "\n";
            my_file_pvert4K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K2_vval(2, j, k, 0).real() << " " << vertex.densvertex.pvertex.K2_vval(2, j, k, 0).imag()<< "\n";
            my_file_pvert5K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K2_vval(3, j, k, 0).real() << " " << vertex.densvertex.pvertex.K2_vval(3, j, k, 0).imag()<< "\n";
            my_file_pvert13K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K2_vval(4, j, k, 0).real() << " " << vertex.densvertex.pvertex.K2_vval(4, j, k, 0).imag()<< "\n";

            my_file_tvert0K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K2_vval(0, j, k, 0).real() << " " << vertex.densvertex.tvertex.K2_vval(0, j, k, 0).imag()<< "\n";
            my_file_tvert1K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K2_vval(1, j, k, 0).real() << " " << vertex.densvertex.tvertex.K2_vval(1, j, k, 0).imag()<< "\n";
            my_file_tvert2K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K2_vval(2, j, k, 0).real() << " " << vertex.densvertex.tvertex.K2_vval(2, j, k, 0).imag()<< "\n";
            my_file_tvert3K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K2_vval(3, j, k, 0).real() << " " << vertex.densvertex.tvertex.K2_vval(3, j, k, 0).imag()<< "\n";
            my_file_tvert11K2 << bfreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K2_vval(4, j, k, 0).real() << " " << vertex.densvertex.tvertex.K2_vval(4, j, k, 0).imag()<< "\n";


        }
    }

    my_file_avert0K2.close();
    my_file_avert1K2.close();
    my_file_avert2K2.close();
    my_file_avert3K2.close();
    my_file_avert11K2.close();

    my_file_pvert0K2.close();
    my_file_pvert1K2.close();
    my_file_pvert4K2.close();
    my_file_pvert5K2.close();
    my_file_pvert13K2.close();

    my_file_tvert0K2.close();
    my_file_tvert1K2.close();
    my_file_tvert2K2.close();
    my_file_tvert3K2.close();
    my_file_tvert11K2.close();
#endif

#if DIAG_CLASS>=3
    ostringstream avert0K3, avert1K3, avert3K3, avert5K3, avert6K3, avert7K3;
    ostringstream pvert0K3, pvert1K3, pvert3K3, pvert5K3, pvert6K3, pvert7K3;
    ostringstream tvert0K3, tvert1K3, tvert3K3, tvert5K3, tvert6K3, tvert7K3;


    avert0K3 << "avert0K3"<<i<<".dat";
    avert1K3 << "avert1K3"<<i<<".dat";
    avert3K3 << "avert3K3"<<i<<".dat";
    avert5K3 << "avert5K3"<<i<<".dat";
    avert6K3 << "avert6K3"<<i<<".dat";
    avert7K3 << "avert7K3"<<i<<".dat";

    pvert0K3 << "pvert0K3"<<i<<".dat";
    pvert1K3 << "pvert1K3"<<i<<".dat";
    pvert3K3 << "pvert3K3"<<i<<".dat";
    pvert5K3 << "pvert5K3"<<i<<".dat";
    pvert6K3 << "pvert6K3"<<i<<".dat";
    pvert7K3 << "pvert7K3"<<i<<".dat";

    tvert0K3 << "tvert0K3"<<i<<".dat";
    tvert1K3 << "tvert1K3"<<i<<".dat";
    tvert3K3 << "tvert3K3"<<i<<".dat";
    tvert5K3 << "tvert5K3"<<i<<".dat";
    tvert6K3 << "tvert6K3"<<i<<".dat";
    tvert7K3 << "tvert7K3"<<i<<".dat";

    ofstream  my_file_avert0K3, my_file_avert1K3, my_file_avert3K3, my_file_avert5K3, my_file_avert6K3, my_file_avert7K3;
    ofstream  my_file_pvert0K3, my_file_pvert1K3, my_file_pvert3K3, my_file_pvert5K3, my_file_pvert6K3, my_file_pvert7K3;
    ofstream  my_file_tvert0K3, my_file_tvert1K3, my_file_tvert3K3, my_file_tvert5K3, my_file_tvert6K3, my_file_tvert7K3;

    my_file_avert0K3.open(avert0K3.str());
    my_file_avert1K3.open(avert1K3.str());
    my_file_avert3K3.open(avert3K3.str());
    my_file_avert5K3.open(avert5K3.str());
    my_file_avert6K3.open(avert6K3.str());
    my_file_avert7K3.open(avert7K3.str());

    my_file_pvert0K3.open(pvert0K3.str());
    my_file_pvert1K3.open(pvert1K3.str());
    my_file_pvert3K3.open(pvert3K3.str());
    my_file_pvert5K3.open(pvert5K3.str());
    my_file_pvert6K3.open(pvert6K3.str());
    my_file_pvert7K3.open(pvert7K3.str());

    my_file_tvert0K3.open(tvert0K3.str());
    my_file_tvert1K3.open(tvert1K3.str());
    my_file_tvert3K3.open(tvert3K3.str());
    my_file_tvert5K3.open(tvert5K3.str());
    my_file_tvert6K3.open(tvert6K3.str());
    my_file_tvert7K3.open(tvert7K3.str());


    for(int j=0; i<ffreqs.size(); ++i){
        for(int k=0; j<ffreqs.size(); j++){

            int l = fconv_bos(0.);

            my_file_avert0K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K3_vval(0, l, j, k, 0).real() << " " << vertex.densvertex.avertex.K3_vval(0, l, j, k, 0).imag()<< "\n";
            my_file_avert1K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K3_vval(1, l, j, k, 0).real() << " " << vertex.densvertex.avertex.K3_vval(1, l, j, k, 0).imag()<< "\n";
            my_file_avert3K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K3_vval(2, l, j, k, 0).real() << " " << vertex.densvertex.avertex.K3_vval(2, l, j, k, 0).imag()<< "\n";
            my_file_avert5K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K3_vval(3, l, j, k, 0).real() << " " << vertex.densvertex.avertex.K3_vval(3, l, j, k, 0).imag()<< "\n";
            my_file_avert6K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K3_vval(4, l, j, k, 0).real() << " " << vertex.densvertex.avertex.K3_vval(4, l, j, k, 0).imag()<< "\n";
            my_file_avert7K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.avertex.K3_vval(5, l, j, k, 0).real() << " " << vertex.densvertex.avertex.K3_vval(5, l, j, k, 0).imag()<< "\n";

            my_file_pvert0K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K3_vval(0, l, j, k, 0).real() << " " << vertex.densvertex.pvertex.K3_vval(0, l, j, k, 0).imag()<< "\n";
            my_file_pvert1K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K3_vval(1, l, j, k, 0).real() << " " << vertex.densvertex.pvertex.K3_vval(1, l, j, k, 0).imag()<< "\n";
            my_file_pvert3K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K3_vval(2, l, j, k, 0).real() << " " << vertex.densvertex.pvertex.K3_vval(2, l, j, k, 0).imag()<< "\n";
            my_file_pvert5K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K3_vval(3, l, j, k, 0).real() << " " << vertex.densvertex.pvertex.K3_vval(3, l, j, k, 0).imag()<< "\n";
            my_file_pvert6K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K3_vval(4, l, j, k, 0).real() << " " << vertex.densvertex.pvertex.K3_vval(4, l, j, k, 0).imag()<< "\n";
            my_file_pvert7K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.pvertex.K3_vval(5, l, j, k, 0).real() << " " << vertex.densvertex.pvertex.K3_vval(5, l, j, k, 0).imag()<< "\n";

            my_file_tvert0K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K3_vval(0, l, j, k, 0).real() << " " << vertex.densvertex.tvertex.K3_vval(0, l, j, k, 0).imag()<< "\n";
            my_file_tvert1K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K3_vval(1, l, j, k, 0).real() << " " << vertex.densvertex.tvertex.K3_vval(1, l, j, k, 0).imag()<< "\n";
            my_file_tvert3K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K3_vval(2, l, j, k, 0).real() << " " << vertex.densvertex.tvertex.K3_vval(2, l, j, k, 0).imag()<< "\n";
            my_file_tvert5K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K3_vval(3, l, j, k, 0).real() << " " << vertex.densvertex.tvertex.K3_vval(3, l, j, k, 0).imag()<< "\n";
            my_file_tvert6K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K3_vval(4, l, j, k, 0).real() << " " << vertex.densvertex.tvertex.K3_vval(4, l, j, k, 0).imag()<< "\n";
            my_file_tvert7K3 << ffreqs[j] << " " << ffreqs[k] << " " << vertex.densvertex.tvertex.K3_vval(5, l, j, k, 0).real() << " " << vertex.densvertex.tvertex.K3_vval(5, l, j, k, 0).imag()<< "\n";
        }
    }

    my_file_avert0K3.close();
    my_file_avert1K3.close();
    my_file_avert3K3.close();
    my_file_avert5K3.close();
    my_file_avert6K3.close();
    my_file_avert7K3.close();

    my_file_pvert0K3.close();
    my_file_pvert1K3.close();
    my_file_pvert3K3.close();
    my_file_pvert5K3.close();
    my_file_pvert6K3.close();
    my_file_pvert7K3.close();

    my_file_tvert0K3.close();
    my_file_tvert1K3.close();
    my_file_tvert3K3.close();
    my_file_tvert5K3.close();
    my_file_tvert6K3.close();
    my_file_tvert7K3.close();
#endif

    cout << "Wrote out. ";
    get_time(t_write);
}


template<typename Q>
void sopt(State<Q>& dPsi, double Lambda, State<Q> &state) {

    Propagator g = propag(Lambda, state.selfenergy, state.selfenergy, 'g');
    cout << "G calculated" << endl;

    bool diff = false;

    cout << "bubble started" << endl;
    double t2 = get_time();
    //Lines 7-9
    double ta = get_time();
    bubble_function(dPsi.vertex, state.vertex, state.vertex, g, g, 'a', diff, '.');
    cout<<  "a - Bubble:";
    get_time(ta);

    double tp = get_time();
    bubble_function(dPsi.vertex, state.vertex, state.vertex, g, g, 'p', diff, '.');
    cout<<  "p - Bubble:";
    get_time(tp);

    double tt = get_time();
    bubble_function(dPsi.vertex, state.vertex, state.vertex, g, g, 't', diff, '.');
    cout<<  "t - Bubble:";
    get_time(tt);

    cout << "bubble finished. ";
    get_time(t2);

    Propagator s = propag(Lambda, state.selfenergy, state.selfenergy, 's');

    double tloop = get_time();
    SelfEnergy<comp> Sigma_std = loop(dPsi.vertex, s);
    cout << "loop calculated ";
    get_time(tloop);

    Sigma_std*= dL;

    state.selfenergy += Sigma_std;

}

void writeOutSOPT(double Lambda, Propagator& propagator, SelfEnergy<comp>& selfEnergy, Vertex<fullvert<comp> >& vertex)
{
    int i = fconv_Lambda(Lambda);

    ostringstream self_energyR, self_energyA, self_energyK, propR, propA, propK;
    self_energyR << "Output/self_energyR"<<i<<".dat";
    self_energyA << "Output/self_energyA"<<i<<".dat";
    self_energyK << "Output/self_energyK"<<i<<".dat";

    propR << "Output/propagatorR"<<i<<".dat";
    propA << "Output/propagatorA"<<i<<".dat";
    propK << "Output/propagatorK"<<i<<".dat";


    ofstream my_file_sigmaR, my_file_sigmaA, my_file_sigmaK, my_file_propR, my_file_propA, my_file_propK;
    my_file_sigmaR.open(self_energyR.str());
    my_file_sigmaA.open(self_energyA.str());
    my_file_sigmaK.open(self_energyK.str());

    my_file_propR.open(propR.str());
    my_file_propA.open(propA.str());
    my_file_propK.open(propK.str());


    for (int j = 0; j < ffreqs.size(); j++) {
        my_file_sigmaR << ffreqs[j] << " " << selfEnergy.sval(0, j).real() << " " <<  selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaA << ffreqs[j] << " " << selfEnergy.sval(0, j).real() << " " << -selfEnergy.sval(0, j).imag() << "\n";
        my_file_sigmaK << ffreqs[j] << " " << selfEnergy.sval(1, j).real() << " " <<  selfEnergy.sval(1, j).imag() << "\n";

        my_file_propR << ffreqs[j] << " " << propagator.pval(0, j).real() << " " <<  propagator.pval(0, j).imag() << "\n";
        my_file_propA << ffreqs[j] << " " << propagator.pval(0, j).real() << " " << -propagator.pval(0, j).imag() << "\n";
        my_file_propK << ffreqs[j] << " " << propagator.pval(1, j).real() << " " <<  propagator.pval(1, j).imag() << "\n";
    }

#if DIAG_CLASS >=1
    ostringstream avert1, pvert1, tvert1, avert3, pvert5, tvert3;

    avert1 << "Output/avert1"<<i<<".dat";
    pvert1 << "Output/pvert1"<<i<<".dat";
    tvert1 << "Output/tvert1"<<i<<".dat";

    avert3 << "Output/avert3"<<i<<".dat";
    pvert5 << "Output/pvert5"<<i<<".dat";
    tvert3 << "Output/tvert3"<<i<<".dat";

    ofstream  my_file_avert1, my_file_pvert1, my_file_tvert1, my_file_avert3, my_file_pvert5, my_file_tvert3;

    my_file_avert1.open(avert1.str());
    my_file_pvert1.open(pvert1.str());
    my_file_tvert1.open(tvert1.str());

    my_file_avert3.open(avert3.str());
    my_file_pvert5.open(pvert5.str());
    my_file_tvert3.open(tvert3.str());

    for (int j = 0; j<bfreqs.size(); j++){
        my_file_avert1 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(0, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(0, j, 0).imag()<< "\n";
        my_file_pvert1 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(0, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(0, j, 0).imag()<< "\n";
        my_file_tvert1 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(0, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(0, j, 0).imag()<< "\n";

        my_file_avert3 << bfreqs[j] << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.avertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_pvert5 << bfreqs[j] << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.pvertex.K1_vval(1, j, 0).imag()<< "\n";
        my_file_tvert3 << bfreqs[j] << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).real() << " " << vertex.densvertex.tvertex.K1_vval(1, j, 0).imag()<< "\n";

    }

    my_file_avert1.close();
    my_file_pvert1.close();
    my_file_tvert1.close();

    my_file_avert3.close();
    my_file_pvert5.close();
    my_file_tvert3.close();
#endif

    my_file_sigmaR.close();
    my_file_sigmaA.close();
    my_file_sigmaK.close();

    my_file_propR.close();
    my_file_propA.close();
    my_file_propK.close();
}