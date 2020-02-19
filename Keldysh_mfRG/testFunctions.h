//
// Created by Sa.Aguirre on 2/19/20.
//

#include "loop.h"

#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

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


    SelfEnergy<comp> Sigma_std = loop(dPsi.vertex, g);
    cout << "loop calculated ";

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

void testBubbles(Propagator& g1, Propagator& g2){
    vector<comp> w (nBOS);
    vector<comp> Pia_odd_even (nBOS);
    vector<comp> Pia_odd_odd (nBOS);

    for(auto w:bfreqs){
        int i = fconv_bos(w);
        IntegrandBubble integrandPia6  (g1, g2, false, w, 6,  'a');
        IntegrandBubble integrandPia9  (g1, g2, false, w, 9,  'a');
        IntegrandBubble integrandPia11 (g1, g2, false, w, 11, 'a');
        IntegrandBubble integrandPia13 (g1, g2, false, w, 13, 'a');
        IntegrandBubble integrandPia15 (g1, g2, false, w, 15, 'a');

        Pia_odd_even[i] = integrator(integrandPia11 , w_lower_b, w_upper_b) + integrator(integrandPia13 , w_lower_b, w_upper_b);
        Pia_odd_odd [i] = integrator(integrandPia6 , w_lower_b, w_upper_b) + integrator(integrandPia9 , w_lower_b, w_upper_b) + integrator(integrandPia15 , w_lower_b, w_upper_b);
    }


    ostringstream Pia_odd_evenfile;
    ostringstream Pia_odd_oddfile;

    Pia_odd_evenfile << "Output/Pia_odd_even.dat";
    Pia_odd_oddfile  << "Output/Pia_odd_odd.dat";


    ofstream my_file_Pia_odd_even ;
    ofstream my_file_Pia_odd_odd;

    my_file_Pia_odd_even.open(Pia_odd_evenfile.str()) ;
    my_file_Pia_odd_odd.open( Pia_odd_oddfile.str());

    for(int i = 0; i<nBOS; i++){
        my_file_Pia_odd_even<< bfreqs[i] << " " << Pia_odd_even[i].real() << " " << Pia_odd_even [i].imag() << "\n";
        my_file_Pia_odd_odd << bfreqs[i] << " " << Pia_odd_odd[i].real()  << " " << Pia_odd_odd[i].imag() << "\n";
    }

    my_file_Pia_odd_even.close() ;
    my_file_Pia_odd_odd.close();
}

void testSelfEnergy(Propagator& g1, State<comp>& state){

    sopt(state, 1.0, state);

    ostringstream selfEnergyR_cxx;
    ostringstream selfEnergyK_cxx;

    selfEnergyR_cxx << "Output/SelfEnergyR_cxx.dat";
    selfEnergyK_cxx << "Output/SelfEnergyK_cxx.dat";

    ofstream my_file_SelfEnergyR_cxx;
    ofstream my_file_SelfEnergyK_cxx;

    my_file_SelfEnergyR_cxx.open(selfEnergyR_cxx.str());
    my_file_SelfEnergyK_cxx.open(selfEnergyK_cxx.str());

    for(int i = 0; i<nFER; i++){
        my_file_SelfEnergyR_cxx<< bfreqs[i] << " " << state.selfenergy.sval(0,i).real() << " " << state.selfenergy.sval(0,i).imag() << "\n";
        my_file_SelfEnergyK_cxx<< bfreqs[i] << " " << state.selfenergy.sval(1,i).real() << " " << state.selfenergy.sval(1,i).imag() << "\n";
    }

    my_file_SelfEnergyR_cxx.close();
    my_file_SelfEnergyK_cxx.close();
}

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
