//
// Created by Sa.Aguirre on 2/19/20.
//

#include "loop.h"

#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

template <typename Q, typename T >
class IntegrandSigma{
    Vertex<T>& vertex;
    int iK;
    double w;
    int i_in;

    Propagator& propagator;
    int prop_iK;

public:
    IntegrandSigma(Vertex<T>& vertex_in, Propagator& prop_in, int iK_in, double w_in, int i_in_in, int prop_iK_in)
            :         vertex(vertex_in), propagator(prop_in), iK(iK_in),     w(w_in), i_in(i_in_in), prop_iK(prop_iK_in) {};

    auto operator()(double wp) -> Q
    {
//        Q aid1 = propagator.pvalsmooth(0, wp);
//        Q aid2 = conj(propagator.pvalsmooth(0, wp));
//        Q aid3 = propagator.pvalsmooth(1, wp);
//
//        Q termR = 2.*vertex.spinvertex.value(3, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(3, w, wp, w, i_in, 1, 'f');
//        Q termA = 2.*vertex.spinvertex.value(6, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(6, w, wp, w, i_in, 1, 'f');
//        Q termK = 2.*vertex.spinvertex.value(7, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(7, w, wp, w, i_in, 1, 'f');
        Q propTerm;
        if(prop_iK==-1)
            propTerm = conj(propagator.pvalsmooth(0, wp));
        else
            propTerm = propagator.pvalsmooth(prop_iK, wp);


        Q vertexTerm = 2.*vertex.spinvertex.value(iK, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(iK, w, wp, w, i_in, 1, 'f');

        return vertexTerm*propTerm;
    }
};

template <typename Q>
class IntegrandSelfEnergy_test{
    Propagator& propagator;
    int prop_iK;

public:
    IntegrandSelfEnergy_test(Propagator& prop_in, int prop_iK_in)
                :         propagator(prop_in), prop_iK(prop_iK_in) {};

    auto operator()(double wp) -> Q
    {
//        Q aid1 = propagator.pvalsmooth(0, wp);
//        Q aid2 = conj(propagator.pvalsmooth(0, wp));
//        Q aid3 = propagator.pvalsmooth(1, wp);
//
//        Q termR = 2.*vertex.spinvertex.value(3, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(3, w, wp, w, i_in, 1, 'f');
//        Q termA = 2.*vertex.spinvertex.value(6, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(6, w, wp, w, i_in, 1, 'f');
//        Q termK = 2.*vertex.spinvertex.value(7, w, wp, w, i_in, 0, 'f') + vertex.spinvertex.value(7, w, wp, w, i_in, 1, 'f');
        Q propTerm;
        if(prop_iK==-1)
            propTerm = conj(propagator.pvalsmooth(0, wp));
        else
            propTerm = propagator.pvalsmooth(prop_iK, wp);

        return -(1./2.)*(1./(2.*pi*im_unit))*propTerm;
    }
};

template <typename Q>
void loop_test(SelfEnergy<Q>& resp, Vertex<fullvert<Q> >& fullvertex, Propagator& prop, bool corr)
{
    vector<int> indicesSigmaR({3,6,7});
    vector<int> indicesSigmaK({1,4,5});


    for(auto iK : indicesSigmaR) {
//#pragma omp parallel for
        for (int iSE = 0; iSE < nSE * n_in; ++iSE) {
            int i = iSE / n_in;
            int i_in = iSE - i * n_in;

            double w = ffreqs[i];
            int prop_iK;
            switch (iK){
                case 3:
                    prop_iK = 0;
                    break;
                case 6:
                    prop_iK = -1;
                    break;
                case 7:
                    prop_iK = 1;
                    break;
                default:
                    prop_iK = 1000;

            }

            IntegrandSigma <Q, fullvert<Q>> integrandSigma(fullvertex, prop, iK, w, i_in, prop_iK);

            comp integratedSigma = 1./(2. * pi * im_unit)*( integrator(integrandSigma, w_lower_f, w_upper_f));
            if(corr) {
                Q vert = 2. * fullvertex.spinvertex.value(iK, w, w_upper_f, w, i_in, 0, 'f') + fullvertex.spinvertex.value(iK, w, w_upper_f, w, i_in, 1, 'f');
                integratedSigma += 1. / (2. * pi * im_unit) * vert * correctionFunctionSelfEnergy(prop_iK, w_upper_f, w_upper_f);
            }

            resp.addself(0, i, integratedSigma);
        }
    }

    for(auto iK : indicesSigmaK) {
//#pragma omp parallel for
        for (int iSE = 0; iSE < nSE * n_in; ++iSE) {
            int i = iSE / n_in;
            int i_in = iSE - i * n_in;

            double w = ffreqs[i];
            int prop_iK;
            switch (iK){
                case 1:
                    prop_iK = 0;
                    break;
                case 4:
                    prop_iK = -1;
                    break;
                case 5:
                    prop_iK = 1;
                    break;
                default:

                    prop_iK = 1000;

            }

            IntegrandSigma <Q, fullvert<Q>> integrandSigma(fullvertex, prop, iK, w, i_in, prop_iK);

            comp integratedSigma = 1./(2. * pi * im_unit)*( integrator(integrandSigma, w_lower_f, w_upper_f));
            if(corr) {
                Q vert = 2.*fullvertex.spinvertex.value(iK, w, w_upper_f, w, i_in, 0, 'f') + fullvertex.spinvertex.value(iK, w, w_upper_f, w, i_in, 1, 'f');
                integratedSigma += 1./(2. * pi * im_unit) * vert * correctionFunctionSelfEnergy(prop_iK, w_upper_f, w_upper_f);
            }
            resp.addself(1, i, integratedSigma);
        }
    }
}

template<typename Q>
void sopt_state(State<Q>& dPsi, double Lambda, State<Q> &state) {

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

//    double tp = get_time();
//    bubble_function(dPsi.vertex, state.vertex, state.vertex, g, g, 'p', diff, '.');
//    cout<<  "p - Bubble:";
//    get_time(tp);
//
//    double tt = get_time();
//    bubble_function(dPsi.vertex, state.vertex, state.vertex, g, g, 't', diff, '.');
//    cout<<  "t - Bubble:";
//    get_time(tt);

    cout << "bubble finished. ";
    get_time(t2);

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

void testBubbles(Propagator& g1, Propagator& g2, State<comp>& state){

    vector<comp> Pia_odd_even (nBOS);
    vector<comp> Pia_odd_odd (nBOS);
    vector<comp> Pia_odd_odd_c (nBOS);

    vector<comp> PiaOE(nBOS);
    vector<comp> PiaOO(nBOS);

    sopt_state(state, 1.0, state);

    for(int i=0; i<nBOS; i++){
        double w = bfreqs[i];

        tvert<comp>* aid = &state.vertex.spinvertex.tvertex;

        IntegrandBubble integrandPia6  (g1, g2, false, w, 6,  'a');     //AR
        IntegrandBubble integrandPia9  (g1, g2, false, w, 9,  'a');     //RA
        IntegrandBubble integrandPia11 (g1, g2, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13 (g1, g2, false, w, 13, 'a');     //RK
        IntegrandBubble integrandPia15 (g1, g2, false, w, 15, 'a');     //KK

        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        auto cont6  = integrator(integrandPia6 , w_lower_b, w_upper_b);
        auto cont9  = integrator(integrandPia9 , w_lower_b, w_upper_b);
        auto cont15 = integrator(integrandPia15 , w_lower_b, w_upper_b);

        Pia_odd_even[i] = 1./2.*(cont11 + cont13);
        Pia_odd_odd [i] = 1./2.*(cont6 + cont9 + cont15);

        Pia_odd_odd_c[i] = Pia_odd_odd[i] + 1./2.*(1./(2.*pi*im_unit)*(correctionFunctionBubbleAT(w, -1., 1., w_upper_b, w_upper_b)+correctionFunctionBubbleAT(w, 1.,-1., w_upper_b, w_upper_b)));

//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(1, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(1, w, 0, 1, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(7, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(7, w, 0, 1, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(8, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(8, w, 0, 1, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(14, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(14, w, 0, 1, *aid);

//        comp test_int = state.vertex.spinvertex.avertex.K1_vvalsmooth(1, w, 0, *aid);

        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);

//        comp test_direct = state.vertex.spinvertex.avertex.K1_vval(0, i, 0);


//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(3, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(3, w, 0, 1, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(5, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(5, w, 0, 1, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(10, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(10, w, 0, 1, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(12, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(12, w, 0, 1, *aid);

        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);

    }


    ostringstream Pia_odd_evenfile;
    ostringstream Pia_odd_oddfile;
    ostringstream Pia_odd_oddcfile;

    ostringstream PiaOEfile;
    ostringstream PiaOOfile;

    Pia_odd_evenfile << "Output/Pia_odd_even.dat";
    Pia_odd_oddfile  << "Output/Pia_odd_odd.dat";
    Pia_odd_oddcfile  << "Output/Pia_odd_odd_c.dat";

    PiaOEfile << "Output/PiaOE.dat";
    PiaOOfile << "Output/PiaOO.dat";

    ofstream my_file_Pia_odd_even ;
    ofstream my_file_Pia_odd_odd;
    ofstream my_file_Pia_odd_odd_c;

    ofstream my_file_PiaOE;
    ofstream my_file_PiaOO;

    my_file_Pia_odd_even.open(Pia_odd_evenfile.str()) ;
    my_file_Pia_odd_odd.open( Pia_odd_oddfile.str());
    my_file_Pia_odd_odd_c.open( Pia_odd_oddcfile.str());

    my_file_PiaOE.open(PiaOEfile.str());
    my_file_PiaOO.open(PiaOOfile.str());

    for(int i = 0; i<nBOS; i++){
        my_file_Pia_odd_even<< bfreqs[i] << " " << Pia_odd_even[i].real() << " " << Pia_odd_even [i].imag() << "\n";
        my_file_Pia_odd_odd << bfreqs[i] << " " << Pia_odd_odd[i].real()  << " " << Pia_odd_odd[i].imag() << "\n";
        my_file_Pia_odd_odd_c << bfreqs[i] << " " << Pia_odd_odd_c[i].real()  << " " << Pia_odd_odd_c[i].imag() << "\n";

        my_file_PiaOE << bfreqs[i] << " " << PiaOE[i].real() << " " << PiaOE[i].imag() << "\n";
        my_file_PiaOO << bfreqs[i] << " " << PiaOO[i].real() << " " << PiaOO[i].imag() << "\n";
    }

    my_file_Pia_odd_even.close() ;
    my_file_Pia_odd_odd.close();
    my_file_Pia_odd_odd_c.close();

    my_file_PiaOE.close();
    my_file_PiaOO.close();
}

void testSelfEnergy(Propagator& g1, State<comp>& state){

    sopt_state(state, 1.0, state);



//    vector<comp> Pia_odd_even (nBOS);
//    vector<comp> Pia_odd_odd (nBOS);
//    vector<comp> Pia_odd_odd_c (nBOS);
//
//    vector<comp> Pi6 (nBOS);
//    vector<comp> Pi9 (nBOS);
//    vector<comp> Pi11(nBOS);
//    vector<comp> Pi13(nBOS);
//    vector<comp> Pi15(nBOS);
//    vector<comp> correction(nBOS);
//
//    for(int i=0; i<nBOS; i++){
//        double w = bfreqs[i];
//
//        IntegrandBubble integrandPia6  (g1, g1, false, w, 6,  'a');     //AR
//        IntegrandBubble integrandPia9  (g1, g1, false, w, 9,  'a');     //RA
//        IntegrandBubble integrandPia11 (g1, g1, false, w, 11, 'a');     //KA
//        IntegrandBubble integrandPia13 (g1, g1, false, w, 13, 'a');     //RK
//        IntegrandBubble integrandPia15 (g1, g1, false, w, 15, 'a');     //KK
//
//        Pi11[i] = integrator(integrandPia11, w_lower_b, w_upper_b);
//        Pi13[i] = integrator(integrandPia13, w_lower_b, w_upper_b);
//
//        Pi6[i]  = integrator(integrandPia6 , w_lower_b, w_upper_b);
//        Pi9[i]  = integrator(integrandPia9 , w_lower_b, w_upper_b);
//        Pi15[i] = integrator(integrandPia15 , w_lower_b, w_upper_b);
//
//
//        Pia_odd_even[i] = 1./2.*(Pi11[i] + Pi13[i]);
//        Pia_odd_odd[i] = 1./2.*(Pi6[i] + Pi9[i] + Pi15[i]);
//
//        correction[i] = 1./2.*(1./(2.*pi*im_unit)*(correctionFunctionBubbleAT(w, -1., 1., w_upper_b, w_upper_b)+correctionFunctionBubbleAT(w, 1.,-1., w_upper_b, w_upper_b)));
//
//        Pia_odd_odd_c[i] = Pia_odd_odd[i] + correction[i];
//
//    }

    SelfEnergy<comp> selfie;
    SelfEnergy<comp> corrected;

    loop_test(selfie, state.vertex, g1, false);
    loop_test(corrected, state.vertex, g1, true);

//    for(int i=0; i<nFER; i++){
//        double w = ffreqs[i];
//
//        IntegrandSelfEnergy_test<comp> keldysh(g1, 1);
//        IntegrandSelfEnergy_test<comp> retarded(g1, 0);
//        IntegrandSelfEnergy_test<comp> advanced(g1, -1);
//
//        selfie.addself(0, i, (Pia_odd_odd_c[i])*integrator(retarded, w_lower_f, w_upper_f));
////        corrected.addself(0,i, selfie.sval(0, i) + correctionFunctionSelfEnergy(0, w_upper_f, w_upper_f));
//
//
//        selfie.addself(0, i, (Pia_odd_even[i])*integrator(keldysh, w_lower_f, w_upper_f));
////        corrected.addself(0,i, selfie.sval(0, i) + correctionFunctionSelfEnergy(1, w_upper_f, w_upper_f));
//
//    }




    ostringstream selfEnergyR_cxx;
    ostringstream selfEnergyK_cxx;
    ostringstream selfEnergyR_c_cxx;
    ostringstream selfEnergyK_c_cxx;

    selfEnergyR_cxx << "Output/SelfEnergyR_cxx.dat";
    selfEnergyK_cxx << "Output/SelfEnergyK_cxx.dat";
    selfEnergyR_c_cxx << "Output/SelfEnergyR_c_cxx.dat";
    selfEnergyK_c_cxx << "Output/SelfEnergyK_c_cxx.dat";

    ofstream my_file_SelfEnergyR_cxx;
    ofstream my_file_SelfEnergyK_cxx;
    ofstream my_file_SelfEnergyR_c_cxx;
    ofstream my_file_SelfEnergyK_c_cxx;

    my_file_SelfEnergyR_cxx.open(selfEnergyR_cxx.str());
    my_file_SelfEnergyK_cxx.open(selfEnergyK_cxx.str());
    my_file_SelfEnergyR_c_cxx.open(selfEnergyR_c_cxx.str());
    my_file_SelfEnergyK_c_cxx.open(selfEnergyK_c_cxx.str());

    for(int i = 0; i<nFER; i++){
        my_file_SelfEnergyR_cxx<< ffreqs[i] << " " << selfie.sval(0,i).real() << " " << selfie.sval(0,i).imag() << "\n";
        my_file_SelfEnergyK_cxx<< ffreqs[i] << " " << selfie.sval(1,i).real() << " " << selfie.sval(1,i).imag() << "\n";
        my_file_SelfEnergyR_c_cxx<< ffreqs[i] << " " << corrected.sval(0,i).real() << " " << corrected.sval(0,i).imag() << "\n";
        my_file_SelfEnergyK_c_cxx<< ffreqs[i] << " " << corrected.sval(1,i).real() << " " << corrected.sval(1,i).imag() << "\n";
    }

    my_file_SelfEnergyR_cxx.close();
    my_file_SelfEnergyK_cxx.close();
    my_file_SelfEnergyR_c_cxx.close();
    my_file_SelfEnergyK_c_cxx.close();
}

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
