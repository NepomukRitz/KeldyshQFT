//
// Created by Sa.Aguirre on 2/19/20.
//

#include "loop.h"

#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

#include "write_data2file.h"

/**
 * Integrand for the SelfEnergy when terms are calculated individually
 * @tparam Q : Type of data, usually comp
 * @tparam T : Type of Vertex obejct, usually fullvert<Q>
 */
template <typename Q, typename T >
class IntegrandSigma{
    const Vertex<T>& vertex;
    int iK;
    double v;
    int i_in;

    const Propagator& propagator;
    int prop_iK;

public:
    /**
     * Constructor
     * @param vertex_in : Vertex object (ref) for the integrand. Only one component is read out
     * @param prop_in   : Propagator object (ref) for the integrand. Only one component is read out
     * @param iK_in     : Keldysh index for the Vertex to be read out
     * @param v_in      : Frequency at which the calculation is carried out
     * @param i_in_in   : Internal index
     * @param prop_iK_in: Keldysh index for the Propagator to be read out
     */
    IntegrandSigma(const Vertex<T>& vertex_in, const Propagator& prop_in, int iK_in, double v_in, int i_in_in, int prop_iK_in)
            :               vertex(vertex_in),       propagator(prop_in), iK(iK_in),     v(v_in), i_in(i_in_in), prop_iK(prop_iK_in) {};

    /**
     * Call operator
     * @param vp : Frecuency over which it's being integrated.
     * @return The value of the iK component of the vertex at (v, v', v) times the propagator component prop_iK at v'.
     *          If SOPT, the value of the vertex changes from (2V+V^) to just V.
     */
    auto operator()(double vp) const-> Q
    {
        Q aid, propTerm;
        if(prop_iK==-1) {
            propTerm = conj(propagator.valsmooth(0, vp));
        }
        else if(prop_iK==0) {
            propTerm = propagator.valsmooth(0, vp);
        }
        else {
            propTerm = propagator.valsmooth(1, vp);
        }

#ifdef FLOW
        Q vertexTerm = 2.*vertex.spinvertex.value(iK, v, vp, v, i_in, 0, 'f')+ vertex.spinvertex.value(iK, v, vp, v, i_in, 1, 'f');
#else
        Q vertexTerm = vertex.spinvertex.value(iK, v, vp, v, i_in, 0, 'f');
#endif

        return vertexTerm*propTerm;
    }
};

/**
 * Function to calculate the individual components of the self energy
 * @tparam Q            : Type of data, usually comp
 * @param ans           : SelfEnergy object to be updated.
 * @param fullvertex    : The Vertex object (ref) needed for the computation of the loop
 * @param prop          : The Propagator object (ref) needed for the computation of the loop
 * @param prop_iK       : Propagator index singled out for the calculation
 * @param self_iK       : Keldyh index the calculation is being carried out for
 * @param vert_iK       : Vertex index singled out for the calculation
 */
template <typename Q>
void loop_test_individual(SelfEnergy<Q>& ans, const Vertex<fullvert<Q> >& fullvertex, const Propagator& prop, int prop_iK, int self_iK, int vert_iK)
{
#pragma omp parallel for
    for (int iSE=0; iSE<nSE*n_in; ++iSE) {
        int i = iSE / n_in;
        int i_in = iSE - i * n_in;
        double v = ffreqs[i];

        //IntegrandSigma ocject created
        IntegrandSigma<Q, fullvert<Q>> integrandSigma(fullvertex, prop, vert_iK, v, i_in, prop_iK); //prop_iK = 0 for R, =-1 for A and =1 for K

        //Integration of the object. Notice the required limits for the integral
        Q integrated = -1./(2.*pi*glb_i)*integrator(integrandSigma, w_lower_f-fabs(v), w_upper_f+fabs(v));

        //The result is updated
        ans.addself(self_iK, i, integrated);
    }
}

/**
 * Function which calculates a SOPT state. Should however toggle off the components not to be computed.
 * @tparam Q    : Data type of the state, usually comp
 * @param Psi   : State whose Vertex is to be calculated
 * @param Lambda: Data structure-needed parameter. Should be set to 1. in all SOPT calculations
 * @param state : State whose Vertex whould be the bare vertex already initialized
 */
template<typename Q>
void sopt_state(State<Q>& Psi, double Lambda, const State<Q> &state) {

    //Calculate a propagator given the SelfEnergy of the initial condition-state
    Propagator g(Lambda, state.selfenergy, state.selfenergy, 'g');
    cout << "G calculated" << endl;

    //Bubbles are non-differentiated i.e. GG
    bool diff = false;

    cout << "bubble started" << endl;
    double t2 = get_time();
    //Lines 7-9
    //These lines calculate an a-Bubble. Toggle off if p-Bubble is on
    double ta = get_time();
    bubble_function(Psi.vertex, state.vertex, state.vertex, g, g, 'a', diff, '.');
    cout<<  "a - Bubble:";
    get_time(ta);

    //These lines calculate a p-Bubble. Toggle off if a-Bubble is on
//    double tp = get_time();
//    bubble_function(Psi.vertex, state.vertex, state.vertex, g, g, 'p', diff, '.');
//    cout<<  "p - Bubble:";
//    get_time(tp);

    //These lines calculate a t-Bubble. Toggle off for SOPT calculations
//    double tt = get_time();
//    bubble_function(Psi.vertex, state.vertex, state.vertex, g, g, 't', diff, '.');
//    cout<<  "t - Bubble:";
//    get_time(tt);

    cout << "bubble finished. ";
    get_time(t2);

}

void writeOutSOPT(double Lambda, const Propagator& propagator, const SelfEnergy<comp>& selfEnergy, const Vertex<fullvert<comp> >& vertex)
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
        my_file_sigmaR << ffreqs[j] << " " << selfEnergy.val(0, j).real() << " " <<  selfEnergy.val(0, j).imag() << "\n";
        my_file_sigmaA << ffreqs[j] << " " << selfEnergy.val(0, j).real() << " " << -selfEnergy.val(0, j).imag() << "\n";
        my_file_sigmaK << ffreqs[j] << " " << selfEnergy.val(1, j).real() << " " <<  selfEnergy.val(1, j).imag() << "\n";

        my_file_propR << ffreqs[j] << " " << propagator.valsmooth(0, ffreqs[j]).real() << " " <<  propagator.valsmooth(0, ffreqs[j]).imag() << "\n";
        my_file_propA << ffreqs[j] << " " << propagator.valsmooth(0, ffreqs[j]).real() << " " << -propagator.valsmooth(0, ffreqs[j]).imag() << "\n";
        my_file_propK << ffreqs[j] << " " << propagator.valsmooth(1, ffreqs[j]).real() << " " <<  propagator.valsmooth(1, ffreqs[j]).imag() << "\n";
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

/**
 * Function to test the OddOdd and OddEven bubbles
 * The implementation requires sopt_state() to have toggled the a-bubble on
 * @param state : State initialized with initial conditions
 */
void testBubbles(State<comp>& state, double Lambda){

    SelfEnergy<comp> zero;
    Propagator g1(Lambda, state.selfenergy, zero, 'g');
    Propagator g2(Lambda, state.selfenergy, zero, 'g');

    cvec PiaOO (nBOS);
    cvec PiaOE (nBOS);

    //Calculate the vertex for comparison
    sopt_state(state, 1.0, state);

    //for(int i=75; i<76; i++){
    for(int i=0; i<nBOS; i++){
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia6 (g1, g2, false, w, 6,  'a');     //AR
        IntegrandBubble integrandPia9 (g1, g2, false, w, 9,  'a');     //RA
        IntegrandBubble integrandPia11(g1, g2, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(g1, g2, false, w, 13, 'a');     //RK
        IntegrandBubble integrandPia15(g1, g2, false, w, 15, 'a');     //KK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        auto cont6  = integrator(integrandPia6 , w_lower_b, w_upper_b);
        auto cont9  = integrator(integrandPia9 , w_lower_b, w_upper_b);
        auto cont15 = integrator(integrandPia15, w_lower_b, w_upper_b);

        //Add the respective contributions to the respective bubble
        PiaOE[i] = 0.;//1./2.*(cont11+ cont13);
        PiaOO[i] = 1./2.*(cont6 + cont9);// + cont15);
        PiaOO[i]+= 1./2.*(1./(2.*pi*glb_i)*(correctionFunctionBubbleAT(w, -1., 1., w_upper_b, w_upper_b)+correctionFunctionBubbleAT(w, 1.,-1., w_upper_b, w_upper_b)));

/*
        //cout << Pia_odd_even[i] << endl;
        //cout << integrandPia11(-2.6) << endl;
        //cout << g1.pvalsmooth(1,-2.6-w/2.)*conj(g2.pvalsmooth(0,-2.6+w/2.))/(2.*pi*glb_i) << endl;
        cout << "GK(v-w/2): " << g1.pvalsmooth(1,-2.6-w/2.) << endl;


        cout << "GA(v+w/2): " << conj(g2.pvalsmooth(0,-2.6+w/2.)) << endl;
        cout << "GR(v-w/2): " << g2.pvalsmooth(0,-2.6-w/2.) << endl;
        cout << "GK_formula(v-w/2): " << (1.-2.*Fermi_distribution(-2.6-w/2.))*(g1.pvalsmooth(0,-2.6-w/2.)-conj(g1.pvalsmooth(0,-2.6-w/2.))) << endl;
        double myv = -0.1;
        cout << "myv: " << myv << endl;
        cout << "GK(myv): " << g1.pvalsmooth(1,myv) << endl;
        cout << "GK_formula(myv): " << (1.-2.*Fermi_distribution(myv))*(g1.pvalsmooth(0,myv)-conj(g1.pvalsmooth(0,myv))) << endl;

        myv = -0.2;
        cout << "myv: " << myv << endl;
        cout << "GK(myv): " << g1.pvalsmooth(1,myv) << endl;
        cout << "GK_formula(myv): " << (1.-2.*Fermi_distribution(myv))*(g1.pvalsmooth(0,myv)-conj(g1.pvalsmooth(0,myv))) << endl;

        myv = 0.;
        cout << "myv: " << myv << endl;
        cout << "GK(myv): " << g1.pvalsmooth(1,myv) << endl;
        cout << "GK_formula(myv): " << (1.-2.*Fermi_distribution(myv))*(g1.pvalsmooth(0,myv)-conj(g1.pvalsmooth(0,myv))) << endl;

        cout << g1.pvalsmooth(0, 0.) << endl;
        cout << g1.pvalsmooth(1, 0.) << endl;
*/

        //The relevant components are read out and added in the correct places directly.
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);


        //One can conversely use the K1_vvalsmooth (interpolating) functions to determine the value, which should not alter the result!
//        tvert<comp>* aid = &state.vertex.spinvertex.tvertex;

//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(1, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(1, w, 0, 1, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(7, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(7, w, 0, 1, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(8, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(8, w, 0, 1, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(14, w, 0, 0, *aid);
//        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(14, w, 0, 1, *aid);

//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);

//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(3, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(3, w, 0, 1, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(5, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(5, w, 0, 1, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(10, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(10, w, 0, 1, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(12, w, 0, 0, *aid);
//        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vvalsmooth(12, w, 0, 1, *aid);
    }

    write_h5_rvecs("SOPT_Bubbles.h5", {"v", "Bench_RePiaOO", "Bench_ImPiaOO", "Bench_RePiaOE", "Bench_ImPiaOE"},
                   {ffreqs, PiaOO.real(), PiaOO.imag(), PiaOE.real(), PiaOE.imag()});
}


/**
 * Function to test the OddOdd and OddEven bubbles
 * The implementation requires sopt_state() to have toggled the a-bubble on
 * @param state : State initialized with initial conditions
 */
void testBubblesFlow(){

    State<comp> test_state;
    //Initial conditions
    for (int i = 0; i < nSE; ++i) {
        test_state.selfenergy.setself(0, i, 0.);
        test_state.selfenergy.setself(1, i, 0.);
    }
    for (auto i:odd_Keldysh) {
        test_state.vertex.densvertex.irred.setvert(i, 0, 0.);
        test_state.vertex.spinvertex.irred.setvert(i, 0, -glb_U/2.);
    }

    Propagator g1(Lambda_ini, test_state.selfenergy, 'g');
    Propagator g2(Lambda_ini, test_state.selfenergy, 'g');

    for(int i=0; i<nBOS; i++) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia6(g1, g2, false, w, 6, 'a');     //AR
        IntegrandBubble integrandPia9(g1, g2, false, w, 9, 'a');     //RA
        IntegrandBubble integrandPia11(g1, g2, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(g1, g2, false, w, 13, 'a');     //RK
        IntegrandBubble integrandPia15(g1, g2, false, w, 15, 'a');     //KK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        auto cont6 = integrator(integrandPia6, w_lower_b, w_upper_b);
        auto cont9 = integrator(integrandPia9, w_lower_b, w_upper_b);
        auto cont15 = integrator(integrandPia15, w_lower_b, w_upper_b);

        //Add the respective contributions to the respective bubble
        test_state.vertex.spinvertex.avertex.K1_setvert(0, i, 0, 0.);//1./2.*(cont11 + cont13));              //11+13 => OE => Keldysh 0
        test_state.vertex.spinvertex.avertex.K1_setvert(1, i, 0, 1./2.*(cont6 + cont9));// + cont15));        //6+9+15=> OO => Keldysh 1
        test_state.vertex.spinvertex.avertex.K1_addvert(1, i, 0, 1./2.*(1./(2.*pi*glb_i)*(correctionFunctionBubbleAT(w, -1., 1., w_upper_b, w_upper_b)+correctionFunctionBubbleAT(w, 1.,-1., w_upper_b, w_upper_b))));
    }


    State<comp> final_state;
    export_data(test_state, 0);

    ODE_solver_RK4(final_state, 1.0, test_state, 2.0, test_rhs_bubbles_flow, nEVO);

}

/**
 * Function to test the loop function and the calculation of the SelfEnergy
 * @param state : State initialized with initial conditions
 */
void testSelfEnergy_and_Bubbles(State<comp>& state, double Lambda){

    Propagator g1(Lambda, state.selfenergy, 'g');

    //Calculate the vertex
    sopt_state(state, Lambda, state);

    //Calculate the SelfEnergy in SOPT
    loop(state.selfenergy, state.vertex, g1);

    //Print results in .g5 format
    cvec SER(nFER);
    cvec SEK(nFER);
    cvec PiaOO(nBOS);
    cvec PiaOE(nBOS);
    for(int i = 0; i<nFER; i++){
        SER[i] = state.selfenergy.val(0, i);
        SEK[i] = state.selfenergy.val(1, i);
    }
    for(int i = 0; i<nBOS; i++){
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);
        PiaOE[i] += state.vertex.spinvertex.avertex.K1_vval(0, i, 0);

        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
        PiaOO[i] += state.vertex.spinvertex.avertex.K1_vval(1, i, 0);
    }

    write_h5_rvecs("SOPT_SE.h5", {"v", "Bench_ReSER", "Bench_ImSER", "Bench_ReSEK", "Bench_ImSEK",
                                  "Bench_RePiaOO", "Bench_ImPiaOO", "Bench_RePiaOE", "Bench_ImPiaOE"},
                                 {ffreqs, SER.real(), SER.imag(), SEK.real(), SEK.imag(),
                                  PiaOO.real(), PiaOO.imag(), PiaOE.real(), PiaOE.imag()});

}

/**
 * Function to test the behaviour of all terms contributing to the self energy and the calculation of it
 * @param state : State initialized with initial conditions
 */
void test_selfEnergyComponents(State<comp>& state)
{
    SelfEnergy<comp> zero;
    Propagator g1(1., state.selfenergy, zero, 'g');

    SelfEnergy<comp> SER_R;
    SelfEnergy<comp> SER_A;
    SelfEnergy<comp> SER_K;
    SelfEnergy<comp> SEK_R;
    SelfEnergy<comp> SEK_A;
    SelfEnergy<comp> SEK_K;

    //Calculate the vertex
    sopt_state(state, 1.0, state);

    //Calculate the respective contributions to the self energy
    loop_test_individual(SER_R, state.vertex, g1, 0, 0, 3);
    loop_test_individual(SER_A, state.vertex, g1,-1, 0, 6);
    loop_test_individual(SER_K, state.vertex, g1, 1, 0, 7);
    loop_test_individual(SEK_R, state.vertex, g1, 0, 1, 1);
    loop_test_individual(SEK_A, state.vertex, g1,-1, 1, 4);
    loop_test_individual(SEK_K, state.vertex, g1, 1, 1, 5);

    //Print out results
    ostringstream selfEnergyR_R_cxx;
    ostringstream selfEnergyR_A_cxx;
    ostringstream selfEnergyR_K_cxx;
    ostringstream selfEnergyK_R_cxx;
    ostringstream selfEnergyK_A_cxx;
    ostringstream selfEnergyK_K_cxx;

    selfEnergyR_R_cxx << "Output/SelfEnergyR_R_cxx.dat";
    selfEnergyR_A_cxx << "Output/SelfEnergyR_A_cxx.dat";
    selfEnergyR_K_cxx << "Output/SelfEnergyR_K_cxx.dat";
    selfEnergyK_R_cxx << "Output/SelfEnergyK_R_cxx.dat";
    selfEnergyK_A_cxx << "Output/SelfEnergyK_A_cxx.dat";
    selfEnergyK_K_cxx << "Output/SelfEnergyK_K_cxx.dat";

    ofstream my_file_SelfEnergyR_R_cxx;
    ofstream my_file_SelfEnergyR_A_cxx;
    ofstream my_file_SelfEnergyR_K_cxx;
    ofstream my_file_SelfEnergyK_R_cxx;
    ofstream my_file_SelfEnergyK_A_cxx;
    ofstream my_file_SelfEnergyK_K_cxx;


    my_file_SelfEnergyR_R_cxx.open(selfEnergyR_R_cxx.str());
    my_file_SelfEnergyR_A_cxx.open(selfEnergyR_A_cxx.str());
    my_file_SelfEnergyR_K_cxx.open(selfEnergyR_K_cxx.str());
    my_file_SelfEnergyK_R_cxx.open(selfEnergyK_R_cxx.str());
    my_file_SelfEnergyK_A_cxx.open(selfEnergyK_A_cxx.str());
    my_file_SelfEnergyK_K_cxx.open(selfEnergyK_K_cxx.str());

    for(int i = 0; i<nFER; i++){
        my_file_SelfEnergyR_R_cxx << ffreqs[i]<< " " << SER_R.val(0, i).real() << " " << SER_R.val(0, i).imag() << "\n";
        my_file_SelfEnergyR_A_cxx << ffreqs[i]<< " " << SER_A.val(0, i).real() << " " << SER_A.val(0, i).imag() << "\n";
        my_file_SelfEnergyR_K_cxx << ffreqs[i]<< " " << SER_K.val(0, i).real() << " " << SER_K.val(0, i).imag() << "\n";
        my_file_SelfEnergyK_R_cxx << ffreqs[i]<< " " << SEK_R.val(1, i).real() << " " << SEK_R.val(1, i).imag() << "\n";
        my_file_SelfEnergyK_A_cxx << ffreqs[i]<< " " << SEK_A.val(1, i).real() << " " << SEK_A.val(1, i).imag() << "\n";
        my_file_SelfEnergyK_K_cxx << ffreqs[i]<< " " << SEK_K.val(1, i).real() << " " << SEK_K.val(1, i).imag() << "\n";

    }

    my_file_SelfEnergyR_R_cxx.close();
    my_file_SelfEnergyR_A_cxx.close();
    my_file_SelfEnergyR_K_cxx.close();
    my_file_SelfEnergyK_R_cxx.close();
    my_file_SelfEnergyK_A_cxx.close();
    my_file_SelfEnergyK_K_cxx.close();
}

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
