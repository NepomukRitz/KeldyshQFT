#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

#include <cmath>              // use M_PI as pi
#include "loop.h"             // self-energy loop
#include "solvers.h"          // ODE solvers
#include "write_data2file.h"  // writing data to txt or hdf5 file

/**
 * Function which calculates a SOPT state. Should however toggle off the components not to be computed.
 * @tparam Q    : Data type of the state, usually comp
 * @param Psi   : State whose Vertex is to be calculated
 * @param Lambda: Data structure-needed parameter. Should be set to 1. in all SOPT calculations
 * @param state : State whose Vertex whould be the bare vertex already initialized
 */
template<typename Q>
void sopt_state(State<Q>& Psi, double Lambda) {

    State<comp> bare;   //Create bare state
    bare.initialize();  //i.e. a state with a bare vertex and a self-energy initialized at the Hartree value

    Propagator g0(Lambda, bare.selfenergy, 'g');    //Bare propagator

    //Calculate the bubbles -> Vertex in SOPT
    bubble_function(Psi.vertex, bare.vertex, bare.vertex, g0, g0, 'a', false, '.');
    bubble_function(Psi.vertex, bare.vertex, bare.vertex, g0, g0, 'p', false, '.');
    bubble_function(Psi.vertex, bare.vertex, bare.vertex, g0, g0, 't', false, '.');

    //Do an a-Bubble for the calculation of the self-energy
    bubble_function(bare.vertex,  bare.vertex, bare.vertex, g0, g0, 'a', false, '.');

    //Calculate the Self-Energy
    loop(Psi.selfenergy, bare.vertex, g0, false);
}


/**
 * Function desgined to test the flow of the vertex in SOPT using the State class
 * @param input     : Input state (not used inside the function)
 * @param Lambda    : Lambda at which to calculate the rhs of the eq.
 * @return          : State carrying in the K1 vertex the results of the computation
 */
auto rhs_bubbles_flow_wstate(const State<comp>& input, double Lambda) -> State<comp>{
    State<comp> ans;    //Initialize the answer-object

    //Calculating propagator objects of the required types
    Propagator g(Lambda, input.selfenergy, 'g');
    Propagator s(Lambda, input.selfenergy, 's');

    bubble_function(ans.vertex, input.vertex, input.vertex, g, s, 'a', true, '.');
    return ans;
}

/**
 * Function to call when testing the rhs of the bubbles flow with the State class
 * @param N_ODE : Number of ODE steps to take between the globally defined Lambda_ini and Lambda_fin
 */
void test_rhs_bubbles_flow_wstate(int N_ODE) {
    bool write_flag = true; // whether to write output in hdf5
    State<comp> state_dir, state_fin, state_ini; // direct, final, initial K1a_1
    state_dir.initialize(); // initialize
    state_ini.initialize(); // initialize

    Propagator G0ini(Lambda_ini, state_ini.selfenergy, 'g'); // initial propagator
    Propagator G0dir(Lambda_fin, state_dir.selfenergy, 'g'); // final propagator

    sopt_state(state_ini, Lambda_ini); // direct calculation of initial K1a
    sopt_state(state_dir, Lambda_fin); // direct calculation of direct K1a

    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_bubbles_flow_wstate, N_ODE); // final K1a from ODE
    cvec K1a_dif = state_dir.vertex.spinvertex.avertex.K1 + ( state_fin.vertex.spinvertex.avertex.K1*(-1.) ); // difference in results
    cout << "Testing ODE for bare K1a_0 with State class. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("rhs_bubbles_flow_wstate.h5",
                                  {"v", "state_dir_R", "state_dir_I", "state_fin_R", "state_fin_I", "state_ini_R", "state_ini_I"},
                                  {bfreqs, state_dir.vertex.spinvertex.avertex.K1.real(), state_dir.vertex.spinvertex.avertex.K1.imag(),
                                                   state_fin.vertex.spinvertex.avertex.K1.real(), state_fin.vertex.spinvertex.avertex.K1.imag(),
                                                   state_ini.vertex.spinvertex.avertex.K1.real(), state_ini.vertex.spinvertex.avertex.K1.imag()});
}

/**
 * Function desgined to test the flow of the vertex in SOPT using the only cvecs
 * @param input     : Input cvec (not used inside the function)
 * @param Lambda    : Lambda at which to calculate the rhs of the eq.
 * @return          : The results of the calculation
 */
auto rhs_bubbles_flow(const cvec& input, double Lambda) -> cvec{
    cvec ans(nw1_wa);   //Initialize the answer

    SelfEnergy<comp> selfini;   //Initialize self energy
    selfini.initialize(glb_U/2., 0.);   //Hartree term

    Propagator g(Lambda, selfini, 'g'); //Regular propagator
    Propagator s(Lambda, selfini, 's'); //Single-scale propagator

    for(int i=0; i<nBOS; i++){
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia11(g, s, true, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(g, s, true, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        //Add the respective contributions to the respective bubble
        ans[i] = pow(-glb_U/2., 2.)*(cont11+ cont13);            //11+13 = OE => Keldysh comp0
    }

    return ans;
}

/**
 * Function to call when testing the rhs of the bubbles flow with cvecs
 * @param N_ODE : Number of ODE steps to take between the globally defined Lambda_ini and Lambda_fin
 */
void test_rhs_bubbles_flow(int N_ODE){
    bool write_flag = true; // whether to write output in hdf5
    cvec K1a_dir(nw1_wa), K1a_fin(nw1_wa), K1a_ini(nw1_wa); // direct, final, initial K1a_1
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    Propagator G0dir(Lambda_fin, SEin, 'g'); // final propagator

    // direct calculation of initial K1a
    for(int i=0; i<nw1_wa; ++i) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia11(G0ini, G0ini, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(G0ini, G0ini, false, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        K1a_ini[i] = pow(-glb_U/2.,2.)*(cont11 + cont13);
    }

    // direct calculation of direct K1a
    for(int i=0; i<nw1_wa; ++i) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia11(G0dir, G0dir, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(G0dir, G0dir, false, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, w_lower_b, w_upper_b);
        auto cont13 = integrator(integrandPia13, w_lower_b, w_upper_b);

        K1a_dir[i] = pow((-glb_U/2.),2.)*(cont11 + cont13);
    }

    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, rhs_bubbles_flow, N_ODE); // final K1a from ODE
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare K1a_0. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("rhs_bubbles_flow.h5",
                                  {"v", "K1a_dir_R", "K1a_dir_I", "K1a_fin_R", "K1a_fin_I", "K1a_ini_R", "K1a_ini_I"},
                                  {bfreqs, K1a_dir.real(), K1a_dir.imag(), K1a_fin.real(), K1a_fin.imag(), K1a_ini.real(), K1a_ini.imag()});
}



/**
 * Function to test the loop function and the calculation of the SelfEnergy
 * @param state : State initialized with initial conditions
 */
void testSelfEnergy_and_K1(State<comp>& state, double Lambda){

    Propagator g(Lambda, state.selfenergy, 'g');

    //Calculate the vertex
    sopt_state(state, Lambda);

    Vertex<fullvert<comp> > temp_vertex_a, temp_vertex_p; //All zeros
    temp_vertex_a.spinvertex.avertex = state.vertex.spinvertex.avertex;
    temp_vertex_p.spinvertex.pvertex = state.vertex.spinvertex.pvertex;

    loop(state.selfenergy, temp_vertex_a, g, false);//Calculate the SelfEnergy in SOPT

    //Print results in .h5 format
    cvec SER(nFER);
    cvec SEK(nFER);

    for(int iv = 0; iv<nFER; iv++){
        SER[iv] = state.selfenergy.val(0, iv, 0);
        SEK[iv] = state.selfenergy.val(1, iv, 0);
    }

    //Print results
    write_h5_rvecs("SOPT_test.h5", {"v", "ReSER", "ImSER", "ReSEK", "ImSEK",
                                  "K1a_R", "K1a_I", "K1p_R", "K1p_I", "K1t_R", "K1t_I"},
                                 {ffreqs, SER.real(), SER.imag(), SEK.real(), SEK.imag(),
                                  state.vertex.spinvertex.avertex.K1.real(),
                                  state.vertex.spinvertex.avertex.K1.imag(),
                                  state.vertex.spinvertex.pvertex.K1.real(),
                                  state.vertex.spinvertex.pvertex.K1.imag(),
                                  state.vertex.spinvertex.tvertex.K1.real(),
                                  state.vertex.spinvertex.tvertex.K1.imag()});

}


/**
 * Function to implement the flow of a State in SOPT.
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
auto rhs_state_flow_SOPT(const State<comp>& Psi, double Lambda) -> State<comp>{
    State<comp> dPsi;   //Answer object
    SelfEnergy<comp> SEhartree; //Hartree self-energy
    SEhartree.initialize(glb_U/2., 0.); //Initialization of the Hartree self-energy

    Propagator g(Lambda, SEhartree, 'g');    //Bare regular propagator at scale Lambda
    Propagator s(Lambda, SEhartree, 's');    //Bare single scale propagator at scale Lambda

    //Self energy flow
    loop(dPsi.selfenergy, Psi.vertex, s, true);  //Loop for the Self-Energy calculation

    //Vertex flow
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, g, s, 'a', true, '.');  //Differentiated bubble in the a-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, g, s, 'p', true, '.');  //Differentiated bubble in the p-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, g, s, 't', true, '.');  //Differentiated bubble in the t-channel

    return dPsi;
}

/**
 * Function to test the correctness of the flow of the State
 * @param N_ODE : Numbres of ODE-solver steps to be taken
 */
void test_rhs_state_flow_SOPT(int N_ODE){
    bool write_flag = true; // whether to write output in hdf5
    State<comp> state_dir, state_fin, state_ini; // direct, final, initial K1a_1
    state_dir.initialize(); // initialize state
    state_ini.initialize(); // initialize state

    Propagator G0ini(Lambda_ini, state_ini.selfenergy, 'g'); // initial propagator
    Propagator G0dir(Lambda_fin, state_dir.selfenergy, 'g'); // final propagator

    sopt_state(state_ini, Lambda_ini); // direct calculation of initial K1a
    sopt_state(state_dir, Lambda_fin); // direct calculation of direct K1a

    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT, N_ODE); // final K1a from ODE

    cvec K1a0_dif(nBOS);
    for(int i=0; i<nBOS; ++i){
        K1a0_dif[i] = state_dir.vertex.spinvertex.avertex.K1_val(0, i, 0) - state_fin.vertex.spinvertex.avertex.K1_val(0, i, 0);
    }
    cvec SER_dif(nFER);
    cvec SEK_dif(nFER);
    for(int i=0; i<nFER; ++i){
        SER_dif[i] = state_dir.selfenergy.val(0, i, 0) - state_fin.selfenergy.val(0, i, 0);
        SEK_dif[i] = state_dir.selfenergy.val(1, i, 0) - state_fin.selfenergy.val(1, i, 0);
    }

    cout << "Confirming correctness of the bubbles. The max diff between direct and final results is " << K1a0_dif.max_norm() << ". \n";
    cout << "Testing ODE for SelfEnergyR. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << SER_dif.max_norm() << "." << endl;
    cout << "Testing ODE for SelfEnergyK. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << SEK_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("rhs_state_flow_SOPT.h5",
                                  {"v", "dir_SE_R", "dir_SE_I", "fin_SE_R", "fin_SE_I", "ini_SE_R", "ini_SE_I",
                                   "dir_K1a_R", "dir_K1a_I", "fin_K1a_R", "fin_K1a_I", "ini_K1a_R", "ini_K1a_I",
                                   "dir_K1p_R", "dir_K1p_I", "fin_K1p_R", "fin_K1p_I", "ini_K1p_R", "ini_K1p_I",
                                   "dir_K1t_R", "dir_K1t_I", "fin_K1t_R", "fin_K1t_I", "ini_K1t_R", "ini_K1t_I"},
                                  {ffreqs,
                                   state_dir.selfenergy.Sigma.real(), state_dir.selfenergy.Sigma.imag(),
                                   state_fin.selfenergy.Sigma.real(), state_fin.selfenergy.Sigma.imag(),
                                   state_ini.selfenergy.Sigma.real(), state_ini.selfenergy.Sigma.imag(),
                                   state_dir.vertex.spinvertex.avertex.K1.real(), state_dir.vertex.spinvertex.avertex.K1.imag(),
                                   state_fin.vertex.spinvertex.avertex.K1.real(), state_fin.vertex.spinvertex.avertex.K1.imag(),
                                   state_ini.vertex.spinvertex.avertex.K1.real(), state_ini.vertex.spinvertex.avertex.K1.imag(),
                                   state_dir.vertex.spinvertex.pvertex.K1.real(), state_dir.vertex.spinvertex.pvertex.K1.imag(),
                                   state_fin.vertex.spinvertex.pvertex.K1.real(), state_fin.vertex.spinvertex.pvertex.K1.imag(),
                                   state_ini.vertex.spinvertex.pvertex.K1.real(), state_ini.vertex.spinvertex.pvertex.K1.imag(),
                                   state_dir.vertex.spinvertex.tvertex.K1.real(), state_dir.vertex.spinvertex.tvertex.K1.imag(),
                                   state_fin.vertex.spinvertex.tvertex.K1.real(), state_fin.vertex.spinvertex.tvertex.K1.imag(),
                                   state_ini.vertex.spinvertex.tvertex.K1.real(), state_ini.vertex.spinvertex.tvertex.K1.imag()});

}



/**
 * Function that prints out a .h5 file with the value of the rhs of a SOPT flow at the given Lambda for both an FFT and an fRG calculation
 * @param Lambda    : Lambda at which the derivatives are to be calculated
 */
void test_derivatives_K1a(double Lambda){
    cvec blah(nw1_wa);
    cvec rhs_SOPT_FFT_K1a = dSOPT_FFT_K1a_rhs(blah, Lambda);
    cvec rhs_flow = rhs_bubbles_flow(blah, Lambda);

    write_h5_rvecs("derivatives_K1a.h5",
                   {"v", "FFT_R", "FFT_I", "SOPT_R", "SOPT_I"},
                   {ffreqs, rhs_SOPT_FFT_K1a.real(), rhs_SOPT_FFT_K1a.imag(), rhs_flow.real(), rhs_flow.imag()});
}

/**
 * Function that prints out a .h5 file with the value of the rhs of a SOPT flow at the given Lambda for both an FFT and an fRG calculation
 * @param Lambda    : Lambda at which the derivatives are to be calculated
 */
void test_derivatives_SE(double Lambda){
    cvec blah(nw1_wa);
    State<comp> sopt;
    sopt.initialize();
    cvec rhs_SOPT_FFT_K1a = dSOPT_FFT_SE_rhs(blah, Lambda);

    sopt_state(sopt, Lambda);
    State<comp> rhs_flow = rhs_state_flow_SOPT(sopt, Lambda);

    cvec SER_dif(nFER);
    for(int iv=0; iv<nFER;++iv){
        SER_dif[iv] = rhs_flow.selfenergy.val(0, iv, 0) -rhs_SOPT_FFT_K1a[iv];
    }

    cout << "Testing derivatives of the Self Energy. Using Lambda=" << Lambda << ", the max difference between fRG and FFT results is " << SER_dif.max_norm() << "." << endl;
    write_h5_rvecs("derivativess_SE.h5",
                   {"v", "FFT_R", "FFT_I", "SOPT_R", "SOPT_I"},
                   {ffreqs, rhs_SOPT_FFT_K1a.real(), rhs_SOPT_FFT_K1a.imag(),
                    rhs_flow.selfenergy.Sigma.real(), rhs_flow.selfenergy.Sigma.imag()});
}

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
