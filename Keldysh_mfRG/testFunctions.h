#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

#include <cmath>              // use M_PI as pi
#include "loop.h"             // self-energy loop
#include "solvers.h"          // ODE solvers
#include "write_data2file.h"  // writing data to txt or hdf5 file
#include "hdf5_routines.h"    // writing states to hdf5 file
#include "perturbation_theory.h"




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
    cvec K1a_dif = state_dir.vertex[0].avertex.K1 + ( state_fin.vertex[0].avertex.K1*(-1.) ); // difference in results
    print("Testing ODE for bare K1a_0 with State class. Using " +to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +to_string(K1a_dif.max_norm())+ ".", true);
    if(write_flag) write_h5_rvecs("rhs_bubbles_flow_wstate.h5",
                                  {"v", "state_dir_R", "state_dir_I", "state_fin_R", "state_fin_I", "state_ini_R", "state_ini_I"},
                                  {bfreqs, state_dir.vertex[0].avertex.K1.real(), state_dir.vertex[0].avertex.K1.imag(),
                                                   state_fin.vertex[0].avertex.K1.real(), state_fin.vertex[0].avertex.K1.imag(),
                                                   state_ini.vertex[0].avertex.K1.real(), state_ini.vertex[0].avertex.K1.imag()});
}

/**
 * Function desgined to test the flow of the vertex in SOPT using the only cvecs
 * @param input     : Input cvec (not used inside the function)
 * @param Lambda    : Lambda at which to calculate the rhs of the eq.
 * @return          : The results of the calculation
 */
auto rhs_bubbles_flow(const cvec& input, double Lambda) -> cvec{
    cvec ans(nw1_a);   //Initialize the answer

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
        auto cont11 = integrator(integrandPia11, glb_w_lower, glb_w_lower);
        auto cont13 = integrator(integrandPia13, glb_w_lower, glb_w_lower);

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
    cvec K1a_dir(nw1_a), K1a_fin(nw1_a), K1a_ini(nw1_a); // direct, final, initial K1a_1
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    Propagator G0dir(Lambda_fin, SEin, 'g'); // final propagator

    // direct calculation of initial K1a
    for(int i=0; i<nw1_a; ++i) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia11(G0ini, G0ini, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(G0ini, G0ini, false, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, glb_w_lower, glb_w_lower);
        auto cont13 = integrator(integrandPia13, glb_w_lower, glb_w_lower);

        K1a_ini[i] = pow(-glb_U/2.,2.)*(cont11 + cont13);
    }

    // direct calculation of direct K1a
    for(int i=0; i<nw1_a; ++i) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble integrandPia11(G0dir, G0dir, false, w, 11, 'a');     //KA
        IntegrandBubble integrandPia13(G0dir, G0dir, false, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator(integrandPia11, glb_w_lower, glb_w_lower);
        auto cont13 = integrator(integrandPia13, glb_w_lower, glb_w_lower);

        K1a_dir[i] = pow((-glb_U/2.),2.)*(cont11 + cont13);
    }

    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, rhs_bubbles_flow, N_ODE); // final K1a from ODE
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) ); // difference in results
    print("Testing ODE for bare K1a_0. Using " +to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +to_string(K1a_dif.max_norm())+ ".", true);
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

    Vertex<comp> temp_vertex_a (1), temp_vertex_p (1); //All zeros
    temp_vertex_a[0].avertex = state.vertex[0].avertex;
    temp_vertex_p[0].pvertex = state.vertex[0].pvertex;

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
                                  state.vertex[0].avertex.K1.real(),
                                  state.vertex[0].avertex.K1.imag(),
                                  state.vertex[0].pvertex.K1.real(),
                                  state.vertex[0].pvertex.K1.imag(),
                                  state.vertex[0].tvertex.K1.real(),
                                  state.vertex[0].tvertex.K1.imag()});

}


/**
 * Function to implement the flow of a State in SOPT.
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
auto rhs_state_flow_SOPT(const State<comp>& Psi, const double Lambda, const int feedback) -> State<comp>{
    State<comp> dPsi;   //Answer object

    State<comp> bare;   //Bare state
    bare.initialize();  //Initialize bare state
    Propagator G(Lambda, bare.selfenergy,'g');    //Initialization of Propagator objects
    Propagator S(Lambda, bare.selfenergy,'s');    //Initialization of Propagator objects

    if(!(feedback==0 || feedback==3)){  //Check whether Self Energy feedback to Propagators is wanted
        G=Propagator(Lambda, Psi.selfenergy, 'g');
        S=Propagator(Lambda, Psi.selfenergy, 's');
    }

    //Self energy loop
    loop(dPsi.selfenergy, Psi.vertex, S, true);  //Loop for the Self-Energy calculation


    if(feedback>=3){    //If feedback>=3, there is vertex feedback
        bare.vertex = Psi.vertex;
    }

    if(feedback==2 || feedback==5) {  //These two options make use of the Katanin substitution
        S=Propagator(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
    }

    //Vertex flow
    bubble_function(dPsi.vertex, bare.vertex, bare.vertex, G, S, 'a', true, '.');  //Differentiated bubble in the a-channel
    bubble_function(dPsi.vertex, bare.vertex, bare.vertex, G, S, 'p', true, '.');  //Differentiated bubble in the p-channel
    bubble_function(dPsi.vertex, bare.vertex, bare.vertex, G, S, 't', true, '.');  //Differentiated bubble in the t-channel

    return dPsi;
}

//No feedback
auto rhs_state_flow_SOPT_0(const State<comp>& Psi, const double Lambda) -> State<comp>{
    return rhs_state_flow_SOPT(Psi, Lambda, 0);
}

//Self-Energy fed back into Propagators
auto rhs_state_flow_SOPT_1(const State<comp>& Psi, const double Lambda) -> State<comp>{
    return rhs_state_flow_SOPT(Psi, Lambda, 1);
}

//As above + Katanin
auto rhs_state_flow_SOPT_2(const State<comp>& Psi, const double Lambda) -> State<comp>{
    return rhs_state_flow_SOPT(Psi, Lambda, 2);
}

//Vertex feedback, free propagators
auto rhs_state_flow_SOPT_3(const State<comp>& Psi, const double Lambda) -> State<comp>{
    return rhs_state_flow_SOPT(Psi, Lambda, 3);
}

//Self_energy fed back into Propagators + Vertex feedback
auto rhs_state_flow_SOPT_4(const State<comp>& Psi, const double Lambda) -> State<comp>{
    return rhs_state_flow_SOPT(Psi, Lambda, 4);
}

//As above + Katanin
auto rhs_state_flow_SOPT_5(const State<comp>& Psi, const double Lambda) -> State<comp>{
    return rhs_state_flow_SOPT(Psi, Lambda, 5);
}

/**
 * Function to test the correctness of the flow of the State
 * @param N_ODE : Numbres of ODE-solver steps to be taken
 */
void test_rhs_state_flow_SOPT(int N_ODE, int feedback){
    bool write_flag = true; // whether to write output in hdf5
    State<comp> state_dir, state_fin, state_ini; // direct, final, initial K1a_1
    state_dir.initialize(); // initialize state
    state_ini.initialize(); // initialize state

    Propagator G0ini(Lambda_ini, state_ini.selfenergy, 'g'); // initial propagator
    Propagator G0dir(Lambda_fin, state_dir.selfenergy, 'g'); // final propagator

    sopt_state(state_ini, Lambda_ini); // direct calculation of initial K1a
    sopt_state(state_dir, Lambda_fin); // direct calculation of direct K1a

    string name;
    switch (feedback){
        case 1:
            ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT_1, N_ODE); // final K1a from ODE
            name = "rhs_state_flow_feedback_1.h5";
            break;
        case 2:
            ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT_2, N_ODE); // final K1a from ODE
            name = "rhs_state_flow_feedback_2.h5";
            break;
        case 3:
            ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT_3, N_ODE); // final K1a from ODE
            name = "rhs_state_flow_feedback_3.h5";
            break;
        case 4:
            ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT_4, N_ODE); // final K1a from ODE
            name = "rhs_state_flow_feedback_4.h5";
            break;
        case 5:
            ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT_5, N_ODE); // final K1a from ODE
            name = "rhs_state_flow_feedback_5.h5";
            break;
        default:
            ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_state_flow_SOPT_0, N_ODE); // final K1a from ODE
            name = "rhs_state_flow_feedback_0.h5";
    }

    cvec K1a0_dif(nBOS);
    for(int i=0; i<nBOS; ++i){
        K1a0_dif[i] = state_dir.vertex[0].avertex.K1_val(0, i, 0) - state_fin.vertex[0].avertex.K1_val(0, i, 0);
    }
    cvec SER_dif(nFER);
    cvec SEK_dif(nFER);
    for(int i=0; i<nFER; ++i){
        SER_dif[i] = state_dir.selfenergy.val(0, i, 0) - state_fin.selfenergy.val(0, i, 0);
        SEK_dif[i] = state_dir.selfenergy.val(1, i, 0) - state_fin.selfenergy.val(1, i, 0);
    }

    if(mpi_world_rank()==0) {
        print("Confirming correctness of the bubbles. The max diff between direct and final results is " +to_string(K1a0_dif.max_norm())+". \n "+
        "Testing ODE for SelfEnergyR. Using " +to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +
        to_string(SER_dif.max_norm())+ ". \n" +
        "Testing ODE for SelfEnergyK. Using " + to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +
        to_string(SEK_dif.max_norm()) +".", true);
    }
    if(write_flag) write_h5_rvecs(name,
                                  {"v", "dir_SE_R", "dir_SE_I", "fin_SE_R", "fin_SE_I", "ini_SE_R", "ini_SE_I",
                                   "dir_K1a_R", "dir_K1a_I", "fin_K1a_R", "fin_K1a_I", "ini_K1a_R", "ini_K1a_I",
                                   "dir_K1p_R", "dir_K1p_I", "fin_K1p_R", "fin_K1p_I", "ini_K1p_R", "ini_K1p_I",
                                   "dir_K1t_R", "dir_K1t_I", "fin_K1t_R", "fin_K1t_I", "ini_K1t_R", "ini_K1t_I"},
                                  {ffreqs,
                                   state_dir.selfenergy.Sigma.real(), state_dir.selfenergy.Sigma.imag(),
                                   state_fin.selfenergy.Sigma.real(), state_fin.selfenergy.Sigma.imag(),
                                   state_ini.selfenergy.Sigma.real(), state_ini.selfenergy.Sigma.imag(),
                                   state_dir.vertex[0].avertex.K1.real(), state_dir.vertex[0].avertex.K1.imag(),
                                   state_fin.vertex[0].avertex.K1.real(), state_fin.vertex[0].avertex.K1.imag(),
                                   state_ini.vertex[0].avertex.K1.real(), state_ini.vertex[0].avertex.K1.imag(),
                                   state_dir.vertex[0].pvertex.K1.real(), state_dir.vertex[0].pvertex.K1.imag(),
                                   state_fin.vertex[0].pvertex.K1.real(), state_fin.vertex[0].pvertex.K1.imag(),
                                   state_ini.vertex[0].pvertex.K1.real(), state_ini.vertex[0].pvertex.K1.imag(),
                                   state_dir.vertex[0].tvertex.K1.real(), state_dir.vertex[0].tvertex.K1.imag(),
                                   state_fin.vertex[0].tvertex.K1.real(), state_fin.vertex[0].tvertex.K1.imag(),
                                   state_ini.vertex[0].tvertex.K1.real(), state_ini.vertex[0].tvertex.K1.imag()});

}



/**
 * Function that prints out a .h5 file with the value of the rhs of a SOPT flow at the given Lambda for both an FFT and an fRG calculation
 * @param Lambda    : Lambda at which the derivatives are to be calculated
 */
void test_derivatives_K1a(double Lambda){
    cvec blah(nw1_a);
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
    cvec blah(nw1_a);
    State<comp> sopt;
    sopt.initialize();
    cvec rhs_SOPT_FFT_K1a = dSOPT_FFT_SE_rhs(blah, Lambda);

    sopt_state(sopt, Lambda);
    State<comp> rhs_flow = rhs_state_flow_SOPT_0(sopt, Lambda);

    cvec SER_dif(nFER);
    for(int iv=0; iv<nFER;++iv){
        SER_dif[iv] = rhs_flow.selfenergy.val(0, iv, 0) -rhs_SOPT_FFT_K1a[iv];
    }

    print("Testing derivatives of the Self Energy. Using Lambda=" +to_string(Lambda)+ ", the max difference between fRG and FFT results is " +to_string(SER_dif.max_norm())+ ".", true);
    write_h5_rvecs("derivativess_SE.h5",
                   {"v", "FFT_R", "FFT_I", "SOPT_R", "SOPT_I"},
                   {ffreqs, rhs_SOPT_FFT_K1a.real(), rhs_SOPT_FFT_K1a.imag(),
                    rhs_flow.selfenergy.Sigma.real(), rhs_flow.selfenergy.Sigma.imag()});
}

#if DIAG_CLASS >= 2
/**
 * This function checks the consistency of the K2 class
 * @param Lambda    : Scale
 * @param r         : Channel to be tested i.e. K_2^r
 */
auto test_K2_consistency(double Lambda, const char r) -> bool{
    State<comp> bare;   //Create a bare state
    bare.initialize();  //Initialize bare state

    Propagator G(Lambda, bare.selfenergy, 'g'); //Bare propagator at scale Lambda>

    //Create and initialize states to save the channel contributions in
    State<comp> K1a;
    State<comp> K1p;
    State<comp> K1t;
    K1a.initialize();
    K1p.initialize();
    K1t.initialize();

    //Calculate K1a, K1p and K1t contributions separately
    bubble_function(K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false, '.');
    bubble_function(K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false, '.');
    bubble_function(K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false, '.');

    //Create and initialize the K2r-objects to test
    State<comp> test_K2r_with_K1a;
    State<comp> test_K2r_with_K1p;
    State<comp> test_K2r_with_K1t;
    test_K2r_with_K1a.initialize();
    test_K2r_with_K1p.initialize();
    test_K2r_with_K1t.initialize();


    if(r=='p' ||  r=='t') {
        //Perform TOPT calculation of K2r with K1a
        bubble_function(test_K2r_with_K1a.vertex, K1a.vertex, bare.vertex, G, G, r, false, 'L');
    }

    if(r=='a' || r=='t') {
        //Perform TOPT calculation of K2r with K1p
        bubble_function(test_K2r_with_K1p.vertex, K1p.vertex, bare.vertex, G, G, r, false, 'L');
    }

    if(r=='a' || r=='p') {
        //Perform TOPT calculation of K2r with K1t
        bubble_function(test_K2r_with_K1t.vertex, K1t.vertex, bare.vertex, G, G, r, false, 'L');
    }


    //TOPT-consistency checks for K2

    bool empty = true;  //Boolean that stays true if everything is zero
    if(r=='a' || r=='p') {
#pragma omp parallel
        //Check parallelized that everything in the K2 vertex is zero.
        for (int index = 0; index < test_K2r_with_K1t.vertex[0].avertex.K2.size(); ++index) {
            if (test_K2r_with_K1t.vertex[0].avertex.K2_acc(index) != 0.) {
                empty = false;  //If any one element is not zero, K2 is not empty
                break;  //Exit the for-loop early
            }
        }
        //Print result of consistency check
        if (empty) print("TOPT-consistency check passed. K2a or K2p with K1t is zero everywhere.", true);
        else print("TOPT-consistency check failed. K2a or K2p with K1t is not zero everywhere.", true);
    }
    else if(r=='t'){
#pragma omp parallel
        //Check parallelized that everything in the K2 vertices is zero.
        for (int index = 0; index < test_K2r_with_K1a.vertex[0].avertex.K2.size(); ++index) {
            if (test_K2r_with_K1a.vertex[0].avertex.K2_acc(index) != 0. ||
                test_K2r_with_K1p.vertex[0].avertex.K2_acc(index) != 0.) {
                empty = false;  //If any one element is not zero, K2 is not empty
                break;  //Exit the for-loop early
            }
        }
        //Print result of consistency check
        if (empty) print("TOPT-consistency check passed. K2t with K1a and K1p are zero everywhere.", true);
        else print("TOPT-consistency check failed. K2t with K1a and K1p are not zero everywhere.", true);
    }

    return empty;
}
#endif
/**
 * Function to test correctness of K2a when calculating a susceptibility (a K1a-object) //Notice that the same calculation
 * can be performed in the p-channel.
 * @param Lambda : Scale at which the calculation is done.
 */
void test_K2_correctness(double Lambda){

    bool write_flag = true; //Write out results in a HDF5 file

    State<comp> bare;   //Bare state
    bare.initialize();  //Initialize bare state

    Propagator G(Lambda, bare.selfenergy, 'g'); //Bare propagator
    Propagator S(Lambda, bare.selfenergy, 's'); //Bare single-scale propagator

    //Create states for K1-calculations
    State<comp> PT2_K1a;
    State<comp> PT2_K1p;
    State<comp> PT2_K1t;

    //Save K1-bubbles in separate objects - SOPT
    double t0 = get_time();
    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false, '.');
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false, '.');
    bubble_function(PT2_K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false, '.');
    get_time(t0);

    State<comp> PT2_SE_a;
    State<comp> PT2_SE_p;
    State<comp> PT2_SE_t;
    State<comp> PT2_SE_p_1;
    State<comp> PT2_SE_p_4;
    State<comp> PT2_SE_p_5;

    loop(PT2_SE_a.selfenergy, PT2_K1a.vertex, S, true);
    loop(PT2_SE_p.selfenergy, PT2_K1p.vertex, S, true);
    loop(PT2_SE_t.selfenergy, PT2_K1t.vertex, S, true);
    loop(PT2_SE_p_1.selfenergy, PT2_K1p.vertex, S, true, 1);
    loop(PT2_SE_p_4.selfenergy, PT2_K1p.vertex, S, true, 4);
    loop(PT2_SE_p_5.selfenergy, PT2_K1p.vertex, S, true, 5);

    State<comp> PT3_K2a;    //Create state for K2a calculation
    State<comp> PT3_K2a_ia;
    State<comp> PT3_K2a_ib;
    State<comp> PT3_K2a_iia;
    State<comp> PT3_K2a_iib;
    State<comp> PT3_K2a_iva;
    State<comp> PT3_K2a_ivb;
    State<comp> PT3_K2a_t;

    State<comp> PT3_K2p;
    State<comp> PT3_K2t;
    State<comp> PT3_K2t_a;
    State<comp> PT3_K2t_p;

    //Do appropriate calculation for K2a with K1p and K1t being fed back into the left vertex. Notice part = 'L' to ensure
    //that the correct contributions are added on both sides. - TOPT
    t0 = get_time();
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L');   // K2a in PT3

    bubble_function(PT3_K2a_ia.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 9, 6);
    bubble_function(PT3_K2a_ib.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 6, 9);
    bubble_function(PT3_K2a_iia.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 11, 14);
    bubble_function(PT3_K2a_iib.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 7, 9);
    bubble_function(PT3_K2a_iva.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 9, 6);
    bubble_function(PT3_K2a_ivb.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 6, 15);

    bubble_function(PT3_K2a_t.vertex, PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L');   // K2a in PT3
    rvec Lambdas {1.};
    //PT3_K2a = read_hdf("PT4_check_of_K2a_K2_switchedcc_adap_m3m9_g501_101_nI1501_state_PT3_K2a", 0, 1, Lambdas);

    bubble_function(PT3_K2p.vertex, PT2_K1a.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'p', false, 'L');    // K2p  in PT3
    //bubble_function(PT3_K2t.vertex, PT2_K1a.vertex + PT2_K1p.vertex, bare.vertex, G, G, 't', false, 'L');    // K2t  in PT3
    bubble_function(PT3_K2t_a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 't', false, 'L');    // K2t  in PT3
    bubble_function(PT3_K2t_p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 't', false, 'L');    // K2t  in PT3
    PT3_K2t.vertex = PT3_K2t_a.vertex + PT3_K2t_p.vertex;
    get_time(t0);

    // full K2 in PT3
    State<comp> PT3_K2;
    PT3_K2.vertex[0].avertex = PT3_K2a_t.vertex[0].avertex;
    //PT3_K2.vertex[0].pvertex = PT3_K2p.vertex[0].pvertex;
    PT3_K2.vertex[0].tvertex = PT3_K2t_a.vertex[0].tvertex;

    State<comp> PT3_K2at;

    // K2 contribution to self-energy flow
    State<comp> PT3_SE;
    State<comp> PT3_SE_a;
    State<comp> PT3_SE_p;
    State<comp> PT3_SE_t;
    State<comp> PT3_SE_t_1;
    State<comp> PT3_SE_t_4;
    State<comp> PT3_SE_t_5;
    State<comp> PT3_SE_t_a;
    State<comp> PT3_SE_t_a_1;
    State<comp> PT3_SE_t_a_4;
    State<comp> PT3_SE_t_a_5;
    State<comp> PT3_SE_t_p;
    State<comp> PT3_SE_t_p_1;
    State<comp> PT3_SE_t_p_4;
    State<comp> PT3_SE_t_p_5;

    loop(PT3_SE.selfenergy, PT3_K2.vertex, S, true);
    loop(PT3_SE_a.selfenergy, PT3_K2a.vertex, S, true);
    loop(PT3_SE_p.selfenergy, PT3_K2p.vertex, S, true);
    loop(PT3_SE_t.selfenergy, PT3_K2t.vertex, S, true);
    loop(PT3_SE_t_1.selfenergy, PT3_K2t.vertex, S, true, 1);
    loop(PT3_SE_t_4.selfenergy, PT3_K2t.vertex, S, true, 4);
    loop(PT3_SE_t_5.selfenergy, PT3_K2t.vertex, S, true, 5);
    loop(PT3_SE_t_a.selfenergy, PT3_K2t_a.vertex, S, true);
    loop(PT3_SE_t_a_1.selfenergy, PT3_K2t_a.vertex, S, true, 1);
    loop(PT3_SE_t_a_4.selfenergy, PT3_K2t_a.vertex, S, true, 4);
    loop(PT3_SE_t_a_5.selfenergy, PT3_K2t_a.vertex, S, true, 5);
    loop(PT3_SE_t_p.selfenergy, PT3_K2t_p.vertex, S, true);
    loop(PT3_SE_t_p_1.selfenergy, PT3_K2t_p.vertex, S, true, 1);
    loop(PT3_SE_t_p_4.selfenergy, PT3_K2t_p.vertex, S, true, 4);
    loop(PT3_SE_t_p_5.selfenergy, PT3_K2t_p.vertex, S, true, 5);

    State<comp> PT3_K1a;    //Create state to compare with K1a
    t0 = get_time();
    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false, '.');
    get_time(t0);

    State<comp> PT123_a = bare + PT2_K1a + PT3_K1a + PT3_K2a;  //Create vertex of the right side of BSE

    State<comp> PT4_K1a22;
    State<comp> PT4_K1a13_1;
    State<comp> PT4_K1a13_2;
    State<comp> PT4_K1a31_2;
    State<comp> PT4_K1a13_2_11e; // A
    State<comp> PT4_K1a13_2_21e; // B
    State<comp> PT4_K1a13_2_11o; // C
    State<comp> PT4_K1a13_2_21o; // D
    State<comp> PT4_K1a13_2_12o; // TST3TC D
    State<comp> PT4_K1a13_2_22o; // F
    State<comp> PT4_K1a13_2_ia;
    State<comp> PT4_K1a13_2_ib;
    State<comp> PT4_K1a13_2_iia;
    State<comp> PT4_K1a13_2_iib;
    State<comp> PT4_K1a13_2_iva;
    State<comp> PT4_K1a13_2_ivb;

    //Calculate a K1a-object to compare with K1a and NRG-results. Notice part='R', suggesting the "weird" vertex is on the }
    //right and, thanks to this, no unnecessary K2 calculation is entered in bubble_function
    t0 = get_time();
    bubble_function(PT4_K1a22.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false, 'R');
    bubble_function(PT4_K1a13_1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false, 'R');
    glb_int_flag = true;
    bubble_function(PT4_K1a13_2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R');
    bubble_function(PT4_K1a31_2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false, 'R');
    get_time(t0);
    t0 = get_time();
    bubble_function(PT4_K1a13_2_11e.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R', 0, 16, 16, 16); // A
    bubble_function(PT4_K1a13_2_21e.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R', 1, 16, 16, 16); // B
    bubble_function(PT4_K1a13_2_11o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R', 2, 16, 16, 16); // C
    bubble_function(PT4_K1a13_2_21o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R', 3, 16, 16, 16); // D
    bubble_function(PT4_K1a13_2_12o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R', 10, 16, 16, 16); // TST3TC D
    bubble_function(PT4_K1a13_2_22o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 'R', 11, 16, 16, 16); // F

    bubble_function(PT4_K1a13_2_ia.vertex, bare.vertex, PT3_K2a_ia.vertex, G, G, 'a', false, 'R', 11, 6, 16, 16);
    bubble_function(PT4_K1a13_2_ib.vertex, bare.vertex, PT3_K2a_ib.vertex, G, G, 'a', false, 'R', 2, 9, 16, 16);
    bubble_function(PT4_K1a13_2_iia.vertex, bare.vertex, PT3_K2a_iia.vertex, G, G, 'a', false, 'R', 11, 6, 16, 16);
    bubble_function(PT4_K1a13_2_iib.vertex, bare.vertex, PT3_K2a_iib.vertex, G, G, 'a', false, 'R', 10, 11, 16, 16);
    bubble_function(PT4_K1a13_2_iva.vertex, bare.vertex, PT3_K2a_iva.vertex, G, G, 'a', false, 'R', 11, 15, 16, 16);
    bubble_function(PT4_K1a13_2_ivb.vertex, bare.vertex, PT3_K2a_ivb.vertex, G, G, 'a', false, 'R', 2, 9, 16, 16);
    get_time(t0);

    cvec K1a_diff(nBOS);
    for(int iw=0; iw<nBOS; ++iw){
        K1a_diff[iw] = PT4_K1a22.vertex[0].avertex.K1_val(0, iw, 0) - PT2_K1a.vertex[0].avertex.K1_val(0, iw, 0);
    }

    print("Testing correctness of K2a. Using U=" +to_string(glb_U)+ " and Lambda="+to_string(Lambda)+", the maximal difference between direct K1a and K1a over integration of K2a is " +to_string(K1a_diff.max_norm())+"." , true);
    if(write_flag) write_h5_rvecs("PT4_check_of_K2a_K2_switchedcc_t_update_symmrev_new8_SE_symm_full_adap_m3m9_gW10_501_101_nI1501_U1", {"w",
                                                       "PT2_K1a_R", "PT2_K1a_I",
                                                       "PT2_K1p_R", "PT2_K1p_I",
                                                       "PT2_K1t_R", "PT2_K1t_I",
                                                       "PT2_SE_a_R", "PT2_SE_a_I",
                                                       "PT2_SE_p_R", "PT2_SE_p_I",
                                                       "PT2_SE_t_R", "PT2_SE_t_I",
                                                       "PT2_SE_p_1_R", "PT2_SE_p_1_I",
                                                       "PT2_SE_p_4_R", "PT2_SE_p_4_I",
                                                       "PT2_SE_p_5_R", "PT2_SE_p_5_I",
                                                       "PT3_K1a_R", "PT3_K1a_I",
                                                       "PT3_K2a_R", "PT3_K2a_I",
                                                       "PT3_K2a_ia_R", "PT3_K2a_ia_I",
                                                       "PT3_K2a_ib_R", "PT3_K2a_ib_I",
                                                       "PT3_K2a_iia_R", "PT3_K2a_iia_I",
                                                       "PT3_K2a_iib_R", "PT3_K2a_iib_I",
                                                       "PT3_K2a_iva_R", "PT3_K2a_iva_I",
                                                       "PT3_K2a_ivb_R", "PT3_K2a_ivb_I",
                                                       "PT3_K2a_t_R", "PT3_K2a_t_I",
                                                       "PT3_K2p_R", "PT3_K2p_I",
                                                       "PT3_K2t_R", "PT3_K2t_I",
                                                       "PT3_K2t_a_a_R", "PT3_K2t_a_a_I",
                                                       "PT3_K2t_a_p_R", "PT3_K2t_a_p_I",
                                                       "PT3_K2t_a_t_R", "PT3_K2t_a_t_I",
                                                       "PT3_SE_R", "PT3_SE_I",
                                                       "PT3_SE_a_R", "PT3_SE_a_I",
                                                       "PT3_SE_p_R", "PT3_SE_p_I",
                                                       "PT3_SE_t_R", "PT3_SE_t_I",
                                                       "PT3_SE_t_1_R", "PT3_SE_t_1_I",
                                                       "PT3_SE_t_4_R", "PT3_SE_t_4_I",
                                                       "PT3_SE_t_5_R", "PT3_SE_t_5_I",
                                                       "PT3_SE_t_a_R", "PT3_SE_t_a_I",
                                                       "PT3_SE_t_a_1_R", "PT3_SE_t_a_1_I",
                                                       "PT3_SE_t_a_4_R", "PT3_SE_t_a_4_I",
                                                       "PT3_SE_t_a_5_R", "PT3_SE_t_a_5_I",
                                                       "PT3_SE_t_p_R", "PT3_SE_t_p_I",
                                                       "PT3_SE_t_p_1_R", "PT3_SE_t_p_1_I",
                                                       "PT3_SE_t_p_4_R", "PT3_SE_t_p_4_I",
                                                       "PT3_SE_t_p_5_R", "PT3_SE_t_p_5_I",
                                                       "PT4_K1a22_R", "PT4_K1a22_I",
                                                       "PT4_K1a13_1_R", "PT4_K1a13_1_I",
                                                       "PT4_K1a13_2_R", "PT4_K1a13_2_I",
                                                       "PT4_K1a31_2_R", "PT4_K1a31_2_I",
                                                       "PT4_K1a13_2_11e_R", "PT4_K1a13_2_11e_I",
                                                       "PT4_K1a13_2_21e_R", "PT4_K1a13_2_21e_I",
                                                       "PT4_K1a13_2_11o_R", "PT4_K1a13_2_11o_I",
                                                       "PT4_K1a13_2_21o_R", "PT4_K1a13_2_21o_I",
                                                       "PT4_K1a13_2_12o_R", "PT4_K1a13_2_12o_I",
                                                       "PT4_K1a13_2_22o_R", "PT4_K1a13_2_22o_I",
                                                       "PT4_K1a13_2_ia_R", "PT4_K1a13_2_ia_I",
                                                       "PT4_K1a13_2_ib_R", "PT4_K1a13_2_ib_I",
                                                       "PT4_K1a13_2_iia_R", "PT4_K1a13_2_iia_I",
                                                       "PT4_K1a13_2_iib_R", "PT4_K1a13_2_iib_I",
                                                       "PT4_K1a13_2_iva_R", "PT4_K1a13_2_iva_I",
                                                       "PT4_K1a13_2_ivb_R", "PT4_K1a13_2_ivb_I"
                                                       },
                                  {bfreqs,
                                   PT2_K1a.vertex[0].avertex.K1.real(), PT2_K1a.vertex[0].avertex.K1.imag(),
                                   PT2_K1p.vertex[0].pvertex.K1.real(), PT2_K1p.vertex[0].pvertex.K1.imag(),
                                   PT2_K1t.vertex[0].tvertex.K1.real(), PT2_K1t.vertex[0].tvertex.K1.imag(),
                                   PT2_SE_a.selfenergy.Sigma.real(), PT2_SE_a.selfenergy.Sigma.imag(),
                                   PT2_SE_p.selfenergy.Sigma.real(), PT2_SE_p.selfenergy.Sigma.imag(),
                                   PT2_SE_t.selfenergy.Sigma.real(), PT2_SE_t.selfenergy.Sigma.imag(),
                                   PT2_SE_p_1.selfenergy.Sigma.real(), PT2_SE_p_1.selfenergy.Sigma.imag(),
                                   PT2_SE_p_4.selfenergy.Sigma.real(), PT2_SE_p_4.selfenergy.Sigma.imag(),
                                   PT2_SE_p_5.selfenergy.Sigma.real(), PT2_SE_p_5.selfenergy.Sigma.imag(),
                                   PT3_K1a.vertex[0].avertex.K1.real(), PT3_K1a.vertex[0].avertex.K1.imag(),
                                   PT3_K2a.vertex[0].avertex.K2.real(), PT3_K2a.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_ia.vertex[0].avertex.K2.real(), PT3_K2a_ia.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_ib.vertex[0].avertex.K2.real(), PT3_K2a_ib.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_iia.vertex[0].avertex.K2.real(), PT3_K2a_iia.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_iib.vertex[0].avertex.K2.real(), PT3_K2a_iib.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_iva.vertex[0].avertex.K2.real(), PT3_K2a_iva.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_ivb.vertex[0].avertex.K2.real(), PT3_K2a_ivb.vertex[0].avertex.K2.imag(),
                                   PT3_K2a_t.vertex[0].avertex.K2.real(), PT3_K2a_t.vertex[0].avertex.K2.imag(),
                                   PT3_K2p.vertex[0].avertex.K2.real(), PT3_K2p.vertex[0].avertex.K2.imag(),
                                   PT3_K2t.vertex[0].avertex.K2.real(), PT3_K2t.vertex[0].avertex.K2.imag(),
                                   PT3_K2t_a.vertex[0].avertex.K2.real(), PT3_K2t_a.vertex[0].avertex.K2.imag(),
                                   PT3_K2t_a.vertex[0].pvertex.K2.real(), PT3_K2t_a.vertex[0].pvertex.K2.imag(),
                                   PT3_K2t_a.vertex[0].tvertex.K2.real(), PT3_K2t_a.vertex[0].tvertex.K2.imag(),
                                   PT3_SE.selfenergy.Sigma.real(), PT3_SE.selfenergy.Sigma.imag(),
                                   PT3_SE_a.selfenergy.Sigma.real(), PT3_SE_a.selfenergy.Sigma.imag(),
                                   PT3_SE_p.selfenergy.Sigma.real(), PT3_SE_p.selfenergy.Sigma.imag(),
                                   PT3_SE_t.selfenergy.Sigma.real(), PT3_SE_t.selfenergy.Sigma.imag(),
                                   PT3_SE_t_1.selfenergy.Sigma.real(), PT3_SE_t_1.selfenergy.Sigma.imag(),
                                   PT3_SE_t_4.selfenergy.Sigma.real(), PT3_SE_t_4.selfenergy.Sigma.imag(),
                                   PT3_SE_t_5.selfenergy.Sigma.real(), PT3_SE_t_5.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a.selfenergy.Sigma.real(), PT3_SE_t_a.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a_1.selfenergy.Sigma.real(), PT3_SE_t_a_1.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a_4.selfenergy.Sigma.real(), PT3_SE_t_a_4.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a_5.selfenergy.Sigma.real(), PT3_SE_t_a_5.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p.selfenergy.Sigma.real(), PT3_SE_t_p.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p_1.selfenergy.Sigma.real(), PT3_SE_t_p_1.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p_4.selfenergy.Sigma.real(), PT3_SE_t_p_4.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p_5.selfenergy.Sigma.real(), PT3_SE_t_p_5.selfenergy.Sigma.imag(),
                                   PT4_K1a22.vertex[0].avertex.K1.real(), PT4_K1a22.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_1.vertex[0].avertex.K1.real(), PT4_K1a13_1.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2.vertex[0].avertex.K1.real(), PT4_K1a13_2.vertex[0].avertex.K1.imag(),
                                   PT4_K1a31_2.vertex[0].avertex.K1.real(), PT4_K1a31_2.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_11e.vertex[0].avertex.K1.real(), PT4_K1a13_2_11e.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_21e.vertex[0].avertex.K1.real(), PT4_K1a13_2_21e.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_11o.vertex[0].avertex.K1.real(), PT4_K1a13_2_11o.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_21o.vertex[0].avertex.K1.real(), PT4_K1a13_2_21o.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_12o.vertex[0].avertex.K1.real(), PT4_K1a13_2_12o.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_22o.vertex[0].avertex.K1.real(), PT4_K1a13_2_22o.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_ia.vertex[0].avertex.K1.real(), PT4_K1a13_2_ia.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_ib.vertex[0].avertex.K1.real(), PT4_K1a13_2_ib.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_iia.vertex[0].avertex.K1.real(), PT4_K1a13_2_iia.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_iib.vertex[0].avertex.K1.real(), PT4_K1a13_2_iib.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_iva.vertex[0].avertex.K1.real(), PT4_K1a13_2_iva.vertex[0].avertex.K1.imag(),
                                   PT4_K1a13_2_ivb.vertex[0].avertex.K1.real(), PT4_K1a13_2_ivb.vertex[0].avertex.K1.imag()});

    //write_hdf("PT4_check_of_K2a_K2_switchedcc_t_update_SE_symm_full_adap_m3m9_gW10_501_101_nI1501_U1_state_PT3_K2a", 0, 1, PT3_K2a);
}

/**
 * Master function to test both consistency and correctness of K2-class
 * @param Lambda
 */
//void test_K2(double Lambda, bool test_consistency){
//
//
//    //First test consistency
//    if(test_consistency) {
//        bool K2a = test_K2_consistency(Lambda, 'a');    //Consistency of a-channel
//        bool K2p = test_K2_consistency(Lambda, 'p');    //Consistency of p-channel
//        bool K2t = test_K2_consistency(Lambda, 't');    //Consistency of t-channel
//
//        if(K2a&&K2p&&K2t)
//            test_K2_correctness(Lambda);
//    }
//
//    test_K2_correctness(Lambda);
//
//}
//#endif

#ifdef STATIC_FEEDBACK
/**
 * Compute the right hand side of the flow equations according to Severin Jakobs' channel decomposition with
 * approximated channel feedback and modified self-energy feedback (only static level shift to avoid
 * overbroadening of spectral features)
 * @param Psi    : state at which to compute right hand side
 * @param Lambda : Lambda at which to compute right hand side
 * @return       : dPsi (right hand side of flow equation)
 */
auto rhs_channel_decomposition(const State<comp>& Psi, const double Lambda) -> State<comp> {
    State<comp> dPsi; // result

    SelfEnergy<comp> selfEnergy;
    comp static_shift = real(Psi.selfenergy.valsmooth(0, glb_mu, 0));  // only use a static level shift as self-energy
    selfEnergy.initialize(static_shift, 0.);

    Propagator G(Lambda, selfEnergy, 'g');    //Initialization of Propagator objects
    Propagator S(Lambda, selfEnergy, 's');    //Initialization of Propagator objects

    // Self-energy flow
    loop(dPsi.selfenergy, Psi.vertex, S, true);  // self-energy loop

    // Vertex flow
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'a', true, '.'); // diff. bubble in the a-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'p', true, '.'); // diff. bubble in the p-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 't', true, '.'); // diff. bubble in the t-channel

    return dPsi;
}

/**
 * FRG flow according to Severin Jakobs' channel decomposition with approximated channel feedback
 * and modified self-energy feedback (only static level shift to avoid overbroadening of spectral features).
 * Only correct if parameter STATIC_FEEDBACK is defined.
 * @param N_ODE : number of Runge-Kutta ODE iterations
 */
void test_channel_decomposition(int N_ODE) {
    State<comp> state_ini, state_fin;   // create initial and final state
    state_ini.initialize();             // initialize initial state

//    sopt_state(state_ini, Lambda_ini);  // set initial state to SOPT (necessary if Lambda_ini is too small)

    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_channel_decomposition,
                    log_substitution, log_resubstitution, N_ODE); // compute flow

    string name = "test_channel_decomposition.h5";
    write_h5_rvecs(name, {"v", "Sigma_re", "Sigma_im", "Sigma_ini_re", "Sigma_ini_im"},
                         {ffreqs,
                          state_fin.selfenergy.Sigma.real(), state_fin.selfenergy.Sigma.imag(),
                          state_ini.selfenergy.Sigma.real(), state_ini.selfenergy.Sigma.imag()});

    write_hdf("channel_decomposition.h5", 0, 1, state_fin);
}
#endif

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
