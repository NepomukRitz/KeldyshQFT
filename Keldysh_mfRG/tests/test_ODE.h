#ifndef FPP_MFRG_TEST_ODE_H
#define FPP_MFRG_TEST_ODE_H


#include "../state.h"                  // State class
#include "../loop.h"                   // self-energy loop
#include "../bubbles.h"                // bubble function
#include "../ODE_solvers.h"                // ODE solvers
#include "../right_hand_sides.h"       // compute the right hand sides of flow equations
#include "../utilities/write_data2file.h"        // writing data to txt or hdf5 file
#include "../utilities/hdf5_routines.h"          // writing states to hdf5 file
#include "../perturbation_theory.h"
#include "../integrator/integrator.h"


/**
 * to be used as rhs in test below --> test_ODE_solvers()
 * @param y
 * @param x
 * @return
 */
double test_rhs_ODE_exp(const double& y, const double x) {
    return y;
}

/**
 * test ODE solvers in solving dy/dx = y from x=0 to x=1 with y(0)=1; solution is y(x)=e^x, y(1)=e;
 */
void test_ODE_solvers() {
    double y_ini, y_fin_Euler, y_fin_RK4, x_ini, x_fin; // necessary variables
    y_ini = 1.; x_ini = 0.; x_fin = 1.; // boundary values
    const int N_ODE_Euler = 100; const int N_ODE_RK4 = 10; // number of steps in ODE solver
    ODE_solver_Euler(y_fin_Euler,  x_fin, y_ini, x_ini, test_rhs_ODE_exp, N_ODE_Euler);
    ODE_solver_RK4(y_fin_RK4,  x_fin, y_ini, x_ini, test_rhs_ODE_exp, N_ODE_RK4);
    cout << "Exact result is " << exp(1.) << "." << endl;
    cout << "Using " << N_ODE_Euler << " steps, Euler gives " << y_fin_Euler << "." << endl;
    cout << "Using " << N_ODE_RK4 << " steps, RK4 gives " << y_fin_RK4 << "." << endl;
}

double x_t_trafo_quadr(const double t) { // variable transformation x = x(t)
    return t*t;
}

double dx_dt_trafo_quadr(const double t) { // variable transformation dx/dt = \dot{x}(t)
    return 2.*t;
}

template <typename T> // evaluate any right hand side with any variable transformation
T rhs_trafo (const T& y, const double t, T rhs (const T& y, const double x), double x_t_trafo(double t), double dx_dt_trafo(double t)) {
    return rhs(y, x_t_trafo(t)) * dx_dt_trafo(t);
}

double test_rhs_ODE_exp_quadr(const double& y, const double t) { // test ODE with quadratic variable transformation
    return rhs_trafo(y, t, test_rhs_ODE_exp, x_t_trafo_quadr, dx_dt_trafo_quadr);
}

void test_ODE_solvers_quadr() { // test ODE solvers in solving dy/dt = y*2t from t=0 to t=1 with y(0)=1; solution is y(x)=e^{t^2}, y(1)=e;
    double y_ini, y_fin_Euler, y_fin_RK4, t_ini, t_fin; // necessary variables
    y_ini = 1.; t_ini = 0.; t_fin = 1.; // boundary values
    const int N_ODE_Euler = 100; const int N_ODE_RK4 = 10; // number of steps in ODE solver
    ODE_solver_Euler(y_fin_Euler,  t_fin, y_ini, t_ini, test_rhs_ODE_exp_quadr, N_ODE_Euler);
    ODE_solver_RK4(y_fin_RK4,  t_fin, y_ini, t_ini, test_rhs_ODE_exp_quadr, N_ODE_RK4);
    cout << "Exact result is " << exp(1.) << "." << endl;
    cout << "Using " << N_ODE_Euler << " steps, Euler gives " << y_fin_Euler << "." << endl;
    cout << "Using " << N_ODE_RK4 << " steps, RK4 gives " << y_fin_RK4 << "." << endl;
}


double test_rhs_SCE_sqrt(const double& y, const double x) {
    return -x/(1.-y);
}

void test_SCE_solver() { // test SCE solvers in solving y=-x/(1-y); solution is y(x)=1/2*(1\pm\sqrt{1+4x}), stable solution: y(x)=1/2*(1-\sqrt{1+4x}), y(1) = 1/2*(1-\sqrt{5})
    double y_ini, y_fin, x; // necessary variables
    y_ini = 0.; x = 1.; // initial y and fixed x
    const int N_SCE = 100; // number of steps in ODE solver
    SCE_solver(y_fin, y_ini, x, test_rhs_SCE_sqrt, N_SCE, 0.);
    cout << "Exact result is " << (1.-sqrt(5.))/2. << ". Using " << N_SCE << " iterations yields " << y_fin << "." << endl;
}



/**
 * Function desgined to test the flow of the vertex in SOPT using the State class
 * @param input     : Input state (not used inside the function)
 * @param Lambda    : Lambda at which to calculate the rhs of the eq.
 * @return          : State carrying in the K1 vertex the results of the computation
 */
template <typename Q>
auto rhs_bubbles_flow_wstate(const State<Q>& input, double Lambda) -> State<Q>{
    State<Q> ans (Lambda);    //Initialize the answer-object

    //Calculating propagator objects of the required types
    Propagator<Q> g(Lambda, input.selfenergy, 'g');
    Propagator<Q> s(Lambda, input.selfenergy, 's');

    bubble_function(ans.vertex, input.vertex, input.vertex, g, s, 'a', true);
    return ans;
}

/**
 * Function to call when testing the rhs of the bubbles flow with the State class
 * @param N_ODE : Number of ODE steps to take between the globally defined Lambda_ini and Lambda_fin
 * @param Lambda_i : initial Lambda value of flow
 * @param Lambda_f : final Lambda value of flow
 * @param write_flag : whether to write output in hdf5
 */
template <typename Q>
void test_rhs_bubbles_flow_wstate(int N_ODE, double Lambda_i, double Lambda_f, bool write_flag = true) {
    State<Q> state_dir (Lambda_f), state_fin (Lambda_f), state_ini (Lambda_i); // direct, final, initial K1a_1
    state_dir.initialize(); // initialize
    state_ini.initialize(); // initialize

    Propagator<Q> G0ini(Lambda_i, state_ini.selfenergy, 'g'); // initial propagator
    Propagator<Q> G0dir(Lambda_f, state_dir.selfenergy, 'g'); // final propagator

    sopt_state(state_ini, Lambda_i); // direct calculation of initial K1a
    sopt_state(state_dir, Lambda_f); // direct calculation of direct K1a

    ODE_solver_RK4(state_fin, Lambda_f, state_ini, Lambda_i, rhs_bubbles_flow_wstate, N_ODE); // final K1a from ODE
    cvec K1a_dif = state_dir.vertex[0].avertex().K1 + ( state_fin.vertex[0].avertex().K1*(-1.) ); // difference in results
    print("Testing ODE for bare K1a_0 with State class. Using " +to_string(N_ODE)+ " ODE steps, the maximal difference between direct and ODE-final result is " +to_string(K1a_dif.max_norm())+ ".", true);
    if(write_flag) write_h5_rvecs("rhs_bubbles_flow_wstate.h5",
                                  {"v", "state_dir_R", "state_dir_I", "state_fin_R", "state_fin_I", "state_ini_R", "state_ini_I"},
                                  {bfreqs, state_dir.vertex[0].avertex().K1.real(), state_dir.vertex[0].avertex().K1.imag(),
                                   state_fin.vertex[0].avertex().K1.real(), state_fin.vertex[0].avertex().K1.imag(),
                                   state_ini.vertex[0].avertex().K1.real(), state_ini.vertex[0].avertex().K1.imag()});
}

/**
 * Function desgined to test the flow of the vertex in SOPT using the only cvecs
 * @param input     : Input cvec (not used inside the function)
 * @param Lambda    : Lambda at which to calculate the rhs of the eq.
 * @return          : The results of the calculation
 */
template <typename Q>
auto rhs_bubbles_flow(const cvec& input, double Lambda) -> cvec{
    cvec ans(nw1_a);   //Initialize the answer

    SelfEnergy<Q> selfini (Lambda);   //Initialize self energy
    selfini.initialize(glb_U/2., 0.);   //Hartree term

    Propagator<Q> g(Lambda, selfini, 'g'); //Regular propagator
    Propagator<Q> s(Lambda, selfini, 's'); //Single-scale propagator

    for(int i=0; i<nBOS; i++){
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble<Q> integrandPia11(g, s, true, w, 11, 'a');     //KA
        IntegrandBubble<Q> integrandPia13(g, s, true, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator<Q>(integrandPia11, glb_w_lower, glb_w_lower);
        auto cont13 = integrator<Q>(integrandPia13, glb_w_lower, glb_w_lower);

        //Add the respective contributions to the respective bubble
        ans[i] = pow(-glb_U/2., 2.)*(cont11+ cont13);            //11+13 = OE => Keldysh comp0
    }

    return ans;
}

/**
 * Function to call when testing the rhs of the bubbles flow with cvecs
 * @param N_ODE : Number of ODE steps to take between the globally defined Lambda_ini and Lambda_fin
 */
template <typename Q>
void test_rhs_bubbles_flow(int N_ODE){
    bool write_flag = true; // whether to write output in hdf5
    vec<Q> K1a_dir(nw1_a), K1a_fin(nw1_a), K1a_ini(nw1_a); // direct, final, initial K1a_1
    SelfEnergy<Q> SEin (Lambda_ini); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<Q> G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    Propagator<Q> G0dir(Lambda_fin, SEin, 'g'); // final propagator

    // direct calculation of initial K1a
    for(int i=0; i<nw1_a; ++i) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble<Q> integrandPia11(G0ini, G0ini, false, w, 11, 'a');     //KA
        IntegrandBubble<Q> integrandPia13(G0ini, G0ini, false, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator<Q>(integrandPia11, glb_w_lower, glb_w_lower);
        auto cont13 = integrator<Q>(integrandPia13, glb_w_lower, glb_w_lower);

        K1a_ini[i] = pow(-glb_U/2.,2.)*(cont11 + cont13);
    }

    // direct calculation of direct K1a
    for(int i=0; i<nw1_a; ++i) {
        double w = bfreqs[i];

        //Create the objects explicitly designed to return the determined Keldysh component needed
        IntegrandBubble<Q> integrandPia11(G0dir, G0dir, false, w, 11, 'a');     //KA
        IntegrandBubble<Q> integrandPia13(G0dir, G0dir, false, w, 13, 'a');     //RK

        //Calculate the contributions
        auto cont11 = integrator<Q>(integrandPia11, glb_w_lower, glb_w_lower);
        auto cont13 = integrator<Q>(integrandPia13, glb_w_lower, glb_w_lower);

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
 * Function to implement the flow of a State in SOPT.
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
template <typename Q>
auto rhs_state_flow_SOPT(const State<Q>& Psi, const double Lambda, const int feedback) -> State<Q>{
    State<Q> dPsi (Lambda);   //Answer object

    State<Q> bare (Lambda);   //Bare state
    bare.initialize();  //Initialize bare state
    Propagator<Q> G(Lambda, bare.selfenergy,'g');    //Initialization of Propagator objects
    Propagator<Q> S(Lambda, bare.selfenergy,'s');    //Initialization of Propagator objects

    if(!(feedback==0 || feedback==3)){  //Check whether Self Energy feedback to Propagators is wanted
        G=Propagator<Q>(Lambda, Psi.selfenergy, 'g');
        S=Propagator<Q>(Lambda, Psi.selfenergy, 's');
    }

    //Self energy loop
    loop(dPsi.selfenergy, Psi.vertex, S, true);  //Loop for the Self-Energy calculation


    if(feedback>=3){    //If feedback>=3, there is vertex feedback
        bare.vertex = Psi.vertex;
    }

    if(feedback==2 || feedback==5) {  //These two options make use of the Katanin substitution
        S=Propagator<Q>(Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
    }

    //Vertex flow
    bubble_function(dPsi.vertex, bare.vertex, bare.vertex, G, S, 'a', true);  //Differentiated bubble in the a-channel
    bubble_function(dPsi.vertex, bare.vertex, bare.vertex, G, S, 'p', true);  //Differentiated bubble in the p-channel
    bubble_function(dPsi.vertex, bare.vertex, bare.vertex, G, S, 't', true);  //Differentiated bubble in the t-channel

    return dPsi;
}

//No feedback
template <typename Q>
auto rhs_state_flow_SOPT_0(const State<Q>& Psi, const double Lambda) -> State<Q>{
    return rhs_state_flow_SOPT(Psi, Lambda, 0);
}

//Self-Energy fed back into Propagators
template <typename Q>
auto rhs_state_flow_SOPT_1(const State<Q>& Psi, const double Lambda) -> State<Q>{
    return rhs_state_flow_SOPT(Psi, Lambda, 1);
}

//As above + Katanin
template <typename Q>
auto rhs_state_flow_SOPT_2(const State<Q>& Psi, const double Lambda) -> State<Q>{
    return rhs_state_flow_SOPT(Psi, Lambda, 2);
}

//Vertex feedback, free propagators
template <typename Q>
auto rhs_state_flow_SOPT_3(const State<Q>& Psi, const double Lambda) -> State<Q>{
    return rhs_state_flow_SOPT(Psi, Lambda, 3);
}

//Self_energy fed back into Propagators + Vertex feedback
template <typename Q>
auto rhs_state_flow_SOPT_4(const State<Q>& Psi, const double Lambda) -> State<Q>{
    return rhs_state_flow_SOPT(Psi, Lambda, 4);
}

//As above + Katanin
template <typename Q>
auto rhs_state_flow_SOPT_5(const State<Q>& Psi, const double Lambda) -> State<Q>{
    return rhs_state_flow_SOPT(Psi, Lambda, 5);
}

/**
 * Function to test the correctness of the flow of the State
 * @param N_ODE : Numbres of ODE-solver steps to be taken
 */
template <typename Q>
void test_rhs_state_flow_SOPT(int N_ODE, int feedback){
    bool write_flag = true; // whether to write output in hdf5
    State<Q> state_dir (Lambda_fin), state_fin (Lambda_fin), state_ini (Lambda_ini); // direct, final, initial K1a_1
    state_dir.initialize(); // initialize state
    state_ini.initialize(); // initialize state

    Propagator<Q> G0ini(Lambda_ini, state_ini.selfenergy, 'g'); // initial propagator
    Propagator<Q> G0dir(Lambda_fin, state_dir.selfenergy, 'g'); // final propagator

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
        K1a0_dif[i] = state_dir.vertex[0].avertex().K1_val(0, i, 0) - state_fin.vertex[0].avertex().K1_val(0, i, 0);
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
                                   state_dir.vertex[0].avertex().K1.real(), state_dir.vertex[0].avertex().K1.imag(),
                                   state_fin.vertex[0].avertex().K1.real(), state_fin.vertex[0].avertex().K1.imag(),
                                   state_ini.vertex[0].avertex().K1.real(), state_ini.vertex[0].avertex().K1.imag(),
                                   state_dir.vertex[0].pvertex().K1.real(), state_dir.vertex[0].pvertex().K1.imag(),
                                   state_fin.vertex[0].pvertex().K1.real(), state_fin.vertex[0].pvertex().K1.imag(),
                                   state_ini.vertex[0].pvertex().K1.real(), state_ini.vertex[0].pvertex().K1.imag(),
                                   state_dir.vertex[0].tvertex().K1.real(), state_dir.vertex[0].tvertex().K1.imag(),
                                   state_fin.vertex[0].tvertex().K1.real(), state_fin.vertex[0].tvertex().K1.imag(),
                                   state_ini.vertex[0].tvertex().K1.real(), state_ini.vertex[0].tvertex().K1.imag()});

}

// compute differentiated K1a non-ladder diagram in PT4, used to compute flow of this specific diagram
// (as a check for GeneralVertex class)
auto rhs_PT4_K1a_nonladder_flow(const State<state_datatype>& Psi, const double Lambda) -> State<state_datatype> {
    State<state_datatype> bare (Lambda); // bare state
    bare.initialize();         // initialize bare state

    Propagator<state_datatype> G (Lambda, bare.selfenergy, 'g'); // bare propagator
    Propagator<state_datatype> S (Lambda, bare.selfenergy, 's'); // bare single-scale propagator

    // compute K1p in PT2
    State<state_datatype> PT2_K1p (Lambda);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false);

    // compute K2a in PT3
    State<state_datatype> PT3_K2a (Lambda);
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'a', false);   // K2a in PT3 (PT2_K1t = 0)

    // contribution to \dot{K1a} where leftmost bubble is differentiated
    State<state_datatype> PT4_K1a_dot_left (Lambda);
    bubble_function(PT4_K1a_dot_left.vertex, bare.vertex, PT3_K2a.vertex, G, S, 'a', true);

    // compute differentiated K2a in PT3
    State<state_datatype> PT3_K2a_dot_right (Lambda);
    bubble_function(PT3_K2a_dot_right.vertex, PT2_K1p.vertex, bare.vertex, G, S, 'a', true);   // \dot{K2a} in PT3 (PT2_K1t = 0)

    // contribution to \dot{K1a} where rightmost bubble is differentiated
    State<state_datatype> PT4_K1a_dot_right (Lambda);
    bubble_function(PT4_K1a_dot_right.vertex, bare.vertex, PT3_K2a_dot_right.vertex, G, G, 'a', false);

    // compute differentiated K1p in PT2
    State<state_datatype> PT2_K1p_dot (Lambda);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, S, 'p', true);

    // compute K2a/K2ab with differentiated p bubble in PT3, as an ingredient of the contribution to \dot{K1a} in PT4
    // with differentiated center bubble
    State<state_datatype> PT3_K2a_dot_half1 (Lambda);
    State<state_datatype> PT3_K2a_dot_half2 (Lambda);
    bubble_function(PT3_K2a_dot_half1.vertex, PT2_K1p_dot.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT3_K2a_dot_half2.vertex, bare.vertex, PT2_K1p_dot.vertex, G, G, 'a', false);

    // construct non-symmetric K2a with differentiated p bubble
    GeneralVertex<state_datatype, non_symmetric> PT3_K2a_dot (n_spin);
    PT3_K2a_dot[0].half1() = PT3_K2a_dot_half1.vertex[0].half1();
    PT3_K2a_dot[0].half2() = PT3_K2a_dot_half2.vertex[0].half1();

    // construct non-symmetric K2ab with differentiated p bubble
    GeneralVertex<state_datatype, non_symmetric> PT3_K2ab_dot (n_spin);
    PT3_K2ab_dot[0].half1() = PT3_K2a_dot_half2.vertex[0].half1();
    PT3_K2ab_dot[0].half2() = PT3_K2a_dot_half1.vertex[0].half1();

    // contribution to \dot{K1a} where center bubble is differentiated
    State<state_datatype> PT4_K1a_dot_center (Lambda);
    State<state_datatype> PT4_K1a_dot_center_left (Lambda);
    State<state_datatype> PT4_K1a_dot_center_right (Lambda);

    bubble_function(PT4_K1a_dot_center_left.vertex, bare.vertex, PT3_K2a_dot, G, G, 'a', false);
    bubble_function(PT4_K1a_dot_center_right.vertex, PT3_K2ab_dot, bare.vertex, G, G, 'a', false);
    PT4_K1a_dot_center = (PT4_K1a_dot_center_left + PT4_K1a_dot_center_right) * 0.5;

    // return sum of contributions with left, center, and right differentiated bubble
    return PT4_K1a_dot_left + PT4_K1a_dot_center + PT4_K1a_dot_right;
}


// compute K1a non-ladder diagram in PT4 (to compare it to its flow)
auto compute_PT4_K1a_nonladder(const double Lambda) -> State<state_datatype> {
    State<state_datatype> bare (Lambda); // bare state
    bare.initialize();         // initialize bare state

    Propagator<state_datatype> G (Lambda, bare.selfenergy, 'g'); // bare propagator

    // K1p in PT2
    State<state_datatype> PT2_K1p (Lambda);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false);

    // K2a in PT3
    State<state_datatype> PT3_K2a (Lambda);
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'a', false);

    // K1a nonladder in PT4
    State<state_datatype> PT4_K1a_nonladder (Lambda);
    bubble_function(PT4_K1a_nonladder.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false);

    return PT4_K1a_nonladder;
}

// compute K1a non-ladder diagram in PT4 directly and via its flow, as a test for GeneralVertex class
void test_PT4_K1a_nonladder_flow(const double Lambda_i, const double Lambda_f, const string filename) {

    // number of Lambda layers in hdf5 file
    int Lambda_size = nODE + U_NRG.size() + 1;

    // compute PT4 K1a non-ladder at initial Lambda
    State<state_datatype> state_ini = compute_PT4_K1a_nonladder(Lambda_i);
    // save result to 0-th layer in hdf5 file (for both flow and direct computation)
    write_hdf(filename, 0, Lambda_size, state_ini);
    write_hdf(filename + "_dir", 0, Lambda_size, state_ini);

    // generate Lambda grid
    rvec Lambdas = construct_flow_grid(Lambda_f, Lambda_i, sq_substitution, sq_resubstitution, nODE);

    for (int i=1; i<Lambdas.size(); ++i) {
        // compute direct result of PT4 K1a non-ladder at each Lambda step
        State<state_datatype> state_dir = compute_PT4_K1a_nonladder(Lambdas[i]);
        // save result to last layer in hdf5 file
        add_hdf(filename + "_dir", i, Lambda_size, state_dir, Lambdas);
    }

    // compute flow of PT4 K1a non-ladder from initial to final Lambda
    State<state_datatype> state_fin (Lambda_f);
    ODE_solver_RK4(state_fin, Lambda_f, state_ini, Lambda_i, rhs_PT4_K1a_nonladder_flow,
                   sq_substitution, sq_resubstitution, nODE, filename);
}



#endif //FPP_MFRG_TEST_ODE_H
