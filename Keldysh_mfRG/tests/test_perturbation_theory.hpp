#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

#include <cmath>                    // use M_PI as pi
#include "../data_structures.hpp"
#include "../correlation_functions/two_point/propagator.hpp"                // propagators
#include "../correlation_functions/state.hpp"                  // State class
#include "../loop/loop.hpp"                   // self-energy loop
#include "../bubble/bubble_function.hpp"                // bubble function
#include "../ODE_solvers/ODE_solvers.hpp"               // ODE solvers
#include "../mfRG_flow/right_hand_sides.hpp"       // compute the right hand sides of flow equations
#include "../utilities/write_data2file.hpp"        // writing data to txt or hdf5 file
#include "../utilities/hdf5_routines.hpp"          // writing states to hdf5 file
#include "../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include <boost/math/special_functions/polygamma.hpp> // Polygamma function
#include "../integrator/integrator.hpp"
#include "../postprocessing/causality_FDT_checks.hpp"   // check causality and FDTs
#include "../utilities/util.hpp"
#include "../perturbation_theory_and_parquet/hartree_term.hpp"


namespace {
    fRG_config stdConfig;
}


/**
 * Function to test the loop function and the calculation of the SelfEnergy in SOPT
 * compute the selfenergy with an a-bubble and with a p-bubble and compare --> Are they identical?
 * @param state : State initialized with initial conditions
 */
template <typename Q>
void testSelfEnergy_and_K1(double Lambda){
    SelfEnergy<Q> SE_via_a(Lambda);
    SelfEnergy<Q> SE_via_p(Lambda);
    SelfEnergy<Q> SE_bare(Lambda);
    SE_via_a.initialize(glb_U*0.5, 0.);
    SE_via_p.initialize(glb_U*0.5, 0.);
    SE_bare.initialize(glb_U*0.5, 0.);
    Propagator<Q> g(Lambda, SE_bare, 'g');

    //Calculate the vertex
    State<Q> SOPT_state(Lambda);
    SOPT_state.initialize();
    sopt_state(SOPT_state, Lambda);

    Vertex<Q,false> temp_vertex_a (Lambda), temp_vertex_p (Lambda); //All zeros
    temp_vertex_a.initialize((KELDYSH_FORMALISM and !CONTOUR_BASIS) ? -glb_U*0.5 : -glb_U);
    temp_vertex_p.initialize((KELDYSH_FORMALISM and !CONTOUR_BASIS) ? -glb_U*0.5 : -glb_U);
    temp_vertex_a.avertex() = SOPT_state.vertex.avertex();
    temp_vertex_p.pvertex() = SOPT_state.vertex.pvertex();

    loop<false,0>(SE_via_a, temp_vertex_a, g);//Calculate the SelfEnergy in SOPT
    loop<false,0>(SE_via_p, temp_vertex_p, g);//Calculate the SelfEnergy in SOPT

    SelfEnergy<Q> SE_diff = SE_via_a - SE_via_p;

    double diff_max = SE_diff.Sigma.get_vec().max_norm();
    double diff_max_vertex = (temp_vertex_a.avertex().K1.get_vec() + temp_vertex_p.pvertex().K1.get_vec()).max_norm();
    utils::print("Maximal difference between SOPT K1 via a-bubble and p-bubble: \t", diff_max_vertex, " \n");
    utils::print("Maximal difference between SOPT selfenergy via a-bubble and p-bubble: \t", diff_max, " \n");

    if (false) {
        /// output both results in HDF5 file
        H5::H5File file_out(data_dir + "comparison_SOPT_via_a_p_channel.h5", H5F_ACC_TRUNC);
        write_to_hdf(file_out, "SOPT_via_a", SE_via_a.Sigma.get_vec(), false);
        write_to_hdf(file_out, "SOPT_via_p", SE_via_p.Sigma.get_vec(), false);
        file_out.close();

    }


}


/**
 * Function that prints out a .h5 file with the value of the rhs of a SOPT flow at the given Lambda for both an FFT and an fRG calculation
 * @param Lambda    : Lambda at which the derivatives are to be calculated
 */
template <typename Q>
void test_derivatives_K1a(double Lambda){
    vec<Q> blah(nw1_a);
    vec<Q> rhs_SOPT_FFT_K1a = dSOPT_FFT_K1a_rhs(blah, Lambda);
    vec<Q> rhs_flow = rhs_bubbles_flow(blah, Lambda);

    write_h5_rvecs("derivatives_K1a.h5",
                   {"v", "FFT_R", "FFT_I", "SOPT_R", "SOPT_I"},
                   {ffreqs, rhs_SOPT_FFT_K1a.real(), rhs_SOPT_FFT_K1a.imag(), rhs_flow.real(), rhs_flow.imag()});
}

/**
 * This function checks the consistency of the K2 class
 * @param Lambda    : Scale
 * @param r         : Channel to be tested i.e. K_2^r
 */
template <typename Q>
auto test_K2_consistency(double Lambda, const char r) -> bool{
    State<Q> bare (Lambda);   //Create a bare state
    bare.initialize();  //Initialize bare state

    Propagator<Q> G(Lambda, bare.selfenergy, 'g'); //Bare propagator at scale Lambda>

    //Create and initialize states to save the channel contributions in
    State<Q> K1a (Lambda);
    State<Q> K1p (Lambda);
    State<Q> K1t (Lambda);
    K1a.initialize();
    K1p.initialize();
    K1t.initialize();

    //Calculate K1a, K1p and K1t contributions separately
    bubble_function(K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false, stdConfig);

    //Create and initialize the K2r-objects to test
    State<Q> test_K2r_with_K1a (Lambda);
    State<Q> test_K2r_with_K1p (Lambda);
    State<Q> test_K2r_with_K1t (Lambda);
    test_K2r_with_K1a.initialize();
    test_K2r_with_K1p.initialize();
    test_K2r_with_K1t.initialize();


    if(r=='p' ||  r=='t') {
        //Perform TOPT calculation of K2r with K1a
        bubble_function(test_K2r_with_K1a.vertex, K1a.vertex, bare.vertex, G, G, r, false, stdConfig);
    }

    if(r=='a' || r=='t') {
        //Perform TOPT calculation of K2r with K1p
        bubble_function(test_K2r_with_K1p.vertex, K1p.vertex, bare.vertex, G, G, r, false, stdConfig);
    }

    if(r=='a' || r=='p') {
        //Perform TOPT calculation of K2r with K1t
        bubble_function(test_K2r_with_K1t.vertex, K1t.vertex, bare.vertex, G, G, r, false, stdConfig);
    }


    //TOPT-consistency checks for K2

    bool empty = true;  //Boolean that stays true if everything is zero
    if(r=='a' || r=='p') {
#pragma omp parallel
        //Check parallelized that everything in the K2 vertex is zero.
        for (int index = 0; index < test_K2r_with_K1t.vertex.avertex().K2.size(); ++index) {
            if (test_K2r_with_K1t.vertex.avertex().K2_acc(index) != 0.) {
                empty = false;  //If any one element is not zero, K2 is not empty
                break;  //Exit the for-loop early
            }
        }
        //Print result of consistency check
        if (empty) utils::print("TOPT-consistency check passed. K2a or K2p with K1t is zero everywhere.", true);
        else utils::print("TOPT-consistency check failed. K2a or K2p with K1t is not zero everywhere.", true);
    }
    else if(r=='t'){
#pragma omp parallel
        //Check parallelized that everything in the K2 vertices is zero.
        for (int index = 0; index < test_K2r_with_K1a.vertex.avertex().K2.size(); ++index) {
            if (test_K2r_with_K1p.vertex.tvertex().K2_acc(index) != 0.) {
                empty = false;  //If any one element is not zero, K2 is not empty
                break;  //Exit the for-loop early
            }
        }
        //Print result of consistency check
        if (empty) utils::print("TOPT-consistency check passed. K2t with K1a and K1p are zero everywhere.", true);
        else utils::print("TOPT-consistency check failed. K2t with K1a and K1p are not zero everywhere.", true);
    }

    return empty;
}


/**
 * Function that computes K1 (and K2, K3) up to PT4, and performs FDT checks
 */
void test_PT4(double Lambda, bool write_flag = false) {
    utils::print("Compute K1 (and K2, K3) up to PT4.", true);
    const int it_spin = 0;
    // Initialize a bare state
    State<state_datatype> bare (Lambda);
    bare.initialize();

    // Initialize a bare propagator
    Propagator<state_datatype> G(Lambda, bare.selfenergy, 'g');

    // Compute K1 in PT2
    State<state_datatype> PT2_K1a (Lambda);
    State<state_datatype> PT2_K1p (Lambda);
    State<state_datatype> PT2_K1t (Lambda);

    double t0 = utils::get_time();
    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT2_K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false, stdConfig);
    utils::print("Computed K1 in PT2.", true);
    utils::get_time(t0);
    if (write_flag) {
        write_state_to_hdf(data_dir + "PT2_K1a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1a);
        write_state_to_hdf(data_dir + "PT2_K1p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1p);
        write_state_to_hdf(data_dir + "PT2_K1t_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1t);
    }

    // Compute K1 in PT3, using K1 in PT2
    State<state_datatype> PT3_K1a (Lambda);
    State<state_datatype> PT3_K1p (Lambda);
    State<state_datatype> PT3_K1t (Lambda);

    t0 = utils::get_time();
    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT3_K1p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    // for K1t in PT3, need a-vertex in PT2 due to a <-> t symmetry
    bubble_function(PT3_K1t.vertex, PT2_K1t.vertex + PT2_K1a.vertex, bare.vertex, G, G, 't', false, stdConfig);
    vec<state_datatype> zerosK2 (PT3_K1t.vertex.tvertex().K2.get_vec().size());
    if (MAX_DIAG_CLASS >= 2) PT3_K1t.vertex.tvertex().K2.set_vec(zerosK2); // set K2 part of this vertex to zero
    utils::print("Computed K1 in PT3.", true);
    utils::get_time(t0);
    if (write_flag) {
        write_state_to_hdf(data_dir + "PT3_K1a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K1a);
        write_state_to_hdf(data_dir + "PT3_K1p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K1p);
        write_state_to_hdf(data_dir + "PT3_K1t_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K1t);
    }

    // Compute K2 in PT3, using K1p, K1t in PT2
    State<state_datatype> PT3_K2a (Lambda);
    State<state_datatype> PT3_K2p (Lambda);
    State<state_datatype> PT3_K2t (Lambda);
    State<state_datatype> PT3_K2t_a (Lambda);
    State<state_datatype> PT3_K2t_p (Lambda);

    t0 = utils::get_time();
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'a', false, stdConfig);   // K2a in PT3 (PT2_K1t = 0)
    bubble_function(PT3_K2p.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'p', false, stdConfig);   // K2p in PT3 (PT2_K1t = 0)
    bubble_function(PT3_K2t_a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 't', false, stdConfig);   // contribution of K2t in PT3 obtained by inserting K1a in PT2
    bubble_function(PT3_K2t_p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 't', false, stdConfig);   // contribution of K2t in PT3 obtained by inserting K1p in PT2
    // PT3_K2t_a should also have a K1-contribution due to a-t symmetry (PT2_K1t implicitly inserted) --> set to zero
    vec<state_datatype> zerosK1 (PT3_K2t_a.vertex.tvertex().K1.get_vec().size());
    PT3_K2t_a.vertex.tvertex().K1.set_vec(zerosK1);
    PT3_K2t.vertex = PT3_K2t_a.vertex + PT3_K2t_p.vertex; // sum of contributions from a- and p-insertions

    // K2' in PT3 would be obtained by flipping the left and right vertex, but since K2' is not saved, these terms would give zero
    utils::print("Computed K2 in PT3.", true);
    utils::get_time(t0);
    if (write_flag) {
        write_state_to_hdf(data_dir + "PT3_K2a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2a);
        write_state_to_hdf(data_dir + "PT3_K2p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2p);
        write_state_to_hdf(data_dir + "PT3_K2t_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2t);
        write_state_to_hdf(data_dir + "PT3_K2t_a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2t_a);
        write_state_to_hdf(data_dir + "PT3_K2t_p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2t_p);
    }


    // Compute K3 in PT4, using K1 and K2 in PT2 and PT3
    State<state_datatype> PT4_31 (Lambda);
    State<state_datatype> PT4_31_a_a1 (Lambda);
    State<state_datatype> PT4_31_a_p1 (Lambda);
    State<state_datatype> PT4_31_a_t1 (Lambda);
    State<state_datatype> PT4_31_a_a2 (Lambda);
    State<state_datatype> PT4_31_a_p2 (Lambda);
    State<state_datatype> PT4_31_a_t2 (Lambda);
    State<state_datatype> PT4_31_p_a1 (Lambda);
    State<state_datatype> PT4_31_p_p1 (Lambda);
    State<state_datatype> PT4_31_p_t1 (Lambda);
    State<state_datatype> PT4_31_p_a2 (Lambda);
    State<state_datatype> PT4_31_p_p2 (Lambda);
    State<state_datatype> PT4_31_p_t2 (Lambda);
    State<state_datatype> PT4_31_t_a1 (Lambda);
    State<state_datatype> PT4_31_t_p1 (Lambda);
    State<state_datatype> PT4_31_t_t1 (Lambda);
    State<state_datatype> PT4_31_t_a2 (Lambda);
    State<state_datatype> PT4_31_t_p2 (Lambda);
    State<state_datatype> PT4_31_t_t2 (Lambda);

    State<state_datatype> PT4_13 (Lambda);
    State<state_datatype> PT4_13_a_a1 (Lambda);
    State<state_datatype> PT4_13_a_p1 (Lambda);
    State<state_datatype> PT4_13_a_t1 (Lambda);
    State<state_datatype> PT4_13_a_a2 (Lambda);
    State<state_datatype> PT4_13_a_p2 (Lambda);
    State<state_datatype> PT4_13_a_t2 (Lambda);
    State<state_datatype> PT4_13_p_a1 (Lambda);
    State<state_datatype> PT4_13_p_p1 (Lambda);
    State<state_datatype> PT4_13_p_t1 (Lambda);
    State<state_datatype> PT4_13_p_a2 (Lambda);
    State<state_datatype> PT4_13_p_p2 (Lambda);
    State<state_datatype> PT4_13_p_t2 (Lambda);
    State<state_datatype> PT4_13_t_a1 (Lambda);
    State<state_datatype> PT4_13_t_p1 (Lambda);
    State<state_datatype> PT4_13_t_t1 (Lambda);
    State<state_datatype> PT4_13_t_a2 (Lambda);
    State<state_datatype> PT4_13_t_p2 (Lambda);
    State<state_datatype> PT4_13_t_t2 (Lambda);

    State<state_datatype> PT4_22 (Lambda);
    State<state_datatype> PT4_22_a_aa (Lambda);
    State<state_datatype> PT4_22_a_ap (Lambda);
    State<state_datatype> PT4_22_a_pa (Lambda);
    State<state_datatype> PT4_22_a_pp (Lambda);
    State<state_datatype> PT4_22_p_aa (Lambda);
    State<state_datatype> PT4_22_p_ap (Lambda);
    State<state_datatype> PT4_22_p_pa (Lambda);
    State<state_datatype> PT4_22_p_pp (Lambda);
    State<state_datatype> PT4_22_t_aa (Lambda);
    State<state_datatype> PT4_22_t_ap (Lambda);
    State<state_datatype> PT4_22_t_pa (Lambda);
    State<state_datatype> PT4_22_t_pp (Lambda);

    t0 = utils::get_time();

    // a-channel:
    // 3-1: insert all possible PT3 vertices on the left, bare vertex on the right
    bubble_function(PT4_31_a_a1.vertex, PT3_K1a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_31_a_p1.vertex, PT3_K1p.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_31_a_t1.vertex, PT3_K1t.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_31_a_a2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_31_a_p2.vertex, PT3_K2p.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_31_a_t2.vertex, PT3_K2t.vertex, bare.vertex, G, G, 'a', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_31_a_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_a1);
        write_state_to_hdf(data_dir + "PT4_31_a_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_p1);
        write_state_to_hdf(data_dir + "PT4_31_a_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_t1);
        write_state_to_hdf(data_dir + "PT4_31_a_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_a2);
        write_state_to_hdf(data_dir + "PT4_31_a_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_p2);
        write_state_to_hdf(data_dir + "PT4_31_a_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_31.vertex.avertex() = PT4_31_a_a1.vertex.avertex()
                             + PT4_31_a_p1.vertex.avertex()
                             + PT4_31_a_t1.vertex.avertex()
                             + PT4_31_a_a2.vertex.avertex()
                             + PT4_31_a_p2.vertex.avertex()
                             + PT4_31_a_t2.vertex.avertex();

    // 1-3: insert bare vertex on the left, all possible PT3 vertices on the right
    // this can only give a K1 contribution, since we do not compute K2'
    bubble_function(PT4_13_a_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false, stdConfig);
    // the following should all give zero
    bubble_function(PT4_13_a_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_13_a_t1.vertex, bare.vertex, PT3_K1t.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_13_a_a2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_13_a_p2.vertex, bare.vertex, PT3_K2p.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_13_a_t2.vertex, bare.vertex, PT3_K2t.vertex, G, G, 'a', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_13_a_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_a1);
        write_state_to_hdf(data_dir + "PT4_13_a_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_p1);
        write_state_to_hdf(data_dir + "PT4_13_a_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_t1);
        write_state_to_hdf(data_dir + "PT4_13_a_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_a2);
        write_state_to_hdf(data_dir + "PT4_13_a_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_p2);
        write_state_to_hdf(data_dir + "PT4_13_a_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_13.vertex.avertex() = PT4_13_a_a1.vertex.avertex()
                             + PT4_13_a_p1.vertex.avertex()
                             + PT4_13_a_t1.vertex.avertex()
                             + PT4_13_a_a2.vertex.avertex()
                             + PT4_13_a_p2.vertex.avertex()
                             + PT4_13_a_t2.vertex.avertex();

    // 2-2: insert all possible PT2 vertices on the left and right (PT2_K1t is always zero -> neglect it)
    bubble_function(PT4_22_a_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_22_a_ap.vertex, PT2_K1a.vertex, PT2_K1p.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_22_a_pa.vertex, PT2_K1p.vertex, PT2_K1a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_22_a_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 'a', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_22_a_aa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_aa);
        write_state_to_hdf(data_dir + "PT4_22_a_ap_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_ap);
        write_state_to_hdf(data_dir + "PT4_22_a_pa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_pa);
        write_state_to_hdf(data_dir + "PT4_22_a_pp_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_pp);
    }

    // sum of contributions obtained from different insertions
    PT4_22.vertex.avertex() = PT4_22_a_aa.vertex.avertex()
                             + PT4_22_a_ap.vertex.avertex()
                             + PT4_22_a_pa.vertex.avertex()
                             + PT4_22_a_pp.vertex.avertex();

    utils::print("Computed a-channel in PT4.", true);

    // p-channel:
    // 3-1: insert all possible PT3 vertices on the left, bare vertex on the right
    bubble_function(PT4_31_p_a1.vertex, PT3_K1a.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_31_p_p1.vertex, PT3_K1p.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_31_p_t1.vertex, PT3_K1t.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_31_p_a2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_31_p_p2.vertex, PT3_K2p.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_31_p_t2.vertex, PT3_K2t.vertex, bare.vertex, G, G, 'p', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_31_p_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_a1);
        write_state_to_hdf(data_dir + "PT4_31_p_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_p1);
        write_state_to_hdf(data_dir + "PT4_31_p_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_t1);
        write_state_to_hdf(data_dir + "PT4_31_p_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_a2);
        write_state_to_hdf(data_dir + "PT4_31_p_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_p2);
        write_state_to_hdf(data_dir + "PT4_31_p_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_31.vertex.pvertex() = PT4_31_p_a1.vertex.pvertex()
                             + PT4_31_p_p1.vertex.pvertex()
                             + PT4_31_p_t1.vertex.pvertex()
                             + PT4_31_p_a2.vertex.pvertex()
                             + PT4_31_p_p2.vertex.pvertex()
                             + PT4_31_p_t2.vertex.pvertex();

    // 1-3: insert bare vertex on the left, all possible PT3 vertices on the right
    // this can only give a K1 contribution, since we do not compute K2'
    bubble_function(PT4_13_p_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 'p', false, stdConfig);
    // the following should all give zero
    bubble_function(PT4_13_p_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_13_p_t1.vertex, bare.vertex, PT3_K1t.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_13_p_a2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_13_p_p2.vertex, bare.vertex, PT3_K2p.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_13_p_t2.vertex, bare.vertex, PT3_K2t.vertex, G, G, 'p', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_13_p_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_a1);
        write_state_to_hdf(data_dir + "PT4_13_p_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_p1);
        write_state_to_hdf(data_dir + "PT4_13_p_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_t1);
        write_state_to_hdf(data_dir + "PT4_13_p_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_a2);
        write_state_to_hdf(data_dir + "PT4_13_p_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_p2);
        write_state_to_hdf(data_dir + "PT4_13_p_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_13.vertex.pvertex() = PT4_13_p_a1.vertex.pvertex()
                             + PT4_13_p_p1.vertex.pvertex()
                             + PT4_13_p_t1.vertex.pvertex()
                             + PT4_13_p_a2.vertex.pvertex()
                             + PT4_13_p_p2.vertex.pvertex()
                             + PT4_13_p_t2.vertex.pvertex();

    // 2-2: insert all possible PT2 vertices on the left and right (PT2_K1t is always zero -> neglect it)
    bubble_function(PT4_22_p_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_22_p_ap.vertex, PT2_K1a.vertex, PT2_K1p.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_22_p_pa.vertex, PT2_K1p.vertex, PT2_K1a.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_22_p_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 'p', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_22_p_aa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_aa);
        write_state_to_hdf(data_dir + "PT4_22_p_ap_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_ap);
        write_state_to_hdf(data_dir + "PT4_22_p_pa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_pa);
        write_state_to_hdf(data_dir + "PT4_22_p_pp_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_pp);
    }

    // sum of contributions obtained from different insertions
    PT4_22.vertex.pvertex() = PT4_22_p_aa.vertex.pvertex()
                             + PT4_22_p_ap.vertex.pvertex()
                             + PT4_22_p_pa.vertex.pvertex()
                             + PT4_22_p_pp.vertex.pvertex();

    utils::print("Computed p-channel in PT4.", true);

    // t-channel:
    // in the t-channel, we need to insert a and t simultaneously due to a <-> t symmetry // TODO: remove?
    // (spin sum in the t-channel makes use of this symmetry)                             // TODO: remove?

    // 3-1: insert all possible PT3 vertices on the left, bare vertex on the right
    bubble_function(PT4_31_t_a1.vertex, PT3_K1a.vertex, bare.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_31_t_p1.vertex, PT3_K1p.vertex, bare.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_31_t_t1.vertex, PT3_K1t.vertex, bare.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_31_t_a2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_31_t_p2.vertex, PT3_K2p.vertex, bare.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_31_t_t2.vertex, PT3_K2t.vertex, bare.vertex, G, G, 't', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_31_t_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_a1);
        write_state_to_hdf(data_dir + "PT4_31_t_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_p1);
        write_state_to_hdf(data_dir + "PT4_31_t_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_t1);
        write_state_to_hdf(data_dir + "PT4_31_t_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_a2);
        write_state_to_hdf(data_dir + "PT4_31_t_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_p2);
        write_state_to_hdf(data_dir + "PT4_31_t_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_31.vertex.tvertex() = PT4_31_t_a1.vertex.tvertex()
                             + PT4_31_t_p1.vertex.tvertex()
                             + PT4_31_t_t1.vertex.tvertex()
                             + PT4_31_t_a2.vertex.tvertex()
                             + PT4_31_t_p2.vertex.tvertex()
                             + PT4_31_t_t2.vertex.tvertex();

    // 1-3: insert bare vertex on the left, all possible PT3 vertices on the right
    bubble_function(PT4_13_t_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_13_t_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_13_t_t1.vertex, bare.vertex, PT3_K1t.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_13_t_a2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_13_t_p2.vertex, bare.vertex, PT3_K2p.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_13_t_t2.vertex, bare.vertex, PT3_K2t.vertex, G, G, 't', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_13_t_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_a1);
        write_state_to_hdf(data_dir + "PT4_13_t_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_p1);
        write_state_to_hdf(data_dir + "PT4_13_t_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_t1);
        write_state_to_hdf(data_dir + "PT4_13_t_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_a2);
        write_state_to_hdf(data_dir + "PT4_13_t_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_p2);
        write_state_to_hdf(data_dir + "PT4_13_t_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_13.vertex.tvertex() = PT4_13_t_a1.vertex.tvertex()
                             + PT4_13_t_p1.vertex.tvertex()
                             + PT4_13_t_t1.vertex.tvertex()
                             + PT4_13_t_a2.vertex.tvertex()
                             + PT4_13_t_p2.vertex.tvertex()
                             + PT4_13_t_t2.vertex.tvertex();

    // 2-2: insert all possible PT2 vertices on the left and right (PT2_K1t is always zero -> neglect it)
    bubble_function(PT4_22_t_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_22_t_ap.vertex, PT2_K1a.vertex, PT2_K1p.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_22_t_pa.vertex, PT2_K1p.vertex, PT2_K1a.vertex, G, G, 't', false, stdConfig);
    bubble_function(PT4_22_t_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 't', false, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_22_t_aa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_aa);
        write_state_to_hdf(data_dir + "PT4_22_t_ap_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_ap);
        write_state_to_hdf(data_dir + "PT4_22_t_pa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_pa);
        write_state_to_hdf(data_dir + "PT4_22_t_pp_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_pp);
    }

    // sum of contributions obtained from different insertions
    PT4_22.vertex.tvertex() = PT4_22_t_aa.vertex.tvertex()
                             + PT4_22_t_ap.vertex.tvertex()
                             + PT4_22_t_pa.vertex.tvertex()
                             + PT4_22_t_pp.vertex.tvertex();


    utils::print("Computed t-channel in PT4.", true);

    if (write_flag) {
        write_state_to_hdf(data_dir + "PT4_31_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31);
        write_state_to_hdf(data_dir + "PT4_13_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13);
        write_state_to_hdf(data_dir + "PT4_22_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22);
    }

    utils::print("Computed K1, K2, K3 in PT4.", true);
    utils::get_time(t0);

    /** Make automated checks of all diagrams: Compute the values of all diagrams at all frequencies equal to zero,
     * and for all pairs of diagrams that should cancel, compute the relative deviation of their sum from zero.
     * Print all results to log.
     * */

    utils::print("--- CHECK RESULTS: ---", true);
    utils::print("--- print relative error of quantities that should be zero ---", true);

#if KELDYSH_FORMALISM
    int iK = 1;
#else
    int iK = 0;
#endif
    // input variables: all frequencies equal to zero
    VertexInput input_a (iK, it_spin, 0., 0., 0., 0, 'a');
    VertexInput input_p (iK, it_spin, 0., 0., 0., 0, 'p');
    VertexInput input_t (iK, it_spin, 0., 0., 0., 0, 't');

#if KELDYSH_FORMALISM
    vec<int> iK2s = {1, 2, 4}; // Keldysh indices of fully retarded components of K2
#else
    vec<int> iK2s = {0}; // Keldysh indices of Matsubara component of K2
#endif

    // labels to be printed to log
    std::string check_labels[] {"PT2: K1a + K1p: ", "PT2: K1t: ",
                           "PT2: K1a - exact: ", "PT2: K1p - exact: ",
                           "PT3: K1a - exact: ", "PT3: K1p - exact: ", "PT3: K1t - exact: ",
                           "PT3: K2a[1] - exact: ", "PT3: K2p[1] - exact: ", "PT3: K2t[1] - exact: ",
#if KELDYSH_FORMALISM
                           "PT3: K2a[2] - exact: ", "PT3: K2p[2] - exact: ", "PT3: K2t[2] - exact: ",
                           "PT3: K2a[4] - exact: ", "PT3: K2p[4] - exact: ", "PT3: K2t[4] - exact: "
                            ,
#endif
                           "PT4: (K2a[1] <- K1p) + (K2p[1] <- K1a): ",
#if KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K1p) + (K2p[2] <- K1a): ",
                           "PT4: (K2a[4] <- K1p) + (K2p[4] <- K1a): ",
#endif
                           "PT4: (K2a[1] <- K1t) + (K2p[1] <- K1t): ",
#if KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K1t) + (K2p[2] <- K1t): ",
                           "PT4: (K2a[4] <- K1t) + (K2p[4] <- K1t): ",
#endif
                           "PT4: (K2a[1] <- K2a) + (K2p[1] <- K2p): ",
#if KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K2a) + (K2p[2] <- K2p): ",
                           "PT4: (K2a[4] <- K2a) + (K2p[4] <- K2p): ",
#endif
                           "PT4: (K2a[1] <- K2p) + (K2p[1] <- K2a): ",
#if KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K2p) + (K2p[2] <- K2a): ",
                           "PT4: (K2a[4] <- K2p) + (K2p[4] <- K2a): ",
#endif
                           "PT4: (K2a[1] <- K2t) + (K2p[1] <- K2t): ",
#if KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K2t) + (K2p[2] <- K2t): ",
                           "PT4: (K2a[4] <- K2t) + (K2p[4] <- K2t): ",
#endif
                           "PT4: (K2t[1] <- K2a) + (K2t[1] <- K2t): ",
#if KELDYSH_FORMALISM
                           "PT4: (K2t[2] <- K2a) + (K2t[2] <- K2t): ",
                           "PT4: (K2t[4] <- K2a) + (K2t[4] <- K2t): ",
#endif
                           "PT4: K3a + K3p: ",
                           "PT4: K3t (aa) + K3t (ap) + K3t (pa): "
                           };

    // K1 in PT2
    state_datatype PT2_K1a_0 = PT2_K1a.vertex.avertex().valsmooth<k1>(input_a, PT2_K1a.vertex.tvertex());
    state_datatype PT2_K1p_0 = PT2_K1p.vertex.pvertex().valsmooth<k1>(input_p, PT2_K1p.vertex.pvertex());
    state_datatype PT2_K1t_0 = PT2_K1t.vertex.tvertex().valsmooth<k1>(input_t, PT2_K1t.vertex.avertex());

    // K1 in PT3
    state_datatype PT3_K1a_0 = PT3_K1a.vertex.avertex().valsmooth<k1>(input_a, PT3_K1a.vertex.tvertex());
    state_datatype PT3_K1p_0 = PT3_K1p.vertex.pvertex().valsmooth<k1>(input_p, PT3_K1p.vertex.pvertex());
    state_datatype PT3_K1t_0 = PT3_K1t.vertex.tvertex().valsmooth<k1>(input_t, PT3_K1t.vertex.avertex());
#if KELDYSH_FORMALISM
    state_datatype PT2_K1_exact = -(1./2.) * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 1);
    state_datatype PT3_K1_exact = -(1./2.) * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
#else
#ifdef ZERO_TEMP
    state_datatype PT2_K1_exact = - glb_U * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 1);
    state_datatype PT3_K1_exact = - glb_U * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
#else
    state_datatype PT2_K1_exact = - glb_U * pow(glb_U * boost::math::polygamma(1, (glb_Gamma + Lambda) / 2. / (glb_T * 2.*M_PI) + 0.5) / (glb_T * 2 * M_PI * M_PI), 1);
    state_datatype PT3_K1_exact = - glb_U * pow(glb_U * boost::math::polygamma(1, (glb_Gamma + Lambda) / 2. / (glb_T * 2.*M_PI) + 0.5) / (glb_T * 2 * M_PI * M_PI), 2);
#endif
#endif
    std::cout << "PT2 K1 exact: " << PT2_K1_exact << "\n";
    std::cout << "Computed value: " << PT2_K1a_0 << "\n";

    // K1 in PT4
    PT4_31_a_a1.vertex.avertex().valsmooth<k1>(input_a, PT4_31_a_a1.vertex.tvertex()); // state_datatype PT4_K1a_0_ladder
    PT4_31_p_p1.vertex.pvertex().valsmooth<k1>(input_p, PT4_31_p_p1.vertex.pvertex()); // state_datatype PT4_K1p_0_ladder
    PT4_13_a_a2.vertex.avertex().valsmooth<k1>(input_a, PT4_13_a_a2.vertex.tvertex()); // state_datatype PT4_K1a_0_nonladder
    PT4_13_p_p2.vertex.pvertex().valsmooth<k1>(input_p, PT4_13_p_p2.vertex.pvertex()); // state_datatype PT4_K1p_0_nonladder
    PT4_13_t_a2.vertex.tvertex().valsmooth<k1>(input_t, PT4_13_t_a2.vertex.avertex()); // state_datatype PT4_K1t_0_nonladder_a
    PT4_13_t_t2.vertex.tvertex().valsmooth<k1>(input_t, PT4_13_t_t2.vertex.avertex()); // state_datatype PT4_K1t_0_nonladder_t

    // K2 in PT3
    vec<state_datatype> PT3_K2a_0 (3);
    vec<state_datatype> PT3_K2p_0 (3);
    vec<state_datatype> PT3_K2t_0 (3);
    // K2 in PT4
    vec<state_datatype> PT4_K2a_0_p1 (3);
    vec<state_datatype> PT4_K2p_0_a1 (3);
    vec<state_datatype> PT4_K2a_0_t1 (3);
    vec<state_datatype> PT4_K2p_0_t1 (3);
    vec<state_datatype> PT4_K2a_0_a2 (3);
    vec<state_datatype> PT4_K2a_0_p2 (3);
    vec<state_datatype> PT4_K2a_0_t2 (3);
    vec<state_datatype> PT4_K2p_0_a2 (3);
    vec<state_datatype> PT4_K2p_0_p2 (3);
    vec<state_datatype> PT4_K2p_0_t2 (3);
    vec<state_datatype> PT4_K2t_0_a2 (3);
    vec<state_datatype> PT4_K2t_0_t2 (3);

#if MAX_DIAG_CLASS >= 2
        for (int iK2 = 0; iK2 < 3; ++iK2) {
            if (!KELDYSH && (iK2 > 0)) break;
            input_a.iK = iK2s[iK2];
            input_p.iK = iK2s[iK2];
            input_t.iK = iK2s[iK2];
            // K2 in PT3
            PT3_K2a_0[iK2] = PT3_K2a.vertex.avertex().valsmooth<k2>(input_a, PT3_K2a.vertex.tvertex());
            PT3_K2p_0[iK2] = PT3_K2p.vertex.pvertex().valsmooth<k2>(input_p, PT3_K2p.vertex.pvertex());
            PT3_K2t_0[iK2] = PT3_K2t.vertex.tvertex().valsmooth<k2>(input_t, PT3_K2t.vertex.avertex());
            // K2 in PT4
            PT4_K2a_0_p1[iK2] = PT4_31_a_p1.vertex.avertex().valsmooth<k2>(input_a, PT4_31_a_p1.vertex.tvertex());
            PT4_K2p_0_a1[iK2] = PT4_31_p_a1.vertex.pvertex().valsmooth<k2>(input_p, PT4_31_p_a1.vertex.pvertex());
            PT4_K2a_0_t1[iK2] = PT4_31_a_t1.vertex.avertex().valsmooth<k2>(input_a, PT4_31_a_t1.vertex.tvertex());
            PT4_K2p_0_t1[iK2] = PT4_31_p_t1.vertex.pvertex().valsmooth<k2>(input_p, PT4_31_p_t1.vertex.pvertex());
            PT4_K2a_0_a2[iK2] = PT4_31_a_a2.vertex.avertex().valsmooth<k2>(input_a, PT4_31_a_a2.vertex.tvertex());
            PT4_K2a_0_p2[iK2] = PT4_31_a_p2.vertex.avertex().valsmooth<k2>(input_a, PT4_31_a_p2.vertex.tvertex());
            PT4_K2a_0_t2[iK2] = PT4_31_a_t2.vertex.avertex().valsmooth<k2>(input_a, PT4_31_a_t2.vertex.tvertex());
            PT4_K2p_0_a2[iK2] = PT4_31_p_a2.vertex.pvertex().valsmooth<k2>(input_p, PT4_31_p_a2.vertex.pvertex());
            PT4_K2p_0_p2[iK2] = PT4_31_p_p2.vertex.pvertex().valsmooth<k2>(input_p, PT4_31_p_p2.vertex.pvertex());
            PT4_K2p_0_t2[iK2] = PT4_31_p_t2.vertex.pvertex().valsmooth<k2>(input_p, PT4_31_p_t2.vertex.pvertex());
            PT4_K2t_0_a2[iK2] = PT4_31_t_a2.vertex.tvertex().valsmooth<k2>(input_t, PT4_31_t_a2.vertex.avertex());
            PT4_K2t_0_t2[iK2] = PT4_31_t_t2.vertex.tvertex().valsmooth<k2>(input_t, PT4_31_t_t2.vertex.avertex());
        }

#if KELDYSH_FORMALISM
        state_datatype PT3_K2_exact =
                -(1. / 2.) * (2. - M_PI * M_PI / 4.) * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
#else
        state_datatype PT3_K2_exact = - (2. - M_PI*M_PI/4.) * glb_U * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
        std::cout << "PT3 K2 exact: " << PT3_K2_exact << "\n";
        std::cout << "Computed value: " << PT3_K2t_0[0] << "\n";
#endif
#endif

    // K3 in PT4
    if (KELDYSH) {
        input_a.iK = 5;
        input_p.iK = 5;
        input_t.iK = 5;
    }
    state_datatype PT4_K3a_0;
    state_datatype PT4_K3p_0;
    state_datatype PT4_K3t_0_aa;
    state_datatype PT4_K3t_0_ap;
    state_datatype PT4_K3t_0_pa;

    if constexpr(MAX_DIAG_CLASS == 3) {
        PT4_K3a_0 = PT4_22_a_pp.vertex.avertex().valsmooth<k3>(input_a, PT4_22_a_pp.vertex.tvertex());
        PT4_K3p_0 = PT4_22_p_aa.vertex.pvertex().valsmooth<k3>(input_p, PT4_22_p_aa.vertex.pvertex());
        PT4_K3t_0_aa = PT4_22_t_aa.vertex.tvertex().valsmooth<k3>(input_t, PT4_22_t_aa.vertex.avertex());
        PT4_K3t_0_ap = PT4_22_t_ap.vertex.tvertex().valsmooth<k3>(input_t, PT4_22_t_ap.vertex.avertex());
        PT4_K3t_0_pa = PT4_22_t_pa.vertex.tvertex().valsmooth<k3>(input_t, PT4_22_t_pa.vertex.avertex());
    }

    // values to be printed to log
    vec<state_datatype> check_values {PT2_K1a_0 + PT2_K1p_0
                            , PT2_K1t_0
                            , (PT2_K1a_0 - PT2_K1_exact)/PT2_K1_exact
                            , (-PT2_K1p_0 - PT2_K1_exact)/PT2_K1_exact
                            , (PT3_K1a_0 - PT3_K1_exact)/PT3_K1_exact
                            , (PT3_K1p_0 - PT3_K1_exact)/PT3_K1_exact
                            , (PT3_K1t_0 - PT3_K1_exact)/PT3_K1_exact
#if MAX_DIAG_CLASS > 1
                            , (PT3_K2a_0[0] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2p_0[0] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2t_0[0] - PT3_K2_exact)/PT3_K2_exact
#if KELDYSH_FORMALISM
                            , (PT3_K2a_0[1] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2p_0[1] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2t_0[1] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2a_0[2] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2p_0[2] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2t_0[2] - PT3_K2_exact)/PT3_K2_exact
#endif
                            , (PT4_K2a_0_p1[0] + PT4_K2p_0_a1[0])/(std::abs(PT4_K2a_0_p1[0]) + std::abs(PT4_K2p_0_a1[0]))
#if KELDYSH_FORMALISM
                            , (PT4_K2a_0_p1[1] + PT4_K2p_0_a1[1])/(std::abs(PT4_K2a_0_p1[1]) + std::abs(PT4_K2p_0_a1[1]))
                            , (PT4_K2a_0_p1[2] + PT4_K2p_0_a1[2])/(std::abs(PT4_K2a_0_p1[2]) + std::abs(PT4_K2p_0_a1[2]))
#endif
                            , (PT4_K2a_0_t1[0] + PT4_K2p_0_t1[0])/(std::abs(PT4_K2a_0_t1[0]) + std::abs(PT4_K2p_0_t1[0]))
#if KELDYSH_FORMALISM
                            , (PT4_K2a_0_t1[1] + PT4_K2p_0_t1[1])/(std::abs(PT4_K2a_0_t1[1]) + std::abs(PT4_K2p_0_t1[1]))
                            , (PT4_K2a_0_t1[2] + PT4_K2p_0_t1[2])/(std::abs(PT4_K2a_0_t1[2]) + std::abs(PT4_K2p_0_t1[2]))
#endif
                            , (PT4_K2a_0_a2[0] + PT4_K2p_0_p2[0])/(std::abs(PT4_K2a_0_a2[0]) + std::abs(PT4_K2p_0_p2[0]))
#if KELDYSH_FORMALISM
                            , (PT4_K2a_0_a2[1] + PT4_K2p_0_p2[1])/(std::abs(PT4_K2a_0_a2[1]) + std::abs(PT4_K2p_0_p2[1]))
                            , (PT4_K2a_0_a2[2] + PT4_K2p_0_p2[2])/(std::abs(PT4_K2a_0_a2[2]) + std::abs(PT4_K2p_0_p2[2]))
#endif
                            , (PT4_K2a_0_p2[0] + PT4_K2p_0_a2[0])/(std::abs(PT4_K2a_0_p2[0]) + std::abs(PT4_K2p_0_a2[0]))
#if KELDYSH_FORMALISM
                            , (PT4_K2a_0_p2[1] + PT4_K2p_0_a2[1])/(std::abs(PT4_K2a_0_p2[1]) + std::abs(PT4_K2p_0_a2[1]))
                            , (PT4_K2a_0_p2[2] + PT4_K2p_0_a2[2])/(std::abs(PT4_K2a_0_p2[2]) + std::abs(PT4_K2p_0_a2[2]))
#endif
                            , (PT4_K2a_0_t2[0] + PT4_K2p_0_t2[0])/(std::abs(PT4_K2a_0_t2[0]) + std::abs(PT4_K2p_0_t2[0]))
#if KELDYSH_FORMALISM
                            , (PT4_K2a_0_t2[1] + PT4_K2p_0_t2[1])/(std::abs(PT4_K2a_0_t2[1]) + std::abs(PT4_K2p_0_t2[1]))
                            , (PT4_K2a_0_t2[2] + PT4_K2p_0_t2[2])/(std::abs(PT4_K2a_0_t2[2]) + std::abs(PT4_K2p_0_t2[2]))
#endif
                            , (PT4_K2t_0_a2[0] + PT4_K2t_0_t2[0])/(std::abs(PT4_K2t_0_a2[0]) + std::abs(PT4_K2t_0_t2[0]))
#if KELDYSH_FORMALISM
                            , (PT4_K2t_0_a2[1] + PT4_K2t_0_t2[1])/(std::abs(PT4_K2t_0_a2[1]) + std::abs(PT4_K2t_0_t2[1]))
                            , (PT4_K2t_0_a2[2] + PT4_K2t_0_t2[2])/(std::abs(PT4_K2t_0_a2[2]) + std::abs(PT4_K2t_0_t2[2]))
#endif
                            , (PT4_K3a_0 + PT4_K3p_0)/(std::abs(PT4_K3a_0) + std::abs(PT4_K3p_0))
                            , (PT4_K3t_0_aa + PT4_K3t_0_ap + PT4_K3t_0_pa)/(std::abs(PT4_K3t_0_aa) + std::abs(PT4_K3t_0_ap) + std::abs(PT4_K3t_0_pa))
#endif
                            };

    // print to log
    for (unsigned int i=0; i<check_values.size(); ++i) {
        utils::print(check_labels[i], check_values[i], true);
    }

    utils::print("----------------------", true);

    /*
    // Compute K1a contributions in PT4, using
    // (22):   K1a in PT2
    // (13_1): K1a in PT3
    // (13_2): K2a in PT3
    // (31_2): K2a in PT3
    State<Q> PT4_K1a22 (Lambda);
    State<Q> PT4_K1a13_1 (Lambda);
    State<Q> PT4_K1a13_2 (Lambda);
    State<Q> PT4_K1a31_2 (Lambda);

    t0 = utils::get_time();
    bubble_function(PT4_K1a22.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_K1a13_1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_K1a13_2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_K1a31_2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    utils::print("Computed K1 in PT4.", true);
    utils::get_time(t0);

    // FDT checks
    utils::print("Check K1a in PT4 (22):", true);
    check_FDTs(PT4_K1a22);
    utils::print("Check K1a in PT4 (13_1):", true);
    check_FDTs(PT4_K1a13_1);
    utils::print("Check K1a in PT4 (13_2):", true);
    check_FDTs(PT4_K1a13_2);
    utils::print("Check K1a in PT4 (31_2):", true);
    check_FDTs(PT4_K1a31_2);
    // */
}

/**
 * Test K3 dynamics by computing SE diagrams in PT4 using different PT4 vertices, which should all give the same result.
 */
template <typename Q>
void test_K3_dynamics_SE_PT4(double Lambda) {
    // Initialize a bare state
    State<Q> bare (Lambda);
    bare.initialize();

    // Initialize a bare propagator
    Propagator<Q> G(Lambda, bare.selfenergy, 'g');

    // Compute K1a and K1p in PT2
    State<Q> PT2_K1a (Lambda);
    State<Q> PT2_K1p (Lambda);

    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false, stdConfig);

    // Compute K1a and K1p in PT3
    State<Q> PT3_K1a (Lambda);
    State<Q> PT3_K1p (Lambda);

    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT3_K1p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'p', false, stdConfig);

    // Compute K3a, K1p (ladder), K3p, K1a (ladder) in PT4
    State<Q> PT4_22_a_pp (Lambda);
    State<Q> PT4_13_p_p1 (Lambda);
    State<Q> PT4_22_p_aa (Lambda);
    State<Q> PT4_13_a_a1 (Lambda);

    bubble_function(PT4_22_a_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_13_p_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_22_p_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT4_13_a_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false, stdConfig);

    // a-channel:
    // close K3a (single diagram: PT2_K1p - a-bubble - PT2_K1p)
    State<Q> SE_K3a (Lambda);
    loop<false,0>(SE_K3a.selfenergy, PT4_22_a_pp.vertex, G);

    // close K1p ladder (use 1-3 vertex in p-channel, since it contains only p-ladder)
    State<Q> SE_K1p_ladder (Lambda);
    loop<false,0>(SE_K1p_ladder.selfenergy, PT4_13_p_p1.vertex, G);

    // p-channel:
    // close K3p (single diagram: PT2_K1a - p-bubble - PT2_K1a)
    State<Q> SE_K3p (Lambda);
    loop<false,0>(SE_K3p.selfenergy, PT4_22_p_aa.vertex, G);

    // close K1a ladder (use 1-3 vertex in a-channel, since it contains only a-ladder)
    State<Q> SE_K1a_ladder (Lambda);
    loop<false,0>(SE_K1a_ladder.selfenergy, PT4_13_a_a1.vertex, G);

    write_state_to_hdf("SE_K3a_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K3a);
    write_state_to_hdf("SE_K1p_ladder_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K1p_ladder);
    write_state_to_hdf("SE_K3p_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K3p);
    write_state_to_hdf("SE_K1a_ladder_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K1a_ladder);
}


/**
 * Function to test correctness of K2a when calculating a susceptibility (a K1a-object) //Notice that the same calculation
 * can be performed in the p-channel.
 * @param Lambda : Scale at which the calculation is done.
 */
template <typename Q>
void test_K2_correctness(double Lambda){

    bool write_flag = true; //Write out results in a HDF5 file

    State<Q> bare (Lambda);   //Bare state
    bare.initialize();  //Initialize bare state

    Propagator<Q> G(Lambda, bare.selfenergy, 'g'); //Bare propagator
    Propagator<Q> S(Lambda, bare.selfenergy, 's'); //Bare single-scale propagator

    //Create states for K1-calculations
    State<Q> PT2_K1a (Lambda);
    State<Q> PT2_K1p (Lambda);
    State<Q> PT2_K1t (Lambda);

    //Save K1-bubbles in separate objects - SOPT
    double t0 = utils::get_time();
    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false, stdConfig);
    bubble_function(PT2_K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false, stdConfig);
    utils::get_time(t0);

    State<Q> PT2_SE_a (Lambda);
    State<Q> PT2_SE_p (Lambda);
    State<Q> PT2_SE_t (Lambda);
    State<Q> PT2_SE_p_1 (Lambda);
    State<Q> PT2_SE_p_4 (Lambda);
    State<Q> PT2_SE_p_5 (Lambda);

    loop<true,0>(PT2_SE_a.selfenergy, PT2_K1a.vertex, S);
    loop<true,0>(PT2_SE_p.selfenergy, PT2_K1p.vertex, S);
    loop<true,0>(PT2_SE_t.selfenergy, PT2_K1t.vertex, S);

    State<Q> PT3_K2a (Lambda);    //Create state for K2a calculation
    State<Q> PT3_K2a_ia (Lambda);
    State<Q> PT3_K2a_ib (Lambda);
    State<Q> PT3_K2a_iia (Lambda);
    State<Q> PT3_K2a_iib (Lambda);
    State<Q> PT3_K2a_iva (Lambda);
    State<Q> PT3_K2a_ivb (Lambda);
    State<Q> PT3_K2a_t (Lambda);

    State<Q> PT3_K2p (Lambda);
    State<Q> PT3_K2t (Lambda);
    State<Q> PT3_K2t_a (Lambda);
    State<Q> PT3_K2t_p (Lambda);

    //Do appropriate calculation for K2a with K1p and K1t being fed back into the left vertex. Notice part = 'L' to ensure
    //that the correct contributions are added on both sides. - TOPT
    t0 = utils::get_time();
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, stdConfig);   // K2a in PT3

#ifdef DEBUG_MODE
    bubble_function(PT3_K2a_ia.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 9, 6);
    bubble_function(PT3_K2a_ib.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 6, 9);
    bubble_function(PT3_K2a_iia.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 11, 14);
    bubble_function(PT3_K2a_iib.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 7, 9);
    bubble_function(PT3_K2a_iva.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 9, 6);
    bubble_function(PT3_K2a_ivb.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 6, 15);
#endif

    bubble_function(PT3_K2a_t.vertex, PT2_K1t.vertex, bare.vertex, G, G, 'a', false, stdConfig);   // K2a in PT3
    //PT3_K2a = read_hdf("PT4_check_of_K2a_K2_switchedcc_adap_m3m9_g501_101_nI1501_state_PT3_K2a", 0, 1);

    bubble_function(PT3_K2p.vertex, PT2_K1a.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'p', false, stdConfig);    // K2p  in PT3
    //bubble_function(PT3_K2t.vertex, PT2_K1a.vertex + PT2_K1p.vertex, bare.vertex, G, G, 't', false, stdConfig);    // K2t  in PT3
    bubble_function(PT3_K2t_a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 't', false, stdConfig);    // K2t  in PT3
    bubble_function(PT3_K2t_p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 't', false, stdConfig);    // K2t  in PT3
    PT3_K2t.vertex = PT3_K2t_a.vertex + PT3_K2t_p.vertex;
    utils::get_time(t0);

    // full K2 in PT3
    State<Q> PT3_K2 (Lambda);
    PT3_K2.vertex.avertex() = PT3_K2a_t.vertex.avertex();
    //PT3_K2.vertex.pvertex() = PT3_K2p.vertex.pvertex();
    PT3_K2.vertex.tvertex() = PT3_K2t_a.vertex.tvertex();

    State<Q> PT3_K2at (Lambda);

    // K2 contribution to self-energy flow
    State<Q> PT3_SE (Lambda);
    State<Q> PT3_SE_a (Lambda);
    State<Q> PT3_SE_p (Lambda);
    State<Q> PT3_SE_t (Lambda);
    State<Q> PT3_SE_t_1 (Lambda);
    State<Q> PT3_SE_t_4 (Lambda);
    State<Q> PT3_SE_t_5 (Lambda);
    State<Q> PT3_SE_t_a (Lambda);
    State<Q> PT3_SE_t_a_1 (Lambda);
    State<Q> PT3_SE_t_a_4 (Lambda);
    State<Q> PT3_SE_t_a_5 (Lambda);
    State<Q> PT3_SE_t_p (Lambda);
    State<Q> PT3_SE_t_p_1 (Lambda);
    State<Q> PT3_SE_t_p_4 (Lambda);
    State<Q> PT3_SE_t_p_5 (Lambda);

    loop<true,0>(PT3_SE.selfenergy, PT3_K2.vertex, S);
    loop<true,0>(PT3_SE_a.selfenergy, PT3_K2a.vertex, S);
    loop<true,0>(PT3_SE_p.selfenergy, PT3_K2p.vertex, S);
    loop<true,0>(PT3_SE_t.selfenergy, PT3_K2t.vertex, S);

    loop<true,0>(PT3_SE_t_a.selfenergy, PT3_K2t_a.vertex, S);
    loop<true,0>(PT3_SE_t_p.selfenergy, PT3_K2t_p.vertex, S);


    State<Q> PT3_K1a (Lambda);    //Create state to compare with K1a
    t0 = utils::get_time();
    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    utils::get_time(t0);

    State<Q> PT123_a = bare + PT2_K1a + PT3_K1a + PT3_K2a;  //Create vertex of the right side of BSE

    State<Q> PT4_K1a22 (Lambda);
    State<Q> PT4_K1a13_1 (Lambda);
    State<Q> PT4_K1a13_2 (Lambda);
    State<Q> PT4_K1a31_2 (Lambda);
    State<Q> PT4_K1a13_2_11e (Lambda); // A
    State<Q> PT4_K1a13_2_21e (Lambda); // B
    State<Q> PT4_K1a13_2_11o (Lambda); // C
    State<Q> PT4_K1a13_2_21o (Lambda); // D
    State<Q> PT4_K1a13_2_12o (Lambda); // TST3TC D
    State<Q> PT4_K1a13_2_22o (Lambda); // F
    State<Q> PT4_K1a13_2_ia (Lambda);
    State<Q> PT4_K1a13_2_ib (Lambda);
    State<Q> PT4_K1a13_2_iia (Lambda);
    State<Q> PT4_K1a13_2_iib (Lambda);
    State<Q> PT4_K1a13_2_iva (Lambda);
    State<Q> PT4_K1a13_2_ivb (Lambda);

    t0 = utils::get_time();
    bubble_function(PT4_K1a22.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_K1a13_1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false, stdConfig);

    bubble_function(PT4_K1a13_2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, stdConfig);
    bubble_function(PT4_K1a31_2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false, stdConfig);
    utils::get_time(t0);
    t0 = utils::get_time();
#ifdef DEBUG_MODE
    bubble_function(PT4_K1a13_2_11e.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 0, 16, 16, 16); // A
    bubble_function(PT4_K1a13_2_21e.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 1, 16, 16, 16); // B
    bubble_function(PT4_K1a13_2_11o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 2, 16, 16, 16); // C
    bubble_function(PT4_K1a13_2_21o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 3, 16, 16, 16); // D
    bubble_function(PT4_K1a13_2_12o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 10, 16, 16, 16); // TST3TC D
    bubble_function(PT4_K1a13_2_22o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 11, 16, 16, 16); // F

    bubble_function(PT4_K1a13_2_ia.vertex, bare.vertex, PT3_K2a_ia.vertex, G, G, 'a', false, 11, 6, 16, 16);
    bubble_function(PT4_K1a13_2_ib.vertex, bare.vertex, PT3_K2a_ib.vertex, G, G, 'a', false, 2, 9, 16, 16);
    bubble_function(PT4_K1a13_2_iia.vertex, bare.vertex, PT3_K2a_iia.vertex, G, G, 'a', false, 11, 6, 16, 16);
    bubble_function(PT4_K1a13_2_iib.vertex, bare.vertex, PT3_K2a_iib.vertex, G, G, 'a', false, 10, 11, 16, 16);
    bubble_function(PT4_K1a13_2_iva.vertex, bare.vertex, PT3_K2a_iva.vertex, G, G, 'a', false, 11, 15, 16, 16);
    bubble_function(PT4_K1a13_2_ivb.vertex, bare.vertex, PT3_K2a_ivb.vertex, G, G, 'a', false, 2, 9, 16, 16);
#endif
    utils::get_time(t0);

    cvec K1a_diff(nBOS);
    for(int iw=0; iw<nBOS; ++iw){
        K1a_diff[iw] = PT4_K1a22.vertex.avertex().K1.val(0, iw,0, 0) - PT2_K1a.vertex.avertex().K1.val(0, iw,0, 0);
    }

    utils::print("Testing correctness of K2a. Using U=" +std::to_string(glb_U)+ " and Lambda="+std::to_string(Lambda)
        +", the maximal difference between direct K1a and K1a over integration of K2a is " +std::to_string(K1a_diff.max_norm())+"." , true);
    if(write_flag) write_h5_rvecs( data_dir + "/PT4_check_of_K2a_cleanup_GL_gW20_51_21_nI1501_U1", {"w",
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
                                                       "PT3_K2t_p_a_R", "PT3_K2t_p_a_I",
                                                       "PT3_K2t_p_p_R", "PT3_K2t_p_p_I",
                                                       "PT3_K2t_p_t_R", "PT3_K2t_p_t_I",
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
                                   PT2_K1a.vertex.avertex().K1.real(), PT2_K1a.vertex.avertex().K1.imag(),
                                   PT2_K1p.vertex.pvertex().K1.real(), PT2_K1p.vertex.pvertex().K1.imag(),
                                   PT2_K1t.vertex.tvertex().K1.real(), PT2_K1t.vertex.tvertex().K1.imag(),
                                   PT2_SE_a.selfenergy.Sigma.real(), PT2_SE_a.selfenergy.Sigma.imag(),
                                   PT2_SE_p.selfenergy.Sigma.real(), PT2_SE_p.selfenergy.Sigma.imag(),
                                   PT2_SE_t.selfenergy.Sigma.real(), PT2_SE_t.selfenergy.Sigma.imag(),
                                   PT2_SE_p_1.selfenergy.Sigma.real(), PT2_SE_p_1.selfenergy.Sigma.imag(),
                                   PT2_SE_p_4.selfenergy.Sigma.real(), PT2_SE_p_4.selfenergy.Sigma.imag(),
                                   PT2_SE_p_5.selfenergy.Sigma.real(), PT2_SE_p_5.selfenergy.Sigma.imag(),
                                   PT3_K1a.vertex.avertex().K1.real(), PT3_K1a.vertex.avertex().K1.imag(),
                                   PT3_K2a.vertex.avertex().K2.real(), PT3_K2a.vertex.avertex().K2.imag(),
                                   PT3_K2a_ia.vertex.avertex().K2.real(), PT3_K2a_ia.vertex.avertex().K2.imag(),
                                   PT3_K2a_ib.vertex.avertex().K2.real(), PT3_K2a_ib.vertex.avertex().K2.imag(),
                                   PT3_K2a_iia.vertex.avertex().K2.real(), PT3_K2a_iia.vertex.avertex().K2.imag(),
                                   PT3_K2a_iib.vertex.avertex().K2.real(), PT3_K2a_iib.vertex.avertex().K2.imag(),
                                   PT3_K2a_iva.vertex.avertex().K2.real(), PT3_K2a_iva.vertex.avertex().K2.imag(),
                                   PT3_K2a_ivb.vertex.avertex().K2.real(), PT3_K2a_ivb.vertex.avertex().K2.imag(),
                                   PT3_K2a_t.vertex.avertex().K2.real(), PT3_K2a_t.vertex.avertex().K2.imag(),
                                   PT3_K2p.vertex.avertex().K2.real(), PT3_K2p.vertex.avertex().K2.imag(),
                                   PT3_K2t.vertex.avertex().K2.real(), PT3_K2t.vertex.avertex().K2.imag(),
                                   PT3_K2t_a.vertex.avertex().K2.real(), PT3_K2t_a.vertex.avertex().K2.imag(),
                                   PT3_K2t_a.vertex.pvertex().K2.real(), PT3_K2t_a.vertex.pvertex().K2.imag(),
                                   PT3_K2t_a.vertex.tvertex().K2.real(), PT3_K2t_a.vertex.tvertex().K2.imag(),
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
                                   PT4_K1a22.vertex.avertex().K1.real(), PT4_K1a22.vertex.avertex().K1.imag(),
                                   PT4_K1a13_1.vertex.avertex().K1.real(), PT4_K1a13_1.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2.vertex.avertex().K1.real(), PT4_K1a13_2.vertex.avertex().K1.imag(),
                                   PT4_K1a31_2.vertex.avertex().K1.real(), PT4_K1a31_2.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_11e.vertex.avertex().K1.real(), PT4_K1a13_2_11e.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_21e.vertex.avertex().K1.real(), PT4_K1a13_2_21e.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_11o.vertex.avertex().K1.real(), PT4_K1a13_2_11o.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_21o.vertex.avertex().K1.real(), PT4_K1a13_2_21o.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_12o.vertex.avertex().K1.real(), PT4_K1a13_2_12o.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_22o.vertex.avertex().K1.real(), PT4_K1a13_2_22o.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_ia.vertex.avertex().K1.real(), PT4_K1a13_2_ia.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_ib.vertex.avertex().K1.real(), PT4_K1a13_2_ib.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_iia.vertex.avertex().K1.real(), PT4_K1a13_2_iia.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_iib.vertex.avertex().K1.real(), PT4_K1a13_2_iib.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_iva.vertex.avertex().K1.real(), PT4_K1a13_2_iva.vertex.avertex().K1.imag(),
                                   PT4_K1a13_2_ivb.vertex.avertex().K1.real(), PT4_K1a13_2_ivb.vertex.avertex().K1.imag()});

    //write_state_to_hdf("PT4_check_of_K2a_K2_switchedcc_t_update_symmrev_new11_SE_symm_full_adap_m3m9_gW10_501_101_nI1501_U1_state_PT3_K2a", 0, 1, PT3_K2a);
}

/**
 * Master function to test both consistency and correctness of K2-class
 * @param Lambda
 */
template<typename Q>
void test_K2(double Lambda, bool test_consistency){


    //First test consistency
    if(test_consistency) {
        bool K2a = test_K2_consistency<Q>(Lambda, 'a');    //Consistency of a-channel
        bool K2p = test_K2_consistency<Q>(Lambda, 'p');    //Consistency of p-channel
        bool K2t = test_K2_consistency<Q>(Lambda, 't');    //Consistency of t-channel

        if(K2a&&K2p&&K2t)
            test_K2_correctness<Q>(Lambda);
    }

    test_K2_correctness<Q>(Lambda);

}

/**
 * test-integrand for below function test_integrate_over_K1()
 */
template <typename Q>
class TestIntegrandK1a{
    public:
        const double Lambda;
        const double w, v, vp;
        const char channel;
        State<Q> SOPTstate;
        TestIntegrandK1a(const double Lambda_in, const char channel_in, const double w, const double v, const double vp)
        : Lambda(Lambda_in), channel(channel_in), w(w), v(v), vp(vp), SOPTstate(Lambda) {
            //SOPTstate = State<Q>(Lambda);
            SOPTstate.initialize();             // initialize state

            sopt_state(SOPTstate, Lambda);
        }

    void save_integrand(double vmax) {
        int npoints = 1e5;
        FrequencyGrid bfreqs('b', 1, Lambda);
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        double spacing = (bfreqs.t_upper - bfreqs.t_lower) / (double)npoints;
        for (int i=0; i<npoints; ++i) {
            double vpp = bfreqs.frequency_from_t(bfreqs.t_lower + i * spacing);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i + 1) * (1 - frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    void save_state() {
            write_state_to_hdf("SOPT_state.h5", 1.8, 1, SOPTstate);
        }

    auto operator() (double vpp) const -> Q {
            VertexInput input(0, vpp, v, vp + vpp, 0, 0, channel);
            return SOPTstate.vertex.half1().avertex.value(input, this->SOPTstate.vertex.half1().avertex) ;
            //return vpp*vpp;
        }
    };


/**
 * test function; can be used to test the integration and the interpolation
 * @tparam Q
 * @param Lambda
 */
template <typename Q>
void test_integrate_over_K1(double Lambda) {
    double v = 1e2;
    double vp = -1e2;
    TestIntegrandK1a<Q> IntegrandK1a(Lambda, 'a', 0., v, vp);
    IntegrandK1a.save_integrand(1e4);
    double exact_result = -glb_U*glb_U/4;

    double vmax = 1e2;
    //Integrand& integrand, const double vmin, const double vmax, double w_half, const vec<double>& freqs, const double Delta, const int num_freqs;
    vec<double> freqs = {};
    int num_freqs = 0;
    double Delta = (glb_Gamma+Lambda)/2.;
    Q result = integrator_Matsubara_T0(IntegrandK1a, -vmax, vmax, 0, freqs, Delta)/ (2*M_PI);
    utils::print("Result of the integration over K1a:", result, true);
    utils::print("relative difference: ", (result-exact_result)/exact_result, true);


    TestIntegrandK1a<Q> IntegrandK1a2(Lambda, 't', (v-vp)*2, v, vp);
    Q result2 = integrator_Matsubara_T0(IntegrandK1a2, -vmax, vmax, (v-vp), freqs, Delta)/ (2*M_PI);
    utils::print("Result of the integration over K1a from the t-channel:", result2, true);
    utils::print("relative difference: ", (result2-exact_result)/exact_result, true);


    TestIntegrandK1a<Q> IntegrandK1a3(Lambda, 'p', (v+vp)*2, v, vp);
    Q result3 = integrator_Matsubara_T0(IntegrandK1a3, -vmax, vmax, (v+vp), freqs, Delta)/ (2*M_PI);
    utils::print("Result of the integration over K1a from p-channel:", result3, true);
    utils::print("relative difference: ", (result3-exact_result)/exact_result, true);


    IntegrandK1a2.save_integrand(vmax);
    IntegrandK1a3.save_integrand(vmax);
    IntegrandK1a.save_state();
}

#if not KELDYSH_FORMALISM and defined(ZERO_TEMP)
auto SOPT_K1a(double w, double Lambda) -> double {
    double Delta = (glb_Gamma + Lambda) / 2.;
    double result;
    if (w == 0.)
        result = - glb_U*glb_U * Delta / M_PI / (Delta*Delta + glb_mu*glb_mu);
    else
        result = - glb_U*glb_U * Delta/M_PI / (std::abs(w)*(2*Delta + std::abs(w))) * log(1 + (std::abs(w)*(2*Delta + std::abs(w)))/(Delta*Delta + glb_mu*glb_mu));
    return result;
}
auto SOPT_K1a_diff(double w, double Lambda) -> double {
    double Delta = (glb_Gamma + Lambda) / 2.;
    if (abs(w) < inter_tol) return - glb_U*glb_U      /2 / M_PI / (Delta*Delta + glb_mu*glb_mu)
                        + glb_U*glb_U * Delta / M_PI / (Delta*Delta + glb_mu*glb_mu) / (Delta*Delta + glb_mu*glb_mu) * Delta ;
    double term1 =  - glb_U*glb_U      / 2/M_PI / (std::abs(w)*(2*Delta + std::abs(w))) * log(1 + (std::abs(w)*(2*Delta + std::abs(w)))/(Delta*Delta + glb_mu*glb_mu));
    double term2 =  + glb_U*glb_U * Delta /M_PI / (std::abs(w)*(2*Delta + std::abs(w))) * log(1 + (std::abs(w)*(2*Delta + std::abs(w)))/(Delta*Delta + glb_mu*glb_mu)) / (2*Delta + std::abs(w));
    double term3 =  - glb_U*glb_U * Delta /M_PI / (std::abs(w)*(2*Delta + std::abs(w)))     /(   1 + (std::abs(w)*(2*Delta + std::abs(w)))/(Delta*Delta + glb_mu*glb_mu))
            *(
                    std::abs(w)                      / (Delta*Delta + glb_mu*glb_mu)
                   -std::abs(w) * (2*Delta + std::abs(w)) / (Delta*Delta + glb_mu*glb_mu) / (Delta*Delta + glb_mu*glb_mu) * Delta
                    );
    return term1 + term2 + term3;
}

/**
 * integrand for selfenergy in TOPT (SIAM)
 */
template <typename Q>
class Integrand_TOPT_SE{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    const char channel = 'a';
    const bool diff;

    Propagator<Q>& barePropagator;
    Integrand_TOPT_SE(double Lambda_in, double w, double v, bool diff, Propagator<Q>& barePropagator)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), barePropagator(barePropagator) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a(- v + vpp, Lambda) * barePropagator.valsmooth(0, vpp, 0) ;
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class Integrand_TOPTK2a{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    Integrand_TOPTK2a(double Lambda_in, double w, double v, bool diff, Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * glb_U ;
        //return vpp*vpp;
    }
};

/**
 * integrand for K3a in FOPT (SIAM)
 */
template <typename Q>
class Integrand_FOPTK3a{
public:
    const double Lambda;
    const double w, v, vp;
    const bool diff;
    const char channel = 'a';
    Bubble<Q> Pi;
    Integrand_FOPTK3a(double Lambda_in, double w, double v, double vp, bool diff, Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) { }

    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=1; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * SOPT_K1a(vp + vpp, Lambda) ;
        //return vpp*vpp;
    }
};
#if REG == 2
template <typename Q>
void test_PT_state(std::string outputFileName, double Lambda, bool diff) {
    const int it_spin = 0;
    double vmax = 100;
    double Delta = (glb_Gamma + Lambda) / 2.;
    State<Q> bareState (Lambda);bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    Bubble<Q> Pi(barePropagator, barePropagator, diff);

    State<Q> state_cpp (Lambda);   // create final and initial state
    state_cpp.initialize();             // initialize state
#ifdef ADAPTIVE_GRID
    const int N_iterations = 1;
    write_state_to_hdf(data_dir+"PTstate_preOpt.h5", Lambda_ini,  N_iterations+1, state_cpp);  // save the initial state to hdf5 file
    write_state_to_hdf(data_dir+"PTstate_postOpt.h5", Lambda_ini,  N_iterations+1, state_cpp);  // save the initial state to hdf5 file

    for (int it = 0; it < N_iterations; it++) {
        State<Q> state_temp (state_cpp, Lambda);   // copy frequency grids //
        state_temp.initialize();             // initialize state
        fopt_state(state_temp, Lambda);


        add_state_to_hdf(data_dir+"PTstate_preOpt.h5", it+1, state_temp);  // save the initial state to hdf5 file
        //state_temp.vertex.half1().check_vertex_resolution();
        state_temp.findBestFreqGrid(true);
        state_temp.vertex.half1().check_vertex_resolution();
        state_temp.analyze_tails();

        state_cpp.set_frequency_grid(state_temp);


        add_state_to_hdf(data_dir + "PTstate_postOpt.h5", it+1, state_temp);  // save the initial state to hdf5 file
    }
#endif

    fopt_state(state_cpp, Lambda);
    //state_cpp.vertex.half1().check_vertex_resolution();
    //state_cpp.analyze_tails();

    write_state_to_hdf(outputFileName + "_cpp", Lambda, 1, state_cpp);

    State<state_datatype> PT_state(state_cpp, Lambda);

    // compute SOPT self-energy (numerically exact)
    my_defs::SE::index_type idx_SE;
    idx_SE[my_defs::SE::keldysh] = 0;
    idx_SE[my_defs::SE::internal] = 0;
    for (int i = 0; i<nFER; i++) {
        idx_SE[my_defs::SE::nu] = i;
        double v = PT_state.selfenergy.Sigma.frequencies.  primary_grid.get_frequency(i);
        Integrand_TOPT_SE<Q> IntegrandSE(Lambda, 0, v, diff, barePropagator);
        Q val_SE = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandSE, -vmax, vmax, std::abs(0.), {v}, Delta, true);
        PT_state.selfenergy.setself(0, i, 0, val_SE);
    }

    // get SOPT K1 (exact)
    my_defs::K1::index_type idx_K1;
    idx_K1[my_defs::K1::keldysh] = 0;
    idx_K1[my_defs::K1::spin] = it_spin;
    idx_K1[my_defs::K1::internal] = 0;
    for (int i = 0; i<nBOS; i++) {
        idx_K1[my_defs::K1::omega] = i;
        double w = PT_state.vertex.avertex().K1.frequencies.  primary_grid.get_frequency(i);
        Q val_K1;
        if (diff) val_K1 = SOPT_K1a_diff(w, Lambda);
        else val_K1 = SOPT_K1a(w, Lambda);
        //PT_state.vertex.avertex().K1.setvert( val_K1 - val_K1*val_K1/glb_U, i, 0, 0);
        PT_state.vertex.avertex().K1.setvert( val_K1, idx_K1);
        w = PT_state.vertex.pvertex().K1.frequencies.  primary_grid.get_frequency(i);
        if (diff) val_K1 = SOPT_K1a_diff(w, Lambda);
        else val_K1 = SOPT_K1a(w, Lambda);
        //PT_state.vertex.pvertex().K1.setvert( -val_K1 - val_K1*val_K1/glb_U, i, 0, 0);
        PT_state.vertex.pvertex().K1.setvert( -val_K1, idx_K1);
        w = PT_state.vertex.tvertex().K1.frequencies.  primary_grid.get_frequency(i);
        if (diff) val_K1 = SOPT_K1a_diff(w, Lambda);
        else val_K1 = SOPT_K1a(w, Lambda);
        PT_state.vertex.tvertex().K1.setvert( -val_K1*val_K1/glb_U, idx_K1);
    }

#if MAX_DIAG_CLASS > 1
    // compute TOPT K2 (eye diagrams) (numerically exact)
    my_defs::K2::index_type idx_K2;
    idx_K2[my_defs::K2::keldysh] = 0;
    idx_K2[my_defs::K2::spin] = it_spin;
    idx_K2[my_defs::K2::internal] = 0;
    for (int i = 0; i<nBOS2; i++) {
        for (int j = 0; j<nFER2; j++) {
            idx_K2[my_defs::K2::omega] = i;
            idx_K2[my_defs::K2::nu] = j;
            double w, v;
            PT_state.vertex.avertex().K2.frequencies.get_freqs_w(w, v, i, j);
            Integrand_TOPTK2a<Q> IntegrandK2(Lambda, w, v, diff, Pi);
            Q val_K2 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK2, -vmax, vmax, std::abs(w/2), {v, w+v, w-v}, Delta, true);
            PT_state.vertex.avertex().K2.setvert(val_K2, idx_K2);
            PT_state.vertex.pvertex().K2.setvert(val_K2, idx_K2);
            PT_state.vertex.tvertex().K2.setvert(val_K2, idx_K2);
        }
    }
#endif
#if MAX_DIAG_CLASS > 2
    // compute FOPT K3 (numerically exact)
    my_defs::K3::index_type idx_K3;
    idx_K2[my_defs::K3::keldysh] = 0;
    idx_K2[my_defs::K3::spin] = it_spin;
    idx_K2[my_defs::K3::internal] = 0;
    for (int i = 0; i<nBOS3; i++) {
        for (int j = 0; j<nFER3; j++) {
            for (int k = 0; k<(GRID!=2 ? nFER3 : (nFER3-1)/2+1); k++) {
                idx_K3[my_defs::K3::omega] = i;
                idx_K3[my_defs::K3::nu] = j;
                idx_K3[my_defs::K3::nup] = k;
                double w, v, vp;
                PT_state.vertex.avertex().K3.frequencies.get_freqs_w(w, v, vp, i, j, k);
                Integrand_FOPTK3a<Q> IntegrandK3(Lambda, w, v, vp, diff, Pi);
                Q val_K3 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK3, -vmax, vmax, std::abs(w/2), {v, vp, w+v, w-v, w+vp, w-vp}, Delta, true);
                PT_state.vertex.avertex().K3.setvert(val_K3, idx_K3);

                PT_state.vertex.pvertex().K3.frequencies.get_freqs_w(w, v, vp, i, j, k);
                Integrand_FOPTK3a<Q> IntegrandK3_2(Lambda, w, v, vp, diff, Pi);
                val_K3 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK3_2, -vmax, vmax, std::abs(w/2), {v, vp, w+v, w-v, w+vp, w-vp}, Delta, true);
                PT_state.vertex.pvertex().K3.setvert(-val_K3, idx_K3);

                PT_state.vertex.tvertex().K3.frequencies.get_freqs_w(w, v, vp, i, j, k);
                Integrand_FOPTK3a<Q> IntegrandK3_3(Lambda, w, v, vp, diff, Pi);
                val_K3 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK3_3, -vmax, vmax, std::abs(w/2), {v, vp, w+v, w-v, w+vp, w-vp}, Delta, true);
                Integrand_FOPTK3a<Q> IntegrandK3_ap(Lambda, w, -v, vp, diff, Pi);
                Q val_K3_ap = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK3_ap, -vmax, vmax, std::abs(w/2), {v, vp, w+v, w-v, w+vp, w-vp}, Delta, true);
                PT_state.vertex.tvertex().K3.setvert(-2.*(val_K3-val_K3_ap), idx_K3);
            }
        }
    }
#endif

    write_state_to_hdf(outputFileName + "_exact", Lambda, 1, PT_state);


    State<Q> state_diff = state_cpp - PT_state;

    write_state_to_hdf(outputFileName + "_diff", Lambda, 1, state_diff);
    utils::print("SE-difference: ", state_diff.selfenergy.Sigma.get_vec().max_norm() / PT_state.selfenergy.Sigma.get_vec().max_norm(), true);
    utils::print("K1a-difference: ", state_diff.vertex.avertex().K1.get_vec().max_norm() / PT_state.vertex.avertex().K1.get_vec().max_norm(), true);
    utils::print("K1p-difference: ", state_diff.vertex.pvertex().K1.get_vec().max_norm() / PT_state.vertex.pvertex().K1.get_vec().max_norm(), true);
    utils::print("K1t-difference: ", state_diff.vertex.tvertex().K1.get_vec().max_norm() / PT_state.vertex.tvertex().K1.get_vec().max_norm(), true);


    if (MAX_DIAG_CLASS > 1) {
        utils::print("K2a-difference: ", state_diff.vertex.avertex().K2.get_vec().max_norm() / PT_state.vertex.avertex().K2.get_vec().max_norm(), true);
        utils::print("K2p-difference: ", state_diff.vertex.pvertex().K2.get_vec().max_norm() / PT_state.vertex.pvertex().K2.get_vec().max_norm(), true);
        utils::print("K2t-difference: ", state_diff.vertex.tvertex().K2.get_vec().max_norm() / PT_state.vertex.tvertex().K2.get_vec().max_norm(), true);
    }
    if (MAX_DIAG_CLASS > 2) {
        utils::print("K3a-difference: ", state_diff.vertex.avertex().K3.get_vec().max_norm() / PT_state.vertex.avertex().K3.get_vec().max_norm(), true);
        utils::print("K3p-difference: ", state_diff.vertex.pvertex().K3.get_vec().max_norm() / PT_state.vertex.pvertex().K3.get_vec().max_norm(), true);
        utils::print("K3t-difference: ", state_diff.vertex.tvertex().K3.get_vec().max_norm() / PT_state.vertex.tvertex().K3.get_vec().max_norm(), true);
    }

}

#endif



/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class K1rdot_PIa_K1p_exact_K2{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    const char channel = 'a';
    const bool diff;

    const Bubble<Q> Pi;
    K1rdot_PIa_K1p_exact_K2(double Lambda_in, double w, double v, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        //if (std::abs(std::abs(vpp)-std::abs(w/2.)) < 1e-10) {
        //    return 0.;
        //}
        //else {
            return SOPT_K1a_diff(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * glb_U ;
        //}
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class K1rdot_PIa_K1p_exact_K3{
public:
    const double Lambda;
    const double w, v, vp;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    K1rdot_PIa_K1p_exact_K3(double Lambda_in, double w, double v, double vp, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        if (isfinite(vpp)) return SOPT_K1a(vp+ vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * SOPT_K1a_diff(v + vpp, Lambda) ;
        else return 0.;
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class K1p_PIa_K1rdot_exact_K3{
public:
    const double Lambda;
    const double w, v, vp;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    K1p_PIa_K1rdot_exact_K3(double Lambda_in, double w, double v, double vp, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        if (isfinite(vpp)) return SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * SOPT_K1a_diff(vp + vpp, Lambda);
        else return 0.;
        //return vpp*vpp;
    }
};


/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class IntegranddGammaC_exact_K1{
public:
    const double Lambda;
    const double w;
    const double v = 0.;
    const double vp = 0.;
    double vmax = 1e3;
    const double Delta = (Lambda * glb_Gamma)/2;
    const char channel = 'a';
    const bool diff;

    const Bubble<Q> Pi;
    IntegranddGammaC_exact_K1(double Lambda_in, double w, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), diff(diff), Pi(Pi) {}


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {

        K1rdot_PIa_K1p_exact_K2<state_datatype> IntegrandK2(Lambda, w, vpp, false, Pi);
        state_datatype val_K2 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK2, -vmax, vmax, std::abs(w/2), {vpp, vp}, Delta);
        return -glb_U * Pi.value(0, w, vpp, 0, 'a') * val_K2;
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class IntegranddGammaC_exact_K2{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    double vmax = 1e3;
    const double Delta = (Lambda * glb_Gamma)/2;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    IntegranddGammaC_exact_K2(double Lambda_in, double w, double v, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), Pi(Pi) {}


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {

        K1rdot_PIa_K1p_exact_K2<state_datatype> IntegrandK2(Lambda, w, vpp, false, Pi);
        state_datatype val_K2 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK2, -vmax, vmax, std::abs(w/2), {vpp}, Delta);
        return -SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * val_K2;
        //return vpp*vpp;
    }
};


/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class IntegranddGammaC_exact_K3{
public:
    const double Lambda;
    const double w, v, vp;
    double vmax = 1e3;
    const double Delta = (Lambda * glb_Gamma)/2;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    IntegranddGammaC_exact_K3(double Lambda_in, double w, double v, double vp, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) {}


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=1; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i);
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i) * (frac) + this->SOPTstate.vertex.half1().avertex.frequencies.  primary_grid.get_frequency(i+1) * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not KELDYSH_FORMALISM
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = data_dir + "integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {

        K1rdot_PIa_K1p_exact_K3<state_datatype> IntegrandK3(Lambda, w, vpp, vp, false, Pi);
        state_datatype val_K3 = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandK3, -vmax, vmax, std::abs(w/2), {vpp, vp, std::abs(vpp)-std::abs(vp)}, Delta);
        return -SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * val_K3;
        //return vpp*vpp;
    }
};



// to check central part of multi-loop flow equations:
// compute diagrams with non-symmetric_full intermediate results
void compute_non_symmetric_diags(const double Lambda, bool write_flag = false, int version=1, bool compute_exact=false) {
    const int it_spin = 0;
    State<state_datatype> bare (Lambda); // bare state
    bare.initialize();         // initialize bare state

    Propagator<state_datatype> G (Lambda, bare.selfenergy, 'g'); // bare propagator
    Propagator<state_datatype> S (Lambda, bare.selfenergy, 's'); // bare differentiated propagator = single scale propagator
    const Bubble<state_datatype> Pi(G, S, false);
    double Delta = (glb_Gamma + Lambda)/2;

    // Psi := K1p in PT2 + bare vertex
    State<state_datatype> Psi (Lambda);
    Psi.initialize();         // initialize bare state
    bubble_function(Psi.vertex, bare.vertex, bare.vertex, G, G, 'p', false, stdConfig); // Psi = Gamma_0 + K1p

    // K1a_dot in PT2
    State<state_datatype> PT2_K1adot (Lambda);
    bubble_function(PT2_K1adot.vertex, bare.vertex, bare.vertex, G, S, 'a', true, stdConfig);
    // K1p_dot in PT2
    State<state_datatype> PT2_K1pdot (Lambda);
    bubble_function(PT2_K1pdot.vertex, bare.vertex, bare.vertex, G, S, 'p', true, stdConfig);

    if (write_flag) {
        write_state_to_hdf(data_dir + "Psi_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, Psi);
        write_state_to_hdf(data_dir + "PT2_K1a_dot_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1adot);
        write_state_to_hdf(data_dir + "PT2_K1p_dot_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1pdot);
    }

    std::vector<State<state_datatype>> central_bubblestates = {PT2_K1adot, PT2_K1pdot};

    //for (int i = 0; i < 2; i++){
    int i = version;
        State<state_datatype> centralstate_dot = central_bubblestates[i];

        // intermediate results
        State<state_datatype,false> K1rdot_PIa_K1p (Lambda);
        bubble_function(K1rdot_PIa_K1p.vertex, centralstate_dot.vertex, Psi.vertex, G, G, 'a', false, stdConfig);


        State<state_datatype> K1p_PIa_K1rdot (Lambda);
        bubble_function(K1p_PIa_K1rdot.vertex, Psi.vertex, centralstate_dot.vertex, G, G, 'a', false, stdConfig);


        if (write_flag) {
            write_state_to_hdf(data_dir + "K1rdot_PIa_K1p_UNREORDERED_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1rdot_PIa_K1p);
            write_state_to_hdf(data_dir + "K1p_PIa_K1rdot_UNREORDERED_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1p_PIa_K1rdot);
        }

        Vertex<state_datatype,false> dGammaL_half1 = K1rdot_PIa_K1p.vertex;
        Vertex<state_datatype,false> dGammaR_half1 = K1p_PIa_K1rdot.vertex;
#if not DEBUG_SYMMETRIES
        dGammaL_half1.half1().reorder_due2antisymmetry(dGammaR_half1.half1());
#endif
        K1rdot_PIa_K1p.vertex = dGammaL_half1;
        K1p_PIa_K1rdot.vertex = dGammaR_half1;

        // create non-symmetric_full vertex with differentiated vertex on the left
        GeneralVertex<state_datatype , non_symmetric_diffleft, false> dGammaL(Lambda);
        dGammaL.half1()  = dGammaL_half1.half1();  // assign half 1 to dGammaL
        dGammaL.half2() = dGammaR_half1.half1();  // assign half 2 as half 1 to dGammaR [symmetric_full -> left()=right()]


        // insert this non-symmetric_full vertex on the right of the bubble
        State<state_datatype> dGammaC_r(Lambda);

        dGammaL.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble
        bubble_function(dGammaC_r.vertex, Psi.vertex, dGammaL, G, G, 'a', false, stdConfig);


        // create non-symmetric_full vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<state_datatype , non_symmetric_diffright,false> dGammaR (Lambda);
        dGammaR.half1() = dGammaR_half1.half1();  // assign half 1
        dGammaR.half2() = dGammaL_half1.half1();  // assign half 2 as half 1 of dGammaL


        // insert this non-symmetric_full vertex on the left of the bubble
        State<state_datatype> dGammaC_l(Lambda);
        dGammaR.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

        bubble_function(dGammaC_l.vertex, dGammaR, Psi.vertex, G, G, 'a', false, stdConfig);


        if (write_flag) {
            write_state_to_hdf(data_dir + "K1rdot_PIa_K1p_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1rdot_PIa_K1p);
            write_state_to_hdf(data_dir + "K1p_PIa_K1rdot_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1p_PIa_K1rdot);
            write_state_to_hdf(data_dir + "dGammaC_r_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, dGammaC_r);
            write_state_to_hdf(data_dir + "dGammaC_l_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, dGammaC_l);
        }
    //}

    if (compute_exact) {
        /// Now computing vertex for version 1 (K1rdot = K1pdot):
        double vmax = 1e3;
        State<state_datatype> K1pdot_exact(Lambda);        // intermediate result: contains K2 and K3

#pragma omp parallel for schedule(dynamic) default(none) shared(K1pdot_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (size_t iflat = 0; iflat < (nBOS ); iflat++) {
            size_t it = iflat;
            //for (int it = 1; it<nBOS2-1; it++) {
            //    for (int j = 1; j<nFER2-1; j++) {
            double w;
            K1pdot_exact.vertex.avertex().K1.frequencies.get_freqs_w(w, it);
            state_datatype val_K1 = -SOPT_K1a_diff(w, Lambda);
            K1pdot_exact.vertex.pvertex().K1.setvert(val_K1, it_spin, it, 0, 0);
            //    }
        }
        write_state_to_hdf(data_dir + "K1rdot_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  K1pdot_exact);

        State<state_datatype> K1rdot_diff = PT2_K1pdot - K1pdot_exact;        // intermediate result: contains K2 and K3
        write_state_to_hdf(data_dir + "K1rdot_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  K1rdot_diff);


        State<state_datatype> K1rdot_PIa_K1p_exact(Lambda);        // intermediate result: contains K2 and K3

#pragma omp parallel for schedule(dynamic) default(none) shared(K1rdot_PIa_K1p_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (int iflat = 0; iflat < (nBOS2 ) * (nFER2 ); iflat++) {
            int i = iflat / (nFER2 );
            int j = iflat - (i ) * (nFER2 );
            double w, v;
            K1rdot_PIa_K1p_exact.vertex.avertex().K2.frequencies.get_freqs_w(w, v, i, j);
            K1rdot_PIa_K1p_exact_K2<state_datatype> IntegrandK2(Lambda, w, v, false, Pi);
            state_datatype val_K2 =
                    1. / (2 * M_PI) * integrator_Matsubara_T0(IntegrandK2, -vmax, vmax, std::abs(w / 2), {v}, Delta, true);
            K1rdot_PIa_K1p_exact.vertex.avertex().K2.setvert(val_K2, it_spin, i, j, 0, 0);
            //    }
        }

#if MAX_DIAG_CLASS>2
#pragma omp parallel for schedule(dynamic) default(none) shared(K1rdot_PIa_K1p_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (int iflat = 0; iflat < (nBOS3 ) * (nFER3 ) * (nFER3); iflat++) {
            int i = iflat / (nFER3) / (nFER3);
            int j = iflat / (nFER3) - (i) * (nFER3);
            int k = iflat - (i ) * (nFER3 ) * (nFER3 ) - (j ) * (nFER3 );
            double w, v, vp;
            K1rdot_PIa_K1p_exact.vertex.avertex().K3.frequencies.get_freqs_w(w, v, vp, i, j, k);
            K1rdot_PIa_K1p_exact_K3<state_datatype> IntegrandK3(Lambda, w, v, vp, false, Pi);
            state_datatype val_K3 = 1. / (2 * M_PI) *
                                    integrator_Matsubara_T0(IntegrandK3, -vmax, vmax, std::abs(w / 2),
                                                                               {v, vp, std::abs(w) - std::abs(vp), std::abs(w) + std::abs(vp),
                                                                                std::abs(w) - std::abs(v), std::abs(w) + std::abs(v)}, Delta, true);
            K1rdot_PIa_K1p_exact.vertex.avertex().K3.setvert(val_K3, it_spin, i, j, k, 0, 0);
            //        }
            //    }
        }
#endif
        write_state_to_hdf(data_dir + "K1rdot_PIa_K1p_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  K1rdot_PIa_K1p_exact);



        State<state_datatype> K1rdot_PIa_K1p_diff =
                K1rdot_PIa_K1p - K1rdot_PIa_K1p_exact;        // intermediate result: contains K2 and K3
        write_state_to_hdf(data_dir + "K1rdot_PIa_K1p_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  K1rdot_PIa_K1p_diff);


        utils::print("Relative maxabs difference in K1rdot_PIa_K1p --> K1: ", K1rdot_PIa_K1p_diff.vertex.half1().norm_K1(0) / K1rdot_PIa_K1p.vertex.half1().norm_K1(0) , true);
        if (MAX_DIAG_CLASS > 1) utils::print("Relative maxabs difference in K1rdot_PIa_K1p --> K2: ", K1rdot_PIa_K1p_diff.vertex.half1().norm_K2(0) / K1rdot_PIa_K1p.vertex.half1().norm_K2(0) , true );
        if (MAX_DIAG_CLASS > 2) utils::print("Relative maxabs difference in K1rdot_PIa_K1p --> K3: ", K1rdot_PIa_K1p_diff.vertex.half1().norm_K3(0) / K1rdot_PIa_K1p.vertex.half1().norm_K3(0) , true );


        State<state_datatype> K1p_PIa_K1rdot_exact(Lambda);        // intermediate result: contains K2 and K3

#if MAX_DIAG_CLASS>2
#pragma omp parallel for schedule(dynamic) default(none) shared(K1p_PIa_K1rdot_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (int iflat = 0; iflat < (nBOS3 ) * (nFER3 ) * (nFER3 ); iflat++) {
            int i = iflat / (nFER3 ) / (nFER3 );
            int j = iflat / (nFER3 ) - (i) * (nFER3 );
            int k = iflat - (i ) * (nFER3 ) * (nFER3) - j * (nFER3 );
            double w , v, vp;
            K1p_PIa_K1rdot_exact.vertex.avertex().K3.frequencies.get_freqs_w(w, v, vp, i, j, k);
            K1p_PIa_K1rdot_exact_K3<state_datatype> IntegrandK3(Lambda, w, v, vp, false, Pi);
            state_datatype val_K3 = 1. / (2 * M_PI) *
                                    integrator_Matsubara_T0(IntegrandK3, -vmax, vmax, std::abs(w / 2),
                                                                               {v, vp, std::abs(w) - std::abs(vp), std::abs(w) + std::abs(vp),
                                                                                std::abs(w) - std::abs(v), std::abs(w) + std::abs(v)}, Delta, true);
            K1p_PIa_K1rdot_exact.vertex.avertex().K3.setvert(val_K3, 0, it_spin, i, j, k, 0);
            //        }
            //    }
        }
#endif
        write_state_to_hdf(data_dir + "K1p_PIa_K1rdot_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  K1p_PIa_K1rdot_exact);

        State<state_datatype> K1p_PIa_K1rdot_diff =
                K1p_PIa_K1rdot - K1p_PIa_K1rdot_exact;        // intermediate result: contains K2 and K3
        write_state_to_hdf(data_dir + "K1p_PIa_K1rdot_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  K1p_PIa_K1rdot_diff);


        utils::print("Relative maxabs difference in K1p_PIa_K1rdot --> K1: ", K1p_PIa_K1rdot_diff.vertex.half1().norm_K1(0) / K1p_PIa_K1rdot.vertex.half1().norm_K1(0) , true);
        if (MAX_DIAG_CLASS > 1) utils::print("Relative maxabs difference in K1p_PIa_K1rdot --> K2: ", K1p_PIa_K1rdot_diff.vertex.half1().norm_K2(0) / K1p_PIa_K1rdot.vertex.half1().norm_K2(0) , true );
        if (MAX_DIAG_CLASS > 2) utils::print("Relative maxabs difference in K1p_PIa_K1rdot --> K3: ", K1p_PIa_K1rdot_diff.vertex.half1().norm_K3(0) / K1p_PIa_K1rdot.vertex.half1().norm_K3(0) , true );

        State<state_datatype> dGammaC_exact(Lambda);        // final state: contains K1, K2 and K3

#pragma omp parallel for schedule(dynamic) default(none) shared(dGammaC_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (int i = 0; i < nBOS; i++) {
            double w;
            dGammaC_exact.vertex.avertex().K1.frequencies.get_freqs_w(w, i);
            IntegranddGammaC_exact_K1<state_datatype> IntegrandK1(Lambda, w, false, Pi);
            state_datatype val_K1 =
                    1. / (2 * M_PI) * integrator_Matsubara_T0(IntegrandK1, -vmax, vmax, std::abs(w / 2), {}, Delta, true);
            dGammaC_exact.vertex.avertex().K1.setvert(val_K1, it_spin, i, 0, 0);
        }


#pragma omp parallel for schedule(dynamic) default(none) shared(dGammaC_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (int iflat = 0; iflat < (nBOS2) * (nFER2 ); iflat++) {
            int i = iflat / (nFER2 );
            int j = iflat - (i ) * (nFER2 );
            double w;
            double v;
            dGammaC_exact.vertex.avertex().K2.frequencies.get_freqs_w(w, v, i, j);
            IntegranddGammaC_exact_K2<state_datatype> IntegrandK2(Lambda, w, v, false, Pi);
            state_datatype val_K2 =
                    1. / (2 * M_PI) * integrator_Matsubara_T0(IntegrandK2, -vmax, vmax, std::abs(w / 2), {v}, Delta, true);
            dGammaC_exact.vertex.avertex().K2.setvert(val_K2, it_spin, i, j, 0, 0);
            //    }
        }

#if MAX_DIAG_CLASS>2
#pragma omp parallel for schedule(dynamic) default(none) shared(dGammaC_exact, vmax, Delta, Lambda, it_spin, Pi)
        for (int iflat = 0; iflat < (nBOS3 ) * (nFER3 ) * (nFER3 ); iflat++) {
            int i = iflat / (nFER3 ) / (nFER3 );
            int j = iflat / (nFER3 ) - (i ) * (nFER3 );
            int k = iflat - (i) * (nFER3 ) * (nFER3 ) - (j ) * (nFER3 );
            double w, v, vp;
            dGammaC_exact.vertex.avertex().K3.frequencies.get_freqs_w(w, v, vp, i, j, k);
            IntegranddGammaC_exact_K3<state_datatype> IntegrandK3(Lambda, w, v, vp, false, Pi);
            state_datatype val_K3 = 1. / (2 * M_PI) *
                                    integrator_Matsubara_T0(IntegrandK3, -vmax, vmax, std::abs(w / 2),
                                                                  {v, vp, std::abs(v) - std::abs(vp)}, Delta, true);
            dGammaC_exact.vertex.avertex().K3.setvert(val_K3, it_spin, i, j, k, 0, 0);
            //    }
            //}
        }
#endif
        write_state_to_hdf(data_dir + "dGammaC_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  dGammaC_exact);

        State<state_datatype> dGammaC_l_diff = dGammaC_l - dGammaC_exact;        // final result: contains K1, K2 and K3
        write_state_to_hdf(data_dir + "dGammaC_l_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  dGammaC_l_diff);

        State<state_datatype> dGammaC_r_diff = dGammaC_r - dGammaC_exact;        // final result: contains K1, K2 and K3
        write_state_to_hdf(data_dir + "dGammaC_r_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  dGammaC_r_diff);

        utils::print("Relative maxabs difference in dGammaC_l --> K1: ", dGammaC_l_diff.vertex.half1().norm_K1(0) / dGammaC_l.vertex.half1().norm_K1(0)  , true);
        if (MAX_DIAG_CLASS > 1) utils::print("Relative maxabs difference in dGammaC_l --> K2: ", dGammaC_l_diff.vertex.half1().norm_K2(0) / dGammaC_l.vertex.half1().norm_K2(0)  , true);
        if (MAX_DIAG_CLASS > 2) utils::print("Relative maxabs difference in dGammaC_l --> K3: ", dGammaC_l_diff.vertex.half1().norm_K3(0) / dGammaC_l.vertex.half1().norm_K3(0)  , true);

        utils::print("Relative maxabs difference in dGammaC_r --> K1: ", dGammaC_r_diff.vertex.half1().norm_K1(0) / dGammaC_r.vertex.half1().norm_K1(0)  , true);
        if (MAX_DIAG_CLASS > 1) utils::print("Relative maxabs difference in dGammaC_r --> K2: ", dGammaC_r_diff.vertex.half1().norm_K2(0) / dGammaC_r.vertex.half1().norm_K2(0)  , true);
        if (MAX_DIAG_CLASS > 2) utils::print("Relative maxabs difference in dGammaC_r --> K3: ", dGammaC_r_diff.vertex.half1().norm_K3(0) / dGammaC_r.vertex.half1().norm_K3(0)  , true);
    }
}
#else
template <typename Q, int type = 2>
auto SOPT_K1a(double w, double Lambda) -> Q {

    double Delta = REG == 2 ? (glb_Gamma + Lambda) / 2. : glb_Gamma * 0.5;
    const double hartree_term = PARTICLE_HOLE_SYMMETRY ? 0.5*glb_U : Hartree_Solver(Lambda).compute_Hartree_term_bracketing();
    const double HF_corr = glb_Vg + hartree_term - 0.5 * glb_U;

    auto advanced = [Delta,HF_corr,Lambda](const double w_) -> Q {
        const Q factor = REG == 4 ? Lambda*Lambda : 1.;
        if (w_ == 0.)
            return - factor * glb_U*glb_U * Delta / M_PI / (Delta*Delta + HF_corr*HF_corr) * 0.5;
        else
            return - factor * glb_U*glb_U * Delta/M_PI / ((glb_i * w_)*(2*Delta + (glb_i * w_))) * log(1. + ((glb_i * w_)*(2*Delta + (glb_i * w_)))/(Delta*Delta + HF_corr*HF_corr)) * 0.5;
    };

    if (type == 2) {
        /// advanced component of K1a in SOPT:
        Q result = advanced(w);
        return result;
    }
    else if (type == 1) {
        /// retarded component of K1a in SOPT:
        Q result = conj(advanced(w));
        return result;
    }
    else {
        assert(type == 0);
        /// Keldysh component of K1a in SOPT:
        const Q adv = advanced(w);
        const Q ret = conj(adv);
        const Q result = sgn(w) * (ret - adv);
        return result;
    }

}
template <typename Q, int type = 2>
auto SOPT_K1a_diff(double w, double Lambda) -> Q {
    double Delta = (glb_Gamma + Lambda) / 2.;
    const double hartree_term = Hartree_Solver(Lambda).compute_Hartree_term_bracketing();
    const double HF_corr = glb_Vg + hartree_term - 0.5 * glb_U;

    auto advanced = [Delta,HF_corr](const double w_) -> Q {
        if (std::abs(w_)  == 0) return (- glb_U*glb_U / M_PI / (Delta*Delta + HF_corr*HF_corr) * 0.5
                                        + glb_U*glb_U / M_PI / (Delta*Delta + HF_corr*HF_corr) / (Delta*Delta + HF_corr*HF_corr) * Delta * Delta )*0.5 ;
        else {
            Q term1 = -glb_U * glb_U * 0.5   / M_PI / ((glb_i * w_) * (2. * Delta + (glb_i * w_))) * log(1. + ((glb_i * w_) * (2 * Delta + (glb_i * w_))) / (Delta * Delta + HF_corr * HF_corr));
            Q term2 = +glb_U * glb_U * Delta / M_PI / ((glb_i * w_) * (2. * Delta + (glb_i * w_))) * log(1. + ((glb_i * w_) * (2 * Delta + (glb_i * w_))) / (Delta * Delta + HF_corr * HF_corr)) /
                         (2. * Delta + (glb_i * w_));
            Q term3 = -glb_U * glb_U * Delta / M_PI / ((glb_i * w_) * (2. * Delta + (glb_i * w_))) /    (1. + ((glb_i * w_) * (2 * Delta + (glb_i * w_))) / (Delta * Delta + HF_corr * HF_corr))
                         * (
                                 (glb_i * w_) / (Delta * Delta + HF_corr * HF_corr)
                                 - (glb_i * w_) * (2 * Delta + (glb_i * w_)) / (Delta * Delta + HF_corr * HF_corr) /
                                   (Delta * Delta + HF_corr * HF_corr) * Delta
                         );
            return (term1 + term2 + term3) * 0.5;
        }
    };


    if (type == 2) {
        /// advanced component of K1a in SOPT:
        Q result = advanced(w);
        return result;
    }
    else if (type == 1) {
        /// retarded component of K1a in SOPT:
        Q result = conj(advanced(w));
        return result;
    }
    else {
        assert(type == 0);
        /// Keldysh component of K1a in SOPT:
        const Q adv = advanced(w);
        const Q ret = conj(adv);
        const Q result = sgn(w) * (ret - adv);
        return result;
    }
}

template <typename Q, int type>
class Integrand_SOPT_SE {
    const double v;
    const double Lambda;
    const double Delta;
    const Propagator<Q>& barePropagator;

public:
    Integrand_SOPT_SE(const double v, const double Lambda, const Propagator<Q>& propagator) :
    v(v), Lambda(Lambda), Delta(REG == 2 ? (Lambda + glb_Gamma) * 0.5 : glb_Gamma * 0.5), barePropagator(propagator) {};

    auto value(const double vp) const -> Q {
        const Q KR = SOPT_K1a<Q,2>(vp - v, Lambda);
        const Q KA = conj(KR);
        const Q KK = SOPT_K1a<Q,0>(vp - v, Lambda);
        const Q K00 = ( KR + KA + KK);// * 0.5;
        const Q K01 = ( KR - KA - KK);// * 0.5;
        const Q K10 = (-KR + KA - KK);// * 0.5;
        const Q K11 = (-KR - KA + KK);// * 0.5;

        using buffertype_propagator = Eigen::Matrix<Q,2,2>;
        const buffertype_propagator G = barePropagator.template valsmooth_vectorized<buffertype_propagator>(vp, 0);
#if CONTOUR_BASIS == 1
        const Q G00 = G(0,0);
        const Q G01 = G(0,1);
        const Q G10 = G(1,0);
        const Q G11 = G(1,1);
#else
        const Q GA_ = G(0,1);
        const Q GR_ = G(1,0);
        const Q GK_ = G(1,1);
        const Q G00 = 0.5*( GA_ + GR_ + GK_);
        const Q G01 = 0.5*( GA_ - GR_ + GK_);
        const Q G10 = 0.5*(-GA_ + GR_ + GK_);
        const Q G11 = 0.5*(-GA_ - GR_ + GK_);
#endif

        const Q Sigma00 = K00 * G00;
        const Q Sigma01 = K01 * G01;
        const Q Sigma10 = K10 * G10;
        const Q Sigma11 = K11 * G11;

        if (type == 1) {
            /// retarded component of selfenergy:
            const Q result = (Sigma00 + Sigma01 - Sigma10 - Sigma11) * 0.5;
            return result * glb_i;
        }
        else if (type == 0){
            /// Keldysh component of selfenergy;
            const Q result = (Sigma00 - Sigma01 - Sigma10 + Sigma11) * 0.5;
            return result * glb_i;
        }
        else {
            assert(type == 3);
            /// vanishing non-causal component
            const Q result = (Sigma00 + Sigma01 + Sigma10 + Sigma11) * 0.5;
            return result * glb_i;

        }

    }



    auto operator() (const double vp) const -> Q {
        return (value(vp) + value(-vp)) * 0.5;
    }


};

template <typename Q>
void test_PT_state(std::string outputFileName, const double Lambda, const bool diff) {
#ifndef ZERO_TEMP
    assert(false);
#endif
    const int it_spin = 0;
    double vmax = 100;
#if REG==2
    double Delta = (glb_Gamma + Lambda) / 2.;
#else
    double Delta = (glb_Gamma) / 2.;
#endif
    State<Q> bareState (Lambda);bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    Bubble<Q> Pi(barePropagator, barePropagator, diff);

    const int N_iterations = 0;
    State<Q> state_cpp (Lambda);   // create final and initial state
    state_cpp.initialize();             // initialize state
    //write_state_to_hdf("PTstate_preOpt", Lambda_ini,  N_iterations+1, state_cpp);  // save the initial state to hdf5 file
    //write_state_to_hdf("PTstate_postOpt", Lambda_ini,  N_iterations+1, state_cpp);  // save the initial state to hdf5 file
    sopt_state(state_cpp, Lambda, stdConfig,diff);
    //state_cpp.vertex.half1().check_vertex_resolution();
    //state_cpp.analyze_tails();

    write_state_to_hdf(outputFileName + "_cpp", Lambda, 1, state_cpp);

    State<state_datatype> PT_state(state_cpp, Lambda);

    // compute SOPT self-energy (numerically exact)
    for (int i = 0; i<nFER; i++) {
        double v = PT_state.selfenergy.Sigma.frequencies.  primary_grid.get_frequency(i);
        Integrand_SOPT_SE<Q,1> IntegrandSE_R(v, Lambda, barePropagator);
        Q val_SE_R = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandSE_R, -vmax, vmax, std::abs(0.), {v}, Delta, true);
        PT_state.selfenergy.setself(0, i, 0, val_SE_R);

        Integrand_SOPT_SE<Q,0> IntegrandSE_K(v, Lambda, barePropagator);
        Q val_SE_K = 1./(2*M_PI) * integrator_Matsubara_T0(IntegrandSE_K, -vmax, vmax, std::abs(0.), {v}, Delta, true);
        PT_state.selfenergy.setself(1, i, 0, val_SE_K);
    }

    // get SOPT K1 (exact)
    utils::print("Compute exact solution for SOPT: \n");
#pragma omp parallel for
    for (int i = 0; i<nBOS; i++) {
        double w = PT_state.vertex.avertex().K1.frequencies.  primary_grid.get_frequency(i);

        Q KR = diff ? SOPT_K1a_diff<Q,1>(w, Lambda) : SOPT_K1a<Q,1>(w, Lambda);
        Q KA = myconj(KR);
        Q KK = diff ? SOPT_K1a_diff<Q,0>(w, Lambda) : SOPT_K1a<Q,0>(w, Lambda);
        Q K00 = ( KR + KA + KK);// * 0.5;
        Q K01 = ( KR - KA - KK);// * 0.5;
        Q K10 = (-KR + KA - KK);// * 0.5;
        Q K11 = (-KR - KA + KK);// * 0.5;
        Q val_K1 = CONTOUR_BASIS != 1 ? KA : K00;
        Q val_K1_K = CONTOUR_BASIS != 1 ? KK : K10;
        //PT_state.vertex.avertex().K1.setvert( val_K1 - val_K1*val_K1/glb_U, 0, i, 0);
        PT_state.vertex.avertex().K1.setvert( val_K1  , it_spin,  i, 0, 0);
        PT_state.vertex.avertex().K1.setvert( val_K1_K, it_spin,  i, 1, 0);

        w = PT_state.vertex.pvertex().K1.frequencies.  primary_grid.get_frequency(i);
        KR = diff ? SOPT_K1a_diff<Q,1>(w, Lambda) : SOPT_K1a<Q,1>(w, Lambda);
        KA = myconj(KR);
        KK = diff ? SOPT_K1a_diff<Q,0>(w, Lambda) : SOPT_K1a<Q,0>(w, Lambda);
        K00 = ( KR + KA + KK);// * 0.5;
        K01 = ( KR - KA - KK);// * 0.5;
        K10 = (-KR + KA - KK);// * 0.5;
        K11 = (-KR - KA + KK);// * 0.5;
        val_K1 = CONTOUR_BASIS != 1 ? KA : K00;
        val_K1_K = CONTOUR_BASIS != 1 ? KK : K01;
        //PT_state.vertex.pvertex().K1.setvert( -val_K1 - val_K1*val_K1/glb_U, 0, i, 0);
        PT_state.vertex.pvertex().K1.setvert( -val_K1  , it_spin, i, 0, 0);
        PT_state.vertex.pvertex().K1.setvert( -val_K1_K, it_spin, i, 1, 0);

        w = PT_state.vertex.tvertex().K1.frequencies.  primary_grid.get_frequency(i);
        KR = diff ? SOPT_K1a_diff<Q,1>(w, Lambda) : SOPT_K1a<Q,1>(w, Lambda);
        KA = myconj(KR);
        KK = diff ? SOPT_K1a_diff<Q,0>(w, Lambda) : SOPT_K1a<Q,0>(w, Lambda);
        K00 = ( KR + KA + KK);// * 0.5;
        K01 = ( KR - KA - KK);// * 0.5;
        K10 = (-KR + KA - KK);// * 0.5;
        K11 = (-KR - KA + KK);// * 0.5;
        val_K1 = CONTOUR_BASIS != 1 ? KA : K00;
        val_K1_K = CONTOUR_BASIS != 1 ? KK : K10;
        val_K1 = CONTOUR_BASIS != 1 ? KA : K00;
        val_K1_K = state_cpp.vertex.tvertex().K1.val(it_spin, i, 1, 0);
        PT_state.vertex.tvertex().K1.setvert( -val_K1*val_K1*2./glb_U    , it_spin, i, 0, 0);
        PT_state.vertex.tvertex().K1.setvert( val_K1_K, it_spin, i, 1, 0);
    }
    write_state_to_hdf(outputFileName + "_exact", Lambda, 1, PT_state);


    State<Q> state_diff = state_cpp - PT_state;

    write_state_to_hdf(outputFileName + "_diff", Lambda, 1, state_diff);
    utils::print("SE-difference: ", state_diff.selfenergy.Sigma.get_vec().max_norm() / PT_state.selfenergy.Sigma.get_vec().max_norm(), true);
    utils::print("K1a-difference: ", state_diff.vertex.avertex().K1.get_vec().max_norm() / PT_state.vertex.avertex().K1.get_vec().max_norm(), true);
    utils::print("K1p-difference: ", state_diff.vertex.pvertex().K1.get_vec().max_norm() / PT_state.vertex.pvertex().K1.get_vec().max_norm(), true);
    utils::print("K1t-difference: ", state_diff.vertex.tvertex().K1.get_vec().max_norm() / PT_state.vertex.tvertex().K1.get_vec().max_norm(), true);

}

#endif // SIAM Matsubara T=0


#ifdef STATIC_FEEDBACK
/**
 * Compute the right hand side of the flow equations according to Severin Jakobs' channel decomposition with
 * approximated channel feedback and modified self-energy feedback (only static level shift to avoid
 * overbroadening of spectral features)
 * @param Psi    : state at which to compute right hand side
 * @param Lambda : Lambda at which to compute right hand side
 * @return       : dPsi (right hand side of flow equation)
 */
template <typename Q>
auto rhs_channel_decomposition(const State<Q>& Psi, const double Lambda) -> State<Q> {
    State<Q> dPsi; // result

    SelfEnergy<Q> selfEnergy;
    comp static_shift = real(Psi.selfenergy.valsmooth(0, glb_mu, 0));  // only use a static level shift as self-energy
    selfEnergy.initialize(static_shift, 0.);

    Propagator<Q> G(Lambda, selfEnergy, 'g');    //Initialization of Propagator objects
    Propagator<Q> S(Lambda, selfEnergy, 's');    //Initialization of Propagator objects

    // Self-energy flow
    loop<true,0>(dPsi.selfenergy, Psi.vertex, S);  // self-energy loop

    // Vertex flow
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'a', true, stdConfig); // diff. bubble in the a-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'p', true, stdConfig); // diff. bubble in the p-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 't', true, stdConfig); // diff. bubble in the t-channel

    return dPsi;
}

/**
 * FRG flow according to Severin Jakobs' channel decomposition with approximated channel feedback
 * and modified self-energy feedback (only static level shift to avoid overbroadening of spectral features).
 * Only correct if parameter STATIC_FEEDBACK is defined.
 * @param N_ODE : number of Runge-Kutta ODE iterations
 */
template <typename Q>
void test_channel_decomposition(int N_ODE) {
    State<Q> state_ini, state_fin;   // create initial and final state
    state_ini.initialize();             // initialize initial state

//    sopt_state(state_ini, Lambda_ini);  // set initial state to SOPT (necessary if Lambda_ini is too small)

    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_channel_decomposition,
                    log_substitution, log_resubstitution, N_ODE); // compute flow

    std::string name = "test_channel_decomposition.h5";
    write_h5_rvecs(name, {"v", "Sigma_re", "Sigma_im", "Sigma_ini_re", "Sigma_ini_im"},
                         {ffreqs,
                          state_fin.selfenergy.Sigma.real(), state_fin.selfenergy.Sigma.imag(),
                          state_ini.selfenergy.Sigma.real(), state_ini.selfenergy.Sigma.imag()});

    write_state_to_hdf("channel_decomposition.h5", 0, 1, state_fin);
}
#endif

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
