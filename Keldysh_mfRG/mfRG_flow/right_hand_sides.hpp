#ifndef RIGHT_HAND_SIDES_H
#define RIGHT_HAND_SIDES_H

#include "../data_structures.hpp"            // real/complex vector classes, imag. unit
#include "../utilities/write_data2file.hpp"  // writing data into text or hdf5 files
#include "../correlation_functions/two_point/propagator.hpp"                 // propagator to perform second-order perturbation theory (SOPT)
#include "../correlation_functions/two_point/selfenergy.hpp"                 // self-energy used in SOPT
#include "../correlation_functions/state.hpp"                      // state to perform full flow
#include "../correlation_functions/four_point/vertex.hpp"                     // Vertices to put into bubbles
#include "../loop/loop.hpp"                       // compute self-energy loop
#include "../bubble/bubble_function.hpp"                    // compute vertex bubbles
#include "../parameters/master_parameters.hpp"                 // system parameters (lengths of vectors etc.)
#include "../ODE_solvers/ODE_solvers.hpp"                // ODE solvers
#include <cassert>
#include "../utilities/hdf5_routines.hpp"
#include "../utilities/util.hpp"
#include "../bubble/precalculated_bubble.hpp"
#include "../perturbation_theory_and_parquet/parquet_solver.hpp"


template <typename Q> auto rhs_n_loop_flow(const State<Q>& Psi, double Lambda) -> State<Q>;
template <typename Q> void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q,false>& PsiVertex, const Propagator<Q>& S);
template <typename Q, class Bubble_Object> void vertexOneLoopFlow(Vertex<Q,true>& dPsiVertex, const Vertex<Q,false>& Psi, const Bubble_Object& dPi);

template <typename Q> void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const GeneralVertex<Q,symmetric_r_irred,true>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G);

template <typename Q, class Bubble_Object> auto calculate_dGammaL(const GeneralVertex<Q, symmetric_r_irred,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q,true>;
template <typename Q, class Bubble_Object> auto calculate_dGammaR(const GeneralVertex<Q, symmetric_r_irred,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q,true>;
template <typename Q, class Bubble_Object> auto calculate_dGammaC_right_insertion(const Vertex<Q,false>& PsiVertex, const GeneralVertex<Q, non_symmetric_diffleft,true>& nonsymVertex,
                                                                                  const Bubble_Object& Pi) -> Vertex<Q,true>;
template <typename Q, class Bubble_Object> auto calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric_diffright,true>& nonsymVertex, const Vertex<Q,false>& PsiVertex,
                                                                                 const Bubble_Object& Pi) -> Vertex<Q,true>;


template <typename Q> bool vertexConvergedInLoops(Vertex<Q,true>& dGamma_T, Vertex<Q,true>&dGamma);
template <typename Q> bool selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG);

/// compute dSigma in SOPT
    template <typename Q> void calculate_dSigma_SOPT(SelfEnergy<Q>& Sigma_out, const State<Q>& Psi, const double Lambda) {
    SelfEnergy<Q> Sigma_H (Lambda);
    Sigma_H.initialize(Psi.selfenergy.asymp_val_R, 0);
    Propagator<Q> G0(Lambda, Sigma_H, 'g');
    Propagator<Q>dG0(Lambda, Sigma_H, 's');
    Bubble<Q> Pi (G0, G0, false);
    Bubble<Q>dPi (G0,dG0, true);

    Vertex<Q,false> Gamma0(Lambda);
    Gamma0.initialize(Psi.vertex.irred().acc(0));

    Vertex<Q,false> bubble_a(Lambda);
    Vertex<Q,false>dbubble_a(Lambda);
    bubble_a.set_frequency_grid(Psi.vertex);
    dbubble_a.set_frequency_grid(Psi.vertex);

    bubble_function( bubble_a, Gamma0, Gamma0, Pi, 'a');
    bubble_function(dbubble_a, Gamma0, Gamma0,dPi, 'a');

    SelfEnergy<Q> result1(Lambda);
    SelfEnergy<Q> result2(Lambda);

    loop<false,0>(result1, bubble_a, dG0);
    loop<false,0>(result2, dbubble_a, G0);

    Sigma_out = result1 + result2;
}

/**
 * Function to implement an n-loop flow (without Katanin substitution).
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
template <typename Q>
auto rhs_n_loop_flow(const State<Q>& Psi, const double Lambda, const int nloops_max, const vec<size_t> opt, const fRG_config& config) -> State<Q>{  //, const bool save_intermediate=false

    Psi.vertex.check_symmetries("Psi");
    assert(nloops_max>=1);
    std::string dir_str = data_dir + "intermediateResults/";
    int iteration=-1;
    int rkStep=-1;
    double t0 = utils::get_time();

    bool save = false;
    if (opt.size() > 1) {
         iteration = opt[0];
         rkStep = opt[1];
         if (rkStep==0 and iteration==0 and config.save_intermediateResults and false) {
             utils::makedir(dir_str);
             save = true;
         }
    }


    // initialize empty state with frequency grids corresponding to those in Psi:
    State<Q,true> dPsi(Psi, Lambda); // result

#ifndef STATIC_FEEDBACK
    Propagator<Q> S (Lambda, Psi.selfenergy, 's');
    Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
#else
    SelfEnergy<Q> bareSelfEnergy (Psi.selfenergy.Sigma.frequencies);
    bareSelfEnergy.initialize(glb_U/2., 0.);

    Propagator<Q> S (Lambda, bareSelfEnergy, 's');
    Propagator<Q> G (Lambda, bareSelfEnergy, 'g');
#endif

    //For flow without self-energy, comment out this line
    selfEnergyOneLoopFlow(dPsi.selfenergy, Psi.vertex, S);
    //calculate_dSigma_SOPT<Q>(dPsi.selfenergy, Psi, Lambda);
    const SelfEnergy<Q> selfenergy_1loop = dPsi.selfenergy;
    //dPsi.selfenergy.check_resolution();

    int counter_Selfenergy_iterations = 0;
    double selfenergy_correction_abs = 1e10;
    double selfenergy_correction_rel = 1.;
#if SELF_ENERGY_FLOW_CORRECTIONS != 0
    do {
#endif

#ifndef STATIC_FEEDBACK
#ifdef KATANIN
    Propagator<Q> dG (Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
#else
    Propagator<Q> dG (Lambda, Psi.selfenergy, 's');
#endif // KATANIN
    //Run alternatively, for no self-energy feedback
//    Propagator<Q> dG (Lambda, Psi.selfenergy, 's');
#else
    Propagator<Q> dG (Lambda, bareSelfEnergy, 's');
#endif // STATIC_FEEDBACK

    // Initialize bubble objects;
#ifdef HUBBARD // Use precalculated bubble in this case
    PrecalculatedBubble<comp> Pi(G, dG, false);
    PrecalculatedBubble<comp> dPi(G, dG, true);
#else // Otherwise use same type of bubble as before, which directly interpolates
    Bubble<Q> Pi(G, dG, false);
    Bubble<Q> dPi(G, dG, true);
#endif // HUBBARD

    State<Q> Psi_comp = Psi; // Input state to do computations with. Cannot be const, because it has to be modified in computations for the Hubbard model.
    // Calculate all possible cross-projections already here to save time later.
    if (HUBBARD_MODEL) Psi_comp.vertex.calculate_all_cross_projections();

    if (VERBOSE) utils::print("Compute 1-loop contribution: ", true);
    vertexOneLoopFlow(dPsi.vertex, Psi_comp.vertex, dPi);
    if(VERBOSE) {
        dPsi.vertex.half1().check_vertex_resolution();
        dPsi.analyze_tails();
        compare_with_FDTs(Psi.vertex, Lambda, iteration, "Psi_RKstep"+std::to_string(rkStep), false, config.nODE_ + U_NRG.size() + 1);
        compare_with_FDTs(dPsi.vertex, Lambda, iteration, "dPsi_RKstep"+std::to_string(rkStep), false, config.nODE_ + U_NRG.size() + 1);
    }



            /// save intermediate states:
    if (save) {
        if (iteration == 0) {
            write_state_to_hdf<Q>(dir_str+ "Psi"+"_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, Psi);
            write_state_to_hdf<Q>(dir_str+"dPsi"+"_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi);

            if constexpr (DEBUG_SYMMETRIES) {
                Psi.vertex.check_symmetries("Psi");
                dPsi.vertex.check_symmetries("dPsi");
            }
        }
        else {
            add_state_to_hdf<Q>(dir_str+ "Psi_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), iteration, Psi, false);
            add_state_to_hdf<Q>(dir_str+"dPsi_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), iteration, dPsi, false);
        }
    }

    if (nloops_max>=2) {
        // Calculate left and right part of 2-loop contribution.
        // The result contains only part of the information (half 1), thus needs to be completed to a non-symmetric_full vertex
        // when inserted in the 3-loop contribution below.

        GeneralVertex<Q, symmetric_r_irred,true> dGamma_1loop(dPsi.vertex.half1());
        if (VERBOSE) utils::print("Compute dGammaL (2-loop): ", true);
        Vertex<Q,true> dGammaL_half1 = calculate_dGammaL(dGamma_1loop, Psi.vertex, Pi);
        if(VERBOSE) {
            dGammaL_half1.half1().check_vertex_resolution();
            compare_with_FDTs(dGammaL_half1, Lambda, iteration, "dGammaL_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(2), false, config.nODE_ + U_NRG.size() + 1);
        }

        if (VERBOSE) utils::print("Compute dGammaR (2-loop):", true);
        Vertex<Q,true> dGammaR_half1 = calculate_dGammaR(dGamma_1loop, Psi.vertex, Pi);
        if(VERBOSE) {
            dGammaR_half1.half1().check_vertex_resolution();
        }
        Vertex<Q,true> dGammaT =
                dGammaL_half1 + dGammaR_half1; // since sum dGammaL + dGammaR is symmetric_full, half 1 is sufficient
        dPsi.vertex += dGammaT;

    if constexpr(not DEBUG_SYMMETRIES) {
        // subdiagrams don't fulfill the full symmetry of the vertex
        // the symmetry-related diagram with a differentiated vertex on the left might be one with differentiated vertex on the right (vice versa)
        // for further evaluation as part of a bigger diagram they need to be reordered to recover the correct dGammaL and dGammaR
        // acc. to symmetry relations (enforce_symmetry() assumes full symmetry)
        dGammaL_half1.half1().reorder_due2antisymmetry(dGammaR_half1.half1());
    }


        /// save intermediate states:
        if (save) {
            int i = 2;
            State<Q,true> dPsi_L(dGammaL_half1, dPsi.selfenergy, Lambda);
            State<Q,true> dPsi_R(dGammaR_half1, dPsi.selfenergy, Lambda);
            State<Q,true> dPsi_T(dGammaT, dPsi.selfenergy, Lambda);
            if (iteration == 0) {
                write_state_to_hdf<Q>(dir_str + "dPsi_L" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" + std::to_string(i),
                             Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_L);
                write_state_to_hdf<Q>(dir_str + "dPsi_R" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" + std::to_string(i),
                             Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_R);
                write_state_to_hdf<Q>(dir_str + "dPsi_T" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" + std::to_string(i),
                             Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_T);
                if constexpr(DEBUG_SYMMETRIES) {
                    // dPsi_L.vertex.check_symmetries("dPsi_L"); // we don't expect these to have the full symmetry of the vertex
                    // dPsi_R.vertex.check_symmetries("dPsi_R"); // we don't expect these to have the full symmetry of the vertex
                    dPsi_T.vertex.check_symmetries("dPsi_T");
                }
            } else {
                add_state_to_hdf<Q>(dir_str + "dPsi_L" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" + std::to_string(i), iteration, dPsi_L, false);
                add_state_to_hdf<Q>(dir_str + "dPsi_R" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" + std::to_string(i), iteration, dPsi_R, false);
                add_state_to_hdf<Q>(dir_str + "dPsi_T" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" + std::to_string(i), iteration, dPsi_T, false);

            }
        }

        if (nloops_max >= 3) {

#if SELF_ENERGY_FLOW_CORRECTIONS != 0
            // initialize central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
            GeneralVertex<Q, symmetric_r_irred,true> dGammaC_tbar(Lambda,Psi.vertex);
            dGammaC_tbar.set_frequency_grid(Psi.vertex);
            //dGammaC_tbar.set_Ir(true);
#endif

            double abs_loop = 1;
            double rel_loop = 1;
            for (int i = 3; i <= nloops_max; i++) {


                // create non-symmetric_full vertex with differentiated vertex on the left (full dGammaL, containing half 1 and 2)
                // assign half 1 to dGammaL
                // assign half 2 as half 1 of dGammaR [symmetric_full -> half1()=half2()]
                GeneralVertex<Q, non_symmetric_diffleft,true> dGammaL(dGammaL_half1.half1(), dGammaR_half1.half1(), Psi.vertex);
                if (VERBOSE) {
                    compare_with_FDTs(dGammaL, Lambda, iteration,
                                      "dGammaL_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(2), false,
                                      config.nODE_ + U_NRG.size() + 1);
                }


                // insert this non-symmetric_full vertex on the right of the bubble
                if (VERBOSE) utils::print("Compute dGammaC (right insertion) ( ", i, "-loop): \n");
                //save_symmetry_expanded_vertex<false>(dGammaL, "dPsiLIn3LoopCRight");
                Vertex<Q,true> dGammaC_r = calculate_dGammaC_right_insertion(Psi.vertex, dGammaL, Pi);
                if (VERBOSE) dGammaC_r.half1().check_vertex_resolution();
                if (VERBOSE) {
                    compare_with_FDTs(dGammaC_r, Lambda, iteration,
                                      "dGammaC_r_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(2),
                                      false, config.nODE_ + U_NRG.size() + 1);
                }

                // create non-symmetric_full vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
                // assign half 1
                // assign half 2 as half 1 of dGammaL
                GeneralVertex<Q, non_symmetric_diffright, true> dGammaR(dGammaR_half1.half1(), dGammaL_half1.half1(), Psi.vertex);
                if (VERBOSE) {
                    compare_with_FDTs(dGammaR, Lambda, iteration,
                                      "dGammaR_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(2), false,
                                      config.nODE_ + U_NRG.size() + 1);
                }


                // insert this non-symmetric_full vertex on the left of the bubble
                //save_symmetry_expanded_vertex<true>(dGammaR, "dPsiRIn3LoopCLeft");
                Vertex<Q,true> dGammaC_l = calculate_dGammaC_left_insertion(dGammaR, Psi.vertex, Pi);
                if (VERBOSE) {
                    compare_with_FDTs(dGammaC_l, Lambda, iteration,
                                      "dGammaC_l_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(2),
                                      false, config.nODE_ + U_NRG.size() + 1);
                }

                // symmetrize by averaging left and right insertion
                for (char r : {'a', 'p', 't'}) {
                    dGammaC_l.get_rvertex(r).K2 = dGammaC_r.get_rvertex(r).K2;
                }
                if constexpr(DEBUG_SYMMETRIES) {
                    for (char r: {'a', 'p', 't'}) {
                        dGammaC_r.get_rvertex(r).K2b = dGammaC_l.get_rvertex(r).K2b;
                    }
                }
                //for (char r : {'a', 'p', 't'}) {
                //    dGammaC_l.get_rvertex(r).K3 *= 0.;
                //    dGammaC_r.get_rvertex(r).K3 *= 0.;
                //}
                Vertex<Q,true> dGammaC = (dGammaC_r + dGammaC_l) * 0.5; //dGammaC_r; //                  /// TODO: Find better solution --> K2 in dGammaC_l is bad, K3 in dGammaC_l can be obtained from dGammaC_r


                /// save intermediate states:
                if (save) {
                    State<Q,true> dPsi_C (dGammaC, dPsi.selfenergy, Lambda);
                    State<Q,true> dPsi_C_left (dGammaC_l, dPsi.selfenergy, Lambda);
                    State<Q,true> dPsi_C_right(dGammaC_r, dPsi.selfenergy, Lambda);
                    if (iteration == 0) {
                        write_state_to_hdf<Q>(dir_str+"dPsi_C"+"_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations)+"_forLoop"+std::to_string(i),       Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_C);
                        write_state_to_hdf<Q>(dir_str+"dPsi_C_left"+"_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations)+"_forLoop"+std::to_string(i),  Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_C_left);
                        write_state_to_hdf<Q>(dir_str+"dPsi_C_right"+"_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations)+"_forLoop"+std::to_string(i), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_C_right);
                        if constexpr(DEBUG_SYMMETRIES) {
                            //dPsi_C_left.vertex.check_symmetries("dPsi_C_left");   // we don't expect these to have the full symmetry of the vertex
                            //dPsi_C_right.vertex.check_symmetries("dPsi_C_right"); // we don't expect these to have the full symmetry of the vertex
                            dPsi_C.vertex.check_symmetries("dPsi_C");
                        }
                    } else {
                        add_state_to_hdf<Q>(dir_str + "dPsi_C" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                            std::to_string(i), iteration, dPsi_C, false);
                        //add_state_to_hdf<Q>(dir_str+"dPsi_C_left"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i),  iteration, dPsi_C_left,  false);
                        //add_state_to_hdf<Q>(dir_str+"dPsi_C_right"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), iteration, dPsi_C_right, false);

                    }
                    if (VERBOSE) {
                        utils::print("Analyze tails of dPsi_C");
                        dPsi_C.analyze_tails();
                    }
                }


                if (VERBOSE) utils::print("Compute dGammaL ( ", i, "-loop): \n");
                //save_symmetry_expanded_vertex<true>(dGamma_1loop, "dPsiIn3LoopLLeft");
                dGammaL_half1 = calculate_dGammaL(GeneralVertex<Q, symmetric_r_irred,true>(dGammaT.half1()), Psi.vertex, Pi);
                if (VERBOSE) dGammaL_half1.half1().check_vertex_resolution();
                if (VERBOSE) utils::print("Compute dGammaR ( ", i, "-loop): \n");
                //save_symmetry_expanded_vertex<true>(dGamma_1loop, "dPsiIn3LoopRRight");
                dGammaR_half1 = calculate_dGammaR(GeneralVertex<Q, symmetric_r_irred,true>(dGammaT.half1()), Psi.vertex, Pi);
                if (VERBOSE) dGammaR_half1.half1().check_vertex_resolution();

                dGammaT = dGammaL_half1 + dGammaC +
                          dGammaR_half1; // since sum dGammaL + dGammaR is symmetric_full, half 1 is sufficient
                dPsi.vertex += dGammaT;
                if constexpr (not DEBUG_SYMMETRIES) {
                    // subdiagrams don't fulfill the full symmetry of the vertex
                    // the symmetry-related diagram with a differentiated vertex on the left might be one with differentiated vertex on the right (vice versa)
                    // for further evaluation as part of a bigger diagram they need to be reordered to recover the correct dGammaL and dGammaR
                    // acc. to symmetry relations (enforce_symmetry() assumes full symmetry)
                    dGammaL_half1.half1().reorder_due2antisymmetry(dGammaR_half1.half1());
                }
#if SELF_ENERGY_FLOW_CORRECTIONS != 0
                // extract central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
                //Vertex<Q,true> dGammaC_ap(Lambda);                   // initialize new vertex
                //dGammaC_ap.set_frequency_grid(Psi.vertex);
                //dGammaC_ap.avertex() = dGammaC.avertex();  // copy results from calculations above
                //dGammaC_ap.pvertex() = dGammaC.pvertex();
                dGammaC_tbar += GeneralVertex<Q, symmetric_r_irred,true>(dGammaC.half1());                      // add the i-loop contribution to the full dGammaC_tbar
#endif
                //if(vertexConvergedInLoops(dGammaT, dPsi.vertex))
                //    break;


                /// save intermediate states:
                if (save) {
                    State<Q,true> dPsi_L(dGammaL_half1, dPsi.selfenergy, Lambda);
                    State<Q,true> dPsi_R(dGammaR_half1, dPsi.selfenergy, Lambda);
                    State<Q,true> dPsi_T(dGammaT, dPsi.selfenergy, Lambda);
                    if (iteration == 0) {
                        write_state_to_hdf<Q>(dir_str + "dPsi_L" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                              std::to_string(i), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_L);
                        write_state_to_hdf<Q>(dir_str + "dPsi_R" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                              std::to_string(i), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_R);
                        write_state_to_hdf<Q>(dir_str + "dPsi_T" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                              std::to_string(i), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_T);
                        if constexpr(DEBUG_SYMMETRIES) {
                            // dPsi_L.vertex.check_symmetries("dPsi_L"); // we don't expect these to have the full symmetry of the vertex
                            // dPsi_R.vertex.check_symmetries("dPsi_R"); // we don't expect these to have the full symmetry of the vertex
                            dPsi_T.vertex.check_symmetries("dPsi_T");
                        }
                    } else {
                        add_state_to_hdf<Q>(dir_str + "dPsi_L" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                            std::to_string(i), iteration, dPsi_L, false);
                        add_state_to_hdf<Q>(dir_str + "dPsi_R" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                            std::to_string(i), iteration, dPsi_R, false);
                        add_state_to_hdf<Q>(dir_str + "dPsi_T" + "_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations) + "_forLoop" +
                                            std::to_string(i), iteration, dPsi_T, false);

                    }
                    if (VERBOSE) {
                        utils::print("Analyze tails of dPsi_L");
                        dPsi_L.analyze_tails();
                        utils::print("Analyze tails of dPsi_R");
                        dPsi_R.analyze_tails();
                        utils::print("Analyze tails of dPsi_T");
                        dPsi_T.analyze_tails();
                    }
                }


                abs_loop = dGammaT.norm();
                rel_loop = dGammaT.norm() / dPsi.vertex.norm();

                utils::print("Rel. contribution from loop ", i, ": ", rel_loop*100, " %\n");
                utils::print("Abs. contribution from loop ", i, ": ", abs_loop, "\n");
                if (abs_loop < loop_tol_abs or rel_loop < loop_tol_rel) break;

            }

#if SELF_ENERGY_FLOW_CORRECTIONS == 1
            // compute multiloop corrections to self-energy flow
            State<Q,true> Psi_SEcorrection(dPsi.vertex, dPsi.selfenergy, Lambda);
            selfEnergyFlowCorrections(Psi_SEcorrection.selfenergy, dGammaC_tbar, Psi, G); // isolated SE correction
            SelfEnergy<Q> selfEnergy_old = dPsi.selfenergy;
            dPsi.selfenergy = selfenergy_1loop + Psi_SEcorrection.selfenergy;
            SelfEnergy<Q> selfEnergy_err = dPsi.selfenergy - selfEnergy_old;

            selfenergy_correction_abs = selfEnergy_err.norm();
            selfenergy_correction_rel = selfEnergy_err.norm() / dPsi.selfenergy.norm();
            State<Q> state_self_diffrel(Lambda);
            state_self_diffrel.selfenergy = selfEnergy_err * (1./ dPsi.selfenergy.norm());


            /// save intermediate states:
            if (save) {
                State<Q,true> dPsi_C_tbar(Vertex<Q,true>(dGammaC_tbar.half1()), dPsi.selfenergy, Lambda);
                if (iteration == 0) {
                    //write_state_to_hdf<Q>(dir_str+"dPsi_C_tbar_RKstep"+std::to_string(rkStep), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_C_tbar);
                    write_state_to_hdf<Q>(dir_str + "SE_correction_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), Psi.Lambda,
                                          config.nODE_ + U_NRG.size() + 1, Psi_SEcorrection);
                    write_state_to_hdf<Q>(dir_str + "SE_correction_rel_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), Psi.Lambda,
                                          config.nODE_ + U_NRG.size() + 1, state_self_diffrel);

                } else {
                    //add_state_to_hdf<Q>(dir_str+"dPsi_C_tbar_RKstep"+std::to_string(rkStep), iteration, dPsi_C_tbar, false);
                    add_state_to_hdf<Q>(dir_str + "SE_correction_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), iteration,
                                        Psi_SEcorrection, false);
                    add_state_to_hdf<Q>(dir_str + "SE_correction_rel_RKstep" + std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), iteration,
                                        state_self_diffrel, false);

#if DEBUG_SYMMETRIES
                    add_state_to_hdf<Q>(dir_str+"SE_correction_RKstep"+std::to_string(rkStep)+std::to_string(counter_Selfenergy_iterations), iteration, Psi_SEcorrection);
#endif

                }
            }
            //TODO(low): Implement self-energy iterations (see lines 37-39 of pseudo-code).
            //if(selfEnergyConverged(Psi.selfenergy, Lambda))
            //    break;
#elif SELF_ENERGY_FLOW_CORRECTIONS == 2
            /// compute new estimate for dSigma via Schwinger-Dyson equation
            utils::print("  ---  Computing dSigma from the SDE  --- \n");
            State<Q,true> Psi_SDE_diff_documentation = dPsi;
            double tol = 1.e-3;
            double diff = 1;
            int iteration_SDE = 0;
            do {
                const SelfEnergy<Q> dSigma_SDE = compute_diff_SDE<Q>(Psi, dPsi.vertex, dPsi.selfenergy)   ;

                const SelfEnergy<Q> selfenergy_diff = dPsi.selfenergy - dSigma_SDE;

                dPsi.selfenergy = dSigma_SDE;
                const double norm_compare =  dPsi.selfenergy.norm(0);
                diff = selfenergy_diff.norm(0) / (norm_compare < 1e-10 ? 1 : norm_compare);

                utils::print("SDE iteration ", iteration_SDE ," \t relative change : ", diff, "\n");
                iteration_SDE++;
            }
            while (diff > tol and iteration_SDE < 4);

            // Compute difference of old dSigma and new dSigma
            Psi_SDE_diff_documentation.selfenergy -= dPsi.selfenergy;

            selfenergy_correction_abs = Psi_SDE_diff_documentation.norm();
            selfenergy_correction_rel = Psi_SDE_diff_documentation.norm() / dPsi.selfenergy.norm();

            /// save intermediate states:
            if (save) {
                if (iteration == 0) {
                    //write_state_to_hdf<Q>(dir_str+"dPsi_C_tbar_RKstep"+std::to_string(rkStep), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, dPsi_C_tbar);
                    write_state_to_hdf<Q>(dir_str+"SE_SDE_diff_RKstep"+std::to_string(rkStep), Psi.Lambda, config.nODE_ + U_NRG.size() + 1, Psi_SDE_diff_documentation);

                }
                else {
                    //add_state_to_hdf<Q>(dir_str+"dPsi_C_tbar_RKstep"+std::to_string(rkStep), iteration, dPsi_C_tbar, false);
                    add_state_to_hdf<Q>(dir_str+"SE_SDE_diff_RKstep"+std::to_string(rkStep), iteration, Psi_SDE_diff_documentation, false);

                }
            }
            //TODO(low): Implement self-energy iterations (see lines 37-39 of pseudo-code).
            //if(selfEnergyConverged(Psi.selfenergy, Lambda))
            //    break;
#endif

            utils::print("Self-energy iteration ", counter_Selfenergy_iterations, "\n");
            utils::print("Self-energy correction (abs.): ", selfenergy_correction_abs, "\n");
            utils::print("Self-energy correction (rel.): ", selfenergy_correction_rel*100, " %\n");

            counter_Selfenergy_iterations++;
        }
//#endif

        }

#if SELF_ENERGY_FLOW_CORRECTIONS != 0
    } while (nloops_max > 2 and counter_Selfenergy_iterations < nmax_Selfenergy_iterations and selfenergy_correction_abs > tol_selfenergy_correction_abs and selfenergy_correction_rel > tol_selfenergy_correction_rel);
#endif

    if (true and mpi_world_rank() == 0) { //
        utils::print("Time needed for evaluation of RHS --> ");
        utils::get_time(t0); // measure time for one iteration
    }


    State<Q,false> dPsi_return(Vertex<Q,false>(dPsi.vertex.half1()), dPsi.selfenergy, Lambda);
    return dPsi_return;

}

template <typename Q>
class rhs_n_loop_flow_t {
public:
    mutable std::size_t rk_step = 0;
    mutable std::size_t iteration = 0;
    const fRG_config& frgConfig;
    const int nloops = frgConfig.nloops;

    rhs_n_loop_flow_t(const fRG_config& config) : frgConfig(config) {};

    void operator() (const State<Q>& Psi, State<Q>& dState_dLambda,  const double Lambda) const {
        dState_dLambda = rhs_n_loop_flow(Psi, Lambda, nloops, vec<size_t>({iteration, rk_step}), frgConfig);
        rk_step++;
    }
};

template <typename Q>
void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q,false>& PsiVertex, const Propagator<Q>& S){
    // Self-energy flow
    loop<true,0>(dPsiSelfEnergy, PsiVertex, S); // Loop for the Self-Energy calculation
}

template <typename Q, class Bubble_Object>
void vertexOneLoopFlow(Vertex<Q,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& dPi){
    // Vertex flow
    if constexpr (SBE_DECOMPOSITION and MAX_DIAG_CLASS==2 and USE_NEW_MFRG_EQS) {
        Vertex<Q,false> One = PsiVertex;
        for (char r: {'a', 'p','t'}) {
            One.get_rvertex(r).K2 *= 0;
        }
        dPsiVertex *= 0;
        for (char r: "apt") {
            bubble_function(dPsiVertex, PsiVertex, One, dPi, r);  // Differentiated bubble in channel r \in {a, p, t}
            bubble_function(dPsiVertex, One, PsiVertex, dPi, r);  // Differentiated bubble in channel r \in {a, p, t}
            Vertex<Q,true> dPsiVertex_temp = dPsiVertex;
            dPsiVertex_temp *= 0.;
            bubble_function(dPsiVertex_temp, One, One, dPi, r);  // Differentiated bubble in channel r \in {a, p, t}
            dPsiVertex -= dPsiVertex_temp;

        }
    }
    else {
        for (char r: "apt") {
            bubble_function(dPsiVertex, PsiVertex, PsiVertex, dPi,
                            r);  // Differentiated bubble in channel r \in {a, p, t}
        }
    }
}



/// for dGammaR in channel r we combine:     dgammaL_r = dgamma_rbar ◯ Pi_r ◯ Gamma
/// Hence dgamma_rbar only returns values of the r-IRreducible sector
template <typename Q, class Bubble_Object>
auto calculate_dGammaL(const GeneralVertex<Q, symmetric_r_irred,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q,true>{
    Vertex<Q,true> dGammaL(Lambda_ini, PsiVertex);
    dGammaL.set_frequency_grid(PsiVertex);

    //Vertex<Q,false> dPsiVertex_calc = dPsiVertex;
    //dPsiVertex_calc.set_Ir(true); // Only use the r-irreducible part
    //Vertex<Q> PsiVertex_calc = PsiVertex;
    //for (char r : {'a', 'p', 't'}) {
    //    PsiVertex_calc.get_rvertex(r).K1 *= 0;
    //    PsiVertex_calc.get_rvertex(r).K2 *= 0;
    //    if (DEBUG_SYMMETRIES) PsiVertex_calc.get_rvertex(r).K2b *= 0;
    //}

    if constexpr (SBE_DECOMPOSITION and MAX_DIAG_CLASS==2 and USE_NEW_MFRG_EQS) {
        Vertex<Q,false> One = PsiVertex;
        for (char r: {'a', 'p','t'}) {
            One.get_rvertex(r).K2 *= 0;
        }
        for (char r: "apt") {
            bubble_function(dGammaL, dPsiVertex, One, Pi, r);
        }
    }
    else {
        for (char r: "apt") {
            bubble_function(dGammaL, dPsiVertex, PsiVertex, Pi, r);
        }
    }

    return dGammaL;
}
/// for dGammaR in channel r we combine:     dgammaR_r = Gamma ◯ Pi_r ◯ dgamma_rbar
/// Hence dgamma_rbar only returns values of the r-IRreducible sector
template <typename Q, class Bubble_Object>
auto calculate_dGammaR(const GeneralVertex<Q, symmetric_r_irred,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q,true>{
    Vertex<Q,true> dGammaR(Lambda_ini, PsiVertex);
    dGammaR.set_frequency_grid(PsiVertex);

    //Vertex<Q,true> dPsiVertex_calc = dPsiVertex;
    //dPsiVertex_calc.set_Ir(true); // Only use the r-irreducible part

    if constexpr (SBE_DECOMPOSITION and MAX_DIAG_CLASS==2 and USE_NEW_MFRG_EQS) {
        Vertex<Q,false> One = PsiVertex;
        for (char r: {'a', 'p','t'}) {
            One.get_rvertex(r).K2 *= 0;
        }
        for (char r: "apt") {
            bubble_function(dGammaR, One, dPsiVertex, Pi, r);
        }
    }
    else {
        for (char r: "apt") {
            bubble_function(dGammaR, PsiVertex, dPsiVertex, Pi, r);
        }
    }

    return dGammaR;
}



template <typename Q>
class oneBubble {

public:
    const Propagator<Q>& g;
    bool diff = false;
    /**
     * Constructor:
     * @param propagatorG : first propagator (always a standard one)
     * @param propagatorS : second propagator (standard or single-scale/differentiated, depending on "diff_in")
     * @param diff_in      : whether to compute standard (false) or differentiated (true) bubble
     */
    oneBubble(const Propagator<Q>& g_in) : g(g_in) {};

    /**
     * Call operator:
     * @param iK    : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param v1    : frequency of first propagator
     * @param v2    : frequency of second propagator
     * @param i_in  : internal structure index
     * @return comp : value of the bubble evaluated at (iK, v1, v2)
     */
    auto value(int iK, double v1, double v2, int i_in) const -> Q{
        return 1.;
    }

    template<char ch_bubble>
    Q value_Matsubara(const double w, const double vpp, int i_in) const {
        return 1.;
    }

    /**
     * Wrapper for value function above, providing the natural arguments for evaluation of the bubble in each channel:
     * @param iK      : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param w       : bubble transfer frequency of the corresponding channel
     * @param vpp     : bubble integration frequency of the corresponding channel
     * @param i_in    : internal structure index
     * @param channel : channel to which the bubble belongs
     * @return Q      : value of the bubble evaluated at the arguments described above (usually comp)
     */
    auto value(int iK, double w, double vpp, int i_in, char channel) const -> Q {
        return 1.;
    }

    template<char ch_bubble>
    auto convert_to_fermionic_frequencies_1(double vpp, double w) const -> double {
        if constexpr(KELDYSH || ZERO_T) {
            if constexpr(ch_bubble == 'p') return vpp + w * 0.5;
            else return vpp - w * 0.5;
        }
        else { // Matsubata zeroT
            if constexpr(ch_bubble == 'p') return vpp + ceil2bfreq(w * 0.5);
            else return vpp - floor2bfreq(w * 0.5);
        }
    }
    template<char ch_bubble>
    auto convert_to_fermionic_frequencies_2(double vpp, double w) const -> double {
        if constexpr(KELDYSH || ZERO_T) {
            if constexpr(ch_bubble == 'p') return w * 0.5 - vpp;
            else return vpp + w * 0.5;
        }
        else { // Matsubata zeroT
            if constexpr(ch_bubble == 'p') return floor2bfreq(w * 0.5) - vpp;
            else return vpp + ceil2bfreq(w * 0.5);
        }
    }


    template<char ch_bubble>
    auto value_vectorized(const double w, const double vpp, const int i_in) const {
        return 1.;

    }

};


/// for dGammaC in channel r we combine:     dgammaC_r = Gamma ◯ Pi_r ◯ dGammaL_r
/// Hence dGammaL only returns values of the r-reducible sector
template <typename Q, class Bubble_Object>
auto calculate_dGammaC_right_insertion(const Vertex<Q,false>& PsiVertex, const GeneralVertex<Q, non_symmetric_diffleft,true>& nonsymVertex,
                                       const Bubble_Object& Pi) -> Vertex<Q,true> {
    Vertex<Q,true> dGammaC (Lambda_ini, PsiVertex);
    dGammaC.set_frequency_grid(PsiVertex);

    //GeneralVertex<Q, non_symmetric_diffleft> nonsymVertex_calc = nonsymVertex;
    //for (char r: {'a', 'p', 't'}) {
    //    nonsymVertex_calc.get_rvertex(r).K2 *= 0;
    //    if (DEBUG_SYMMETRIES) nonsymVertex_calc.get_rvertex(r).K2b *= 0;
    //    if (r != 'p')
    //        nonsymVertex_calc.get_rvertex(r).K3 *= 0;
    //}
    ////nonsymVertex_calc.get_rvertex('p').K3 += 1;
    ////nonsymVertex_calc.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    if constexpr (SBE_DECOMPOSITION and MAX_DIAG_CLASS==2 and USE_NEW_MFRG_EQS) {
        Vertex<Q,false> One = PsiVertex;
        for (char r: {'a', 'p','t'}) {
            One.get_rvertex(r).K2 *= 0;
        }
        for (char r: "apt") {
            bubble_function(dGammaC, One, nonsymVertex, Pi, r);
        }
    }
    else {
        for (char r: "apt") {
            bubble_function(dGammaC, PsiVertex, nonsymVertex, Pi, r);
        }
    }

    return dGammaC;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric_diffright,true>& nonsymVertex, const Vertex<Q,false>& PsiVertex,
                                      const Bubble_Object& Pi) -> Vertex<Q,true> {
    Vertex<Q,true> dGammaC (Lambda_ini, PsiVertex);
    dGammaC.set_frequency_grid(PsiVertex);

    //GeneralVertex<Q, non_symmetric_diffright> nonsymVertex_calc = nonsymVertex;
    //for (char r: {'a', 'p', 't'}) {
    //    nonsymVertex_calc.get_rvertex(r).K2 *= 0;
    //    if (DEBUG_SYMMETRIES) nonsymVertex_calc.get_rvertex(r).K2b *= 0;
    //    if (r != 'p')
    //        nonsymVertex_calc.get_rvertex(r).K3 *= 0;
    //}
    ////nonsymVertex_calc.get_rvertex('p').K3 += 1;
    //Vertex<Q> PsiVertex_calc = PsiVertex;
    //for (char r: {'a', 'p', 't'}) {
    //    if (r != 'a')
    //        PsiVertex_calc.get_rvertex(r).K1 *= 0;
    //    PsiVertex_calc.get_rvertex(r).K2 *= 0;
    //    if (DEBUG_SYMMETRIES) PsiVertex_calc.get_rvertex(r).K2b *= 0;
    //    //if (r != 'p')
    //        //PsiVertex_calc.get_rvertex(r).K3 *= 0;
    //}
    ////PsiVertex_calc.get_rvertex('a').K1 += 1;
    ////nonsymVertex_calc.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    if constexpr (SBE_DECOMPOSITION and MAX_DIAG_CLASS==2 and USE_NEW_MFRG_EQS) {
        Vertex<Q,false> One = PsiVertex;
        for (char r: {'a', 'p','t'}) {
            One.get_rvertex(r).K2 *= 0;
        }
        for (char r: "apt") {
            bubble_function(dGammaC, nonsymVertex, One, Pi, r);
        }
    }
    else {
        for (char r: "apt") {
            bubble_function(dGammaC, nonsymVertex, PsiVertex, Pi, r);
        }
    }

    return dGammaC;
}

// compute multiloop corrections to self-energy flow
/// The first selfenergy correction is given by
///     dSigma_tbar = G * dGamma_tbar     => dGamma_tbar has to be irreducible in the t-channel
/// The second selfenergy correction is
///     dSigma_t = (G dSigma_tbar G) * Gamma
template <typename Q>
void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const GeneralVertex<Q, symmetric_r_irred,true>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G){
    // TODO(low): also implement self-energy flow via differentiated SDE
    // TODO(low): iterate self-energy corrections (feedback of full SE flow into vertex flow etc.)?

    SelfEnergy<Q> dSigma_tbar(Psi.selfenergy.Sigma.frequencies);
    SelfEnergy<Q> dSigma_t(Psi.selfenergy.Sigma.frequencies);

    // compute first multiloop correction to self-energy flow, irreducible in the t channel
    loop<true,0>(dSigma_tbar, dGammaC_tbar, G);

    // compute second multiloop correction to self-energy flow, reducible in the t channel
    Propagator<Q> extension (G.Lambda, Psi.selfenergy, dSigma_tbar, 'e');
    loop<true,0>(dSigma_t, Psi.vertex, extension);

    dPsiSelfEnergy = dSigma_tbar + dSigma_t;

}

template <typename Q>
auto vertexConvergedInLoops(Vertex<Q,true>& dGamma_T, Vertex<Q,true>&dGamma) -> bool {
    return (dGamma_T.norm() / dGamma.norm() < converged_tol);
}
template <typename Q>
auto selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG) -> bool {
    Propagator<Q> compare(dG.Lambda, PsiSelfEnergy, dPsiSelfEnergy, 'k');
    compare += dG*((Q)-1.);

    return (  compare.norm()/ dG.norm() < converged_tol );
}

#endif //RIGHT_HAND_SIDES_H
