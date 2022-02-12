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


template <typename Q> auto rhs_n_loop_flow(const State<Q>& Psi, double Lambda) -> State<Q>;
template <typename Q> void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator<Q>& S);
template <typename Q, class Bubble_Object> void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& Psi, const Bubble_Object& dPi);

template <typename Q> void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G);

template <typename Q, class Bubble_Object> auto calculate_dGammaL(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>;
template <typename Q, class Bubble_Object> auto calculate_dGammaR(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>;
template <typename Q, class Bubble_Object> auto calculate_dGammaC_right_insertion(const Vertex<Q>& PsiVertex, const GeneralVertex<Q, non_symmetric>& nonsymVertex,
                                                                                  const Bubble_Object& Pi) -> Vertex<Q>;
template <typename Q, class Bubble_Object> auto calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric>& nonsymVertex, const Vertex<Q>& PsiVertex,
                                                                                 const Bubble_Object& Pi) -> Vertex<Q>;


template <typename Q> bool vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma);
template <typename Q> bool selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG);


/**
 * Function to implement an n-loop flow (without Katanin substitution).
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
template <typename Q>
auto rhs_n_loop_flow(const State<Q>& Psi, const double Lambda, const vec<size_t> opt) -> State<Q>{  //, const bool save_intermediate=false

    static_assert(N_LOOPS>=1, "");
    std::string dir_str = data_dir + "intermediateResults/";
    int iteration=-1;
    int rkStep=-1;
    bool save_intermediate = false;
    double t0 = get_time();

    if (opt.size() > 1) {
         iteration = opt[0];
         rkStep = opt[1];
         save_intermediate = true;
         if (rkStep==0 and iteration==0 and save_intermediate) makedir(dir_str);
    }


    // initialize empty state with frequency grids corresponding to those in Psi:
    State<Q> dPsi(Psi, Lambda); // result

#ifndef STATIC_FEEDBACK
    Propagator<Q> S (Lambda, Psi.selfenergy, 's');
    Propagator<Q> G (Lambda, Psi.selfenergy, 'g');
#else
    SelfEnergy<Q> bareSelfEnergy (Psi.selfenergy.frequencies);
    bareSelfEnergy.initialize(glb_U/2., 0.);

    Propagator<Q> S (Lambda, bareSelfEnergy, 's');
    Propagator<Q> G (Lambda, bareSelfEnergy, 'g');
#endif

    //For flow without self-energy, comment out this line
    selfEnergyOneLoopFlow(dPsi.selfenergy, Psi.vertex, S);
    //dPsi.selfenergy.check_resolution();

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

    if (VERBOSE) print("Compute 1-loop contribution: ", true);
    vertexOneLoopFlow(dPsi.vertex, Psi_comp.vertex, dPi);
    if(VERBOSE) {
        dPsi.vertex.half1().check_vertex_resolution();
        dPsi.analyze_tails();
        compare_with_FDTs(Psi.vertex, Lambda, iteration, "Psi_RKstep"+std::to_string(rkStep), false, nLambda_layers);
        compare_with_FDTs(dPsi.vertex, Lambda, iteration, "dPsi_RKstep"+std::to_string(rkStep), false, nLambda_layers);
    }



            /// save intermediate states:
    if (save_intermediate) {
        if (iteration == 0) {
            write_state_to_hdf<Q>(dir_str+ "Psi"+"_RKstep"+std::to_string(rkStep), Psi.Lambda, nLambda_layers, Psi);
            write_state_to_hdf<Q>(dir_str+"dPsi"+"_RKstep"+std::to_string(rkStep), Psi.Lambda, nLambda_layers, dPsi);

#ifdef DEBUG_SYMMETRIES
            Psi.vertex.check_symmetries("Psi");
            dPsi.vertex.check_symmetries("dPsi");
#endif
        }
        else {
            add_state_to_hdf<Q>(dir_str+ "Psi_RKstep"+std::to_string(rkStep), iteration, Psi, false);
            add_state_to_hdf<Q>(dir_str+"dPsi_RKstep"+std::to_string(rkStep), iteration, dPsi, false);
        }
    }

    if (N_LOOPS>=2) {
        // Calculate left and right part of 2-loop contribution.
        // The result contains only part of the information (half 1), thus needs to be completed to a non-symmetric vertex
        // when inserted in the 3-loop contribution below.

        if (VERBOSE) print("Compute dGammaL (2-loop): ", true);
        Vertex<Q> dGammaL_half1 = calculate_dGammaL(dPsi.vertex, Psi.vertex, Pi);
        if(VERBOSE) {
            dGammaL_half1.half1().check_vertex_resolution();
            compare_with_FDTs(dGammaL_half1, Lambda, iteration, "dGammaL_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(2), false, nLambda_layers);
        }

        if (VERBOSE) print("Compute dGammaR (2-loop):", true);
        Vertex<Q> dGammaR_half1 = calculate_dGammaR(dPsi.vertex, Psi.vertex, Pi);
        if(VERBOSE) {
            dGammaR_half1.half1().check_vertex_resolution();
        }
        Vertex<Q> dGammaT =
                dGammaL_half1 + dGammaR_half1; // since sum dGammaL + dGammaR is symmetric, half 1 is sufficient
        dPsi.vertex += dGammaT;

#ifndef DEBUG_SYMMETRIES
        // subdiagrams don't fulfill the full symmetry of the vertex
        // the symmetry-related diagram with a differentiated vertex on the left might be one with differentiated vertex on the right (vice versa)
        // for further evaluation as part of a bigger diagram they need to be reordered to recover the correct dGammaL and dGammaR
        // acc. to symmetry relations (enforce_symmetry() assumes full symmetry)
        dGammaL_half1.half1().reorder_due2antisymmetry(dGammaR_half1.half1());
#endif


        /// save intermediate states:
        if (save_intermediate) {
            int i = 2;
            State<Q> dPsi_L(dGammaL_half1, dPsi.selfenergy);
            State<Q> dPsi_R(dGammaR_half1, dPsi.selfenergy);
            State<Q> dPsi_T(dGammaT, dPsi.selfenergy);
            if (iteration == 0) {
                write_state_to_hdf<Q>(dir_str + "dPsi_L" + "_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(i),
                             Psi.Lambda, nLambda_layers, dPsi_L);
                write_state_to_hdf<Q>(dir_str + "dPsi_R" + "_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(i),
                             Psi.Lambda, nLambda_layers, dPsi_R);
                write_state_to_hdf<Q>(dir_str + "dPsi_T" + "_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(i),
                             Psi.Lambda, nLambda_layers, dPsi_T);
#ifdef DEBUG_SYMMETRIES
                // dPsi_L.vertex.check_symmetries("dPsi_L"); // we don't expect these to have the full symmetry of the vertex
                        // dPsi_R.vertex.check_symmetries("dPsi_R"); // we don't expect these to have the full symmetry of the vertex
                        dPsi_T.vertex.check_symmetries("dPsi_T");
#endif
            } else {
                add_state_to_hdf<Q>(dir_str + "dPsi_L" + "_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(i), iteration, dPsi_L, false);
                add_state_to_hdf<Q>(dir_str + "dPsi_R" + "_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(i), iteration, dPsi_R, false);
                add_state_to_hdf<Q>(dir_str + "dPsi_T" + "_RKstep" + std::to_string(rkStep) + "_forLoop" + std::to_string(i), iteration, dPsi_T, false);

            }
        }

        if (N_LOOPS >= 3) {

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
            // initialize central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
            Vertex<Q> dGammaC_tbar(Lambda);
            dGammaC_tbar.set_frequency_grid(Psi.vertex);
            dGammaC_tbar.set_Ir(true);
#endif

            for (int i = 3; i <= N_LOOPS; i++) {


                // create non-symmetric vertex with differentiated vertex on the left (full dGammaL, containing half 1 and 2)
                // assign half 1 to dGammaL
                // assign half 2 as half 1 of dGammaR [symmetric -> half1()=half2()]
                GeneralVertex<Q, non_symmetric> dGammaL(dGammaL_half1.half1(), dGammaR_half1.half1());
                if (VERBOSE) {
                    compare_with_FDTs(dGammaL, Lambda, iteration, "dGammaL_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(2), false, nLambda_layers);
                }


                // insert this non-symmetric vertex on the right of the bubble
                if (VERBOSE) print("Compute dGammaC (right insertion) ( ", i,"-loop): \n");
                Vertex<Q> dGammaC_r = calculate_dGammaC_right_insertion(Psi.vertex, dGammaL, Pi);
                if(VERBOSE) dGammaC_r.half1().check_vertex_resolution();
                if (VERBOSE) {
                    compare_with_FDTs(dGammaC_r, Lambda, iteration, "dGammaC_r_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(2), false, nLambda_layers);
                }

                // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
                // assign half 1
                // assign half 2 as half 1 of dGammaL
                GeneralVertex<Q, non_symmetric> dGammaR (dGammaR_half1.half1(), dGammaL_half1.half1());
                if (VERBOSE) {
                    compare_with_FDTs(dGammaR, Lambda, iteration, "dGammaR_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(2), false, nLambda_layers);
                }


                // insert this non-symmetric vertex on the left of the bubble
                Vertex<Q> dGammaC_l = calculate_dGammaC_left_insertion(dGammaR, Psi.vertex, Pi);
                if (VERBOSE) {
                    compare_with_FDTs(dGammaC_l, Lambda, iteration, "dGammaC_l_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(2), false, nLambda_layers);
                }

                // symmetrize by averaging left and right insertion
                Vertex<Q> dGammaC = (dGammaC_r + dGammaC_l) * 0.5;                  /// TODO: Find better solution --> K2 in dGammaC_l is bad, K3 in dGammaC_l can be obtained from dGammaC_r

                /// save intermediate states:
                if (save_intermediate) {
                    State<Q> dPsi_C (dGammaC, dPsi.selfenergy);
                    State<Q> dPsi_C_left (dGammaC_l, dPsi.selfenergy);
                    State<Q> dPsi_C_right(dGammaC_r, dPsi.selfenergy);
                    if (iteration == 0) {
                        write_state_to_hdf<Q>(dir_str+"dPsi_C"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i),       Psi.Lambda, nLambda_layers, dPsi_C);
                        write_state_to_hdf<Q>(dir_str+"dPsi_C_left"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i),  Psi.Lambda, nLambda_layers, dPsi_C_left);
                        write_state_to_hdf<Q>(dir_str+"dPsi_C_right"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, nLambda_layers, dPsi_C_right);
#ifdef DEBUG_SYMMETRIES
                        //dPsi_C_left.vertex.check_symmetries("dPsi_C_left");   // we don't expect these to have the full symmetry of the vertex
                        //dPsi_C_right.vertex.check_symmetries("dPsi_C_right"); // we don't expect these to have the full symmetry of the vertex
                        dPsi_C.vertex.check_symmetries("dPsi_C");
#endif
                    }
                    else {
                        add_state_to_hdf<Q>(dir_str+"dPsi_C"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i),       iteration, dPsi_C,       false);
                        add_state_to_hdf<Q>(dir_str+"dPsi_C_left"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i),  iteration, dPsi_C_left,  false);
                        add_state_to_hdf<Q>(dir_str+"dPsi_C_right"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), iteration, dPsi_C_right, false);

                    }
                    if (VERBOSE) {
                        print("Analyze tails of dPsi_C");
                        dPsi_C.analyze_tails();
                    }
                }


                if (VERBOSE) print("Compute dGammaL ( ", i,"-loop): \n");
                dGammaL_half1 = calculate_dGammaL(dGammaT, Psi.vertex, Pi);
                if(VERBOSE) dGammaL_half1.half1().check_vertex_resolution();
                if (VERBOSE) print("Compute dGammaR ( ", i,"-loop): \n");
                dGammaR_half1 = calculate_dGammaR(dGammaT, Psi.vertex, Pi);
                if(VERBOSE) dGammaR_half1.half1().check_vertex_resolution();

                dGammaT = dGammaL_half1 + dGammaC +
                          dGammaR_half1; // since sum dGammaL + dGammaR is symmetric, half 1 is sufficient
                dPsi.vertex += dGammaT;
#ifndef DEBUG_SYMMETRIES
                // subdiagrams don't fulfill the full symmetry of the vertex
                // the symmetry-related diagram with a differentiated vertex on the left might be one with differentiated vertex on the right (vice versa)
                // for further evaluation as part of a bigger diagram they need to be reordered to recover the correct dGammaL and dGammaR
                // acc. to symmetry relations (enforce_symmetry() assumes full symmetry)
                dGammaL_half1.half1().reorder_due2antisymmetry(dGammaR_half1.half1());
#endif
#ifdef SELF_ENERGY_FLOW_CORRECTIONS
                // extract central part of the vertex flow in the a and p channels (\bar{t}), needed for self-energy corrections
                //Vertex<Q> dGammaC_ap(Lambda);                   // initialize new vertex
                //dGammaC_ap.set_frequency_grid(Psi.vertex);
                //dGammaC_ap.avertex() = dGammaC.avertex();  // copy results from calculations above
                //dGammaC_ap.pvertex() = dGammaC.pvertex();
                dGammaC_tbar += dGammaC;                      // add the i-loop contribution to the full dGammaC_tbar
#endif
                //if(vertexConvergedInLoops(dGammaT, dPsi.vertex))
                //    break;


                /// save intermediate states:
                if (save_intermediate) {
                    State<Q> dPsi_L(dGammaL_half1, dPsi.selfenergy);
                    State<Q> dPsi_R(dGammaR_half1, dPsi.selfenergy);
                    State<Q> dPsi_T(dGammaT, dPsi.selfenergy);
                    if (iteration == 0) {
                        write_state_to_hdf<Q>(dir_str+"dPsi_L"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, nLambda_layers, dPsi_L);
                        write_state_to_hdf<Q>(dir_str+"dPsi_R"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, nLambda_layers, dPsi_R);
                        write_state_to_hdf<Q>(dir_str+"dPsi_T"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), Psi.Lambda, nLambda_layers, dPsi_T);
#ifdef DEBUG_SYMMETRIES
                        // dPsi_L.vertex.check_symmetries("dPsi_L"); // we don't expect these to have the full symmetry of the vertex
                        // dPsi_R.vertex.check_symmetries("dPsi_R"); // we don't expect these to have the full symmetry of the vertex
                        dPsi_T.vertex.check_symmetries("dPsi_T");
#endif
                    }
                    else {
                        add_state_to_hdf<Q>(dir_str+"dPsi_L"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), iteration, dPsi_L, false);
                        add_state_to_hdf<Q>(dir_str+"dPsi_R"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), iteration, dPsi_R, false);
                        add_state_to_hdf<Q>(dir_str+"dPsi_T"+"_RKstep"+std::to_string(rkStep)+"_forLoop"+std::to_string(i), iteration, dPsi_T, false);

                    }
                    if (VERBOSE) {
                        print("Analyze tails of dPsi_L");
                        dPsi_L.analyze_tails();
                        print("Analyze tails of dPsi_R");
                        dPsi_R.analyze_tails();
                        print("Analyze tails of dPsi_T");
                        dPsi_T.analyze_tails();
                    }
                }

            }

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
            // compute multiloop corrections to self-energy flow
            State<Q> Psi_SEcorrection(dPsi, Lambda);
            selfEnergyFlowCorrections(Psi_SEcorrection.selfenergy, dGammaC_tbar, Psi, G); // isolated SE correction
            dPsi.selfenergy += Psi_SEcorrection.selfenergy;

/// save intermediate states:
            if (save_intermediate) {
                State<Q> dPsi_C_tbar(dGammaC_tbar, dPsi.selfenergy);
                if (iteration == 0) {
                    write_state_to_hdf<Q>(dir_str+"dPsi_C_tbar_RKstep"+std::to_string(rkStep), Psi.Lambda, nODE + U_NRG.size() + 1, dPsi_C_tbar);

//#ifdef DEBUG_SYMMETRIES
                    write_state_to_hdf<Q>(dir_str+"SE_correction_RKstep"+std::to_string(rkStep), Psi.Lambda, nODE + U_NRG.size() + 1, Psi_SEcorrection);
//#endif
                }
                else {
                    add_state_to_hdf<Q>(dir_str+"dPsi_C_tbar_RKstep"+std::to_string(rkStep), iteration, dPsi_C_tbar, false);
                    add_state_to_hdf<Q>(dir_str+"SE_correction_RKstep"+std::to_string(rkStep), iteration, Psi_SEcorrection, false);
#ifdef DEBUG_SYMMETRIES
                    add_state_to_hdf<Q>(dir_str+"SE_correction_RKstep"+std::to_string(rkStep), iteration, Psi_SEcorrection);
#endif

                }
            }
            //TODO(low): Implement self-energy iterations (see lines 37-39 of pseudo-code).
            //if(selfEnergyConverged(Psi.selfenergy, Lambda))
            //    break;
#endif

        }
    }

    if (true and mpi_world_rank() == 0) { //
        print("Time needed for evaluation of RHS --> ");
        get_time(t0); // measure time for one iteration
    }

    return dPsi;

}

template <typename Q>
class rhs_n_loop_flow_t {
public:
    std::size_t rk_step = 0;
    std::size_t iteration = 0;
            void operator() (const State<Q>& Psi, State<Q>& dState_dLambda,  const double Lambda) {
                dState_dLambda = rhs_n_loop_flow(Psi, Lambda, vec<size_t>({iteration, rk_step}));
                rk_step++;
            }
};

template <typename Q>
void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator<Q>& S){
    // Self-energy flow
    loop(dPsiSelfEnergy, PsiVertex, S, true); // Loop for the Self-Energy calculation
}

template <typename Q, class Bubble_Object>
void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& dPi){
    // Vertex flow
    for (char r: "apt") {
        bubble_function(dPsiVertex, PsiVertex, PsiVertex, dPi, r);  // Differentiated bubble in channel r \in {a, p, t}
    }
}



template <typename Q, class Bubble_Object>
auto calculate_dGammaL(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>{
    Vertex<Q> dGammaL(Lambda_ini);
    dGammaL.set_frequency_grid(PsiVertex);

    Vertex<Q> dPsiVertex_calc = dPsiVertex;
    dPsiVertex_calc.set_Ir(true); // Only use the r-irreducible part

    for (char r: "apt") {
        bubble_function(dGammaL, dPsiVertex_calc, PsiVertex, Pi, r);
    }

    return dGammaL;
}
template <typename Q, class Bubble_Object>
auto calculate_dGammaR(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Bubble_Object& Pi) -> Vertex<Q>{
    Vertex<Q> dGammaR(Lambda_ini);
    dGammaR.set_frequency_grid(PsiVertex);

    Vertex<Q> dPsiVertex_calc = dPsiVertex;
    dPsiVertex_calc.set_Ir(true); // Only use the r-irreducible part

    for (char r: "apt") {
        bubble_function(dGammaR, PsiVertex, dPsiVertex_calc, Pi, r);
    }

    return dGammaR;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_right_insertion(const Vertex<Q>& PsiVertex, const GeneralVertex<Q, non_symmetric>& nonsymVertex,
                                       const Bubble_Object& Pi) -> Vertex<Q> {
    Vertex<Q> dGammaC (Lambda_ini);
    dGammaC.set_frequency_grid(PsiVertex);

    GeneralVertex<Q, non_symmetric> nonsymVertex_calc = nonsymVertex;
    nonsymVertex_calc.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    for (char r: "apt") {
        bubble_function(dGammaC, PsiVertex, nonsymVertex_calc, Pi, r);
    }

    return dGammaC;
}

template <typename Q, class Bubble_Object>
auto calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric>& nonsymVertex, const Vertex<Q>& PsiVertex,
                                      const Bubble_Object& Pi) -> Vertex<Q> {
    Vertex<Q> dGammaC (Lambda_ini);
    dGammaC.set_frequency_grid(PsiVertex);

    GeneralVertex<Q, non_symmetric> nonsymVertex_calc = nonsymVertex;
    nonsymVertex_calc.set_only_same_channel(true); // only use channel r of this vertex when computing r bubble

    for (char r: "apt") {
        bubble_function(dGammaC, nonsymVertex_calc, PsiVertex, Pi, r);
    }

    return dGammaC;
}

// compute multiloop corrections to self-energy flow
template <typename Q>
void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G){
    // TODO(low): also implement self-energy flow via differentiated SDE
    // TODO(low): iterate self-energy corrections (feedback of full SE flow into vertex flow etc.)?

    SelfEnergy<Q> dSigma_tbar(Psi.selfenergy.frequencies);
    SelfEnergy<Q> dSigma_t(Psi.selfenergy.frequencies);

    // compute first multiloop correction to self-energy flow, irreducible in the t channel
    loop(dSigma_tbar, dGammaC_tbar, G, true);

    // compute second multiloop correction to self-energy flow, reducible in the t channel
    Propagator<Q> extension (G.Lambda, Psi.selfenergy, dSigma_tbar, 'e');
    loop(dSigma_t, Psi.vertex, extension, true);

    dPsiSelfEnergy = dSigma_tbar + dSigma_t;

}

template <typename Q>
auto vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma) -> bool {
    return (dGamma_T.norm() / dGamma.norm() < converged_tol);
}
template <typename Q>
auto selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator<Q>& dG) -> bool {
    Propagator<Q> compare(dG.Lambda, PsiSelfEnergy, dPsiSelfEnergy, 'k');
    compare += dG*((Q)-1.);

    return (  compare.norm()/ dG.norm() < converged_tol );
}

#endif //RIGHT_HAND_SIDES_H
