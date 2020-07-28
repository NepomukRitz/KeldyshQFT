#ifndef RIGHT_HAND_SIDES_H
#define RIGHT_HAND_SIDES_H

#include "data_structures.h"        // real/complex vector classes, imag. unit
#include "write_data2file.h"        // writing data into text or hdf5 files
#include "propagator.h"             // propagator to perform second-order perturbation theory (SOPT)
#include "selfenergy.h"             // self-energy used in SOPT
#include "state.h"                  // state to perform full flow
#include "loop.h"                   // compute self-energy loop
#include "bubbles.h"                // compute vertex bubbles
#include "parameters.h"             // system parameters (lengths of vectors etc.)
#include "fourier_trafo.h"          // SOPT from Fast Fourier transform (FFT)
#include "solvers.h"                // ODE solver
#include "assert.h"

using namespace std;

template <typename Q> auto rhs_n_loop_flow(const State<Q>& Psi, double Lambda) -> State<Q>;
template <typename Q> void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator& S);
template <typename Q> void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& Psi, const Propagator& G, const Propagator& dG);

template <typename Q> void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammatbar_C, const State<Q>& Psi, const Propagator& G);

template <typename Q> auto calculate_dGammaL(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator& G) -> Vertex<Q>;
template <typename Q> auto calculate_dGammaR(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator& G) -> Vertex<Q>;
template <typename Q> auto calculate_dGammaC_ap(const Vertex<Q>& PsiVertex, const Vertex<Q>& dGammaL, const Propagator& G) -> Vertex<Q>;
template <typename Q> auto calculate_dGammaC_t (const Vertex<Q>& PsiVertex, const Vertex<Q>& dGammaL, const Propagator& G) -> Vertex<Q>;


template <typename Q> bool vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma);
template <typename Q> bool selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator& dG);

/// ------ TEST FUNCTIONS ------ ///
// TODO: Here are also functions that should belong to testFunctions.h. On the other hand side,
//  in testFunctions.h are functions that would rather belong here.

cvec dSOPT_FFT_K1a_rhs(const cvec& K1a, const double Lambda) { // return differentiated K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dK1a_1(nw1_a); // output cvec: bare differentiated K1a, component 1
    diffSOPT_FFT_K1a_R(dK1a_1, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dK1a_1;
}

cvec SOPT_FFT_K1a_rhs(const double Lambda) { // return (Lambda-dependent) K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    cvec K1a_1(nw1_a); // output cvec: bare K1a, component 1
    SOPT_FFT_K1a_R(K1a_1, G0, glb_U, 10000, 80.); // fill ouput cvec
    return K1a_1;
}

void test_ODE_SOPT_FFT_K1a(const int N_ODE) { // test ODE solver using K1a from SOPT_FFT
    bool write_flag = false; // whether to write output in hdf5
    cvec K1a_dir(nw1_a), K1a_fin(nw1_a), K1a_ini(nw1_a); // direct, final, initial K1a_1
    //double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow -> defined in parameters.h
    //SelfEnergy<comp> SEin; // trivial self-energy
    //SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    //Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    K1a_ini = SOPT_FFT_K1a_rhs(Lambda_ini); // direct calculation of initial K1a
    K1a_dir = SOPT_FFT_K1a_rhs(Lambda_fin); // direct calculation of final K1a
    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, dSOPT_FFT_K1a_rhs, N_ODE); // final K1a from ODE
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare K1a_1. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_K1a.h5",
                                  {"v", "K1a_dir_R", "K1a_dir_I", "K1a_fin_R", "K1a_fin_I", "K1a_ini_R", "K1a_ini_I"},
                                  {bfreqs, K1a_dir.real(), K1a_dir.imag(), K1a_fin.real(), K1a_fin.imag(), K1a_ini.real(), K1a_ini.imag()});
}

cvec dSOPT_FFT_SE_rhs(const cvec& SE, const double Lambda) { // return differentiated SE using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dSEout(nSE); // output cvec: bare differentiated SE, retarded component
    diffSOPT_FFT_SE_R(dSEout, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dSEout;
}

cvec SOPT_FFT_SE_rhs(const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G0(Lambda, SEin,'g'); // bare propagator
    cvec SEout(nSE); // output cvec: bare SE, retarded component
    SOPT_FFT_SE_R(SEout, G0, glb_U, 10000, 80.); // fill ouput cvec
    return SEout;
}

SelfEnergy<comp> SOPT_FFT_SE_rhs(const SelfEnergy<comp>& SEin, const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    Propagator G(Lambda, SEin,'g'); // full propagator
    SelfEnergy<comp> SEout; // output self-energy
    SOPT_FFT(SEout, G, glb_U, 100000, 500., false); // fill ouput self-energy
    return SEout;
}

SelfEnergy<comp> SOPTlad_FFT_SE_rhs(const SelfEnergy<comp>& SEin, const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    Propagator G(Lambda, SEin,'g'); // full propagator
    SelfEnergy<comp> SEout; // output self-energy
    SOPT_FFT(SEout, G, glb_U, 100000, 500., true); // fill ouput self-energy
    return SEout;
}

void test_ODE_SOPT_FFT_SE(const int N_ODE) { // test ODE solver using SE from SOPT_FFT
    bool write_flag = true; // whether to write output in hdf5
    cvec SE_dir(nSE), SE_fin(nSE), SE_ini(nSE); // direct, final, initial SE
    //double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow -> defined in parameters.h
    //SelfEnergy<comp> SEin; // trivial self-energy
    //SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    //Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    SE_ini = SOPT_FFT_SE_rhs(Lambda_ini); // direct calculation of initial K1a
    SE_dir = SOPT_FFT_SE_rhs(Lambda_fin); // direct calculation of final K1a
    ODE_solver_RK4(SE_fin, Lambda_fin, SE_ini, Lambda_ini, dSOPT_FFT_SE_rhs, N_ODE); // final K1a from ODE
    cvec SE_dif = SE_dir + ( SE_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare SE_R. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << SE_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_SE.h5",
                                  {"v", "SE_dir_R", "SE_dir_I", "SE_fin_R", "SE_fin_I", "SE_ini_R", "SE_ini_I"},
                                  {ffreqs, SE_dir.real(), SE_dir.imag(), SE_fin.real(), SE_fin.imag(), SE_ini.real(), SE_ini.imag()});
}

void test_SCE_SOPT_FFT_SE(const int N_SCE_SOPT, const int N_SCE_SOPTlad) { // test SCE solver using SE from SOPT(lad)_FFT
    bool write_flag = true; // whether to write output in hdf5
    SelfEnergy<comp> SEini; // trivial self-energy
    SEini.initialize(glb_U/2., 0.); // initialize with Hartree term
    SelfEnergy<comp> SEfin; // final self-energy
    SCE_solver(SEfin, SEini, Lambda_fin, SOPT_FFT_SE_rhs, N_SCE_SOPT, 0.); // final self-energy after self-consistency iterations
    SCE_solver(SEfin, SEfin, Lambda_fin, SOPTlad_FFT_SE_rhs, N_SCE_SOPTlad, 0.); // dito including ladder
    cout << "Testing SCE for SEOPT using " << N_SCE_SOPT << " iterations for self-consistent SOPT and " << N_SCE_SOPTlad << " iterations including ladders." << endl;
    if(write_flag) {
        cvec SEfin_R(nSE), SEfin_K(nSE); // store result in cvecs for writing
        for (int iv=0; iv<nSE; ++iv) {
            SEfin_R[iv] = SEfin.val(0, iv, 0);
            SEfin_K[iv] = SEfin.val(1, iv, 0);
        }
        write_h5_rvecs("SOPT_FFT_SCE_SE.h5",
                {"v", "SE_R_R", "SE_R_I", "SE_K_R", "SE_K_I"},
                {ffreqs, SEfin_R.real(), SEfin_R.imag(), SEfin_K.real(), SEfin_K.imag()});
    }
}

cvec dG_rhs(const cvec& G, const double Lambda) { // return bare differentiated propagator for testing purposes
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator S0(Lambda, SEin,'s'); // bare differentiated propagator from class
    cvec dG(nSE); // bare differentiated propagator as cvec
    for (int i=0; i<nSE; ++i) // fill output cvec
        dG[i] = S0.valsmooth(0, ffreqs[i], 0);
    return dG;
}

void test_ODE_G(const int N_ODE) { // test ODE applied to bare (differentiated) propagator
    bool write_flag = false; // whether to write result in hdf5 file
    cvec G_dir(nw1_a), G_fin(nw1_a), G_ini(nw1_a); // direct, final, initial propagator as cvec
    double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow
    SelfEnergy<comp> SEin; // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator G_ini_prop(Lambda_ini, SEin, 'g'); // initial propagator from class
    Propagator G_dir_prop(Lambda_fin, SEin, 'g'); // final propagator from class
    for (int i=0; i<nSE; ++i) { // fill propagator cvecs
        G_ini[i] = G_ini_prop.valsmooth(0, ffreqs[i], 0);
        G_dir[i] = G_dir_prop.valsmooth(0, ffreqs[i], 0);
    }
    ODE_solver_RK4(G_fin, Lambda_fin, G_ini, Lambda_ini, dG_rhs, N_ODE); // final G from ODE flow
    cvec G_dif = G_dir + ( G_fin*(-1.) ); // difference in results
    cout << "Testing ODE for bare proapgator. Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << G_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("ODE_G.h5",
                                  {"v", "G_dir_R", "G_dir_I", "G_fin_R", "G_fin_I", "G_ini_R", "G_ini_I"},
                                  {ffreqs, G_dir.real(), G_dir.imag(), G_fin.real(), G_fin.imag(), G_ini.real(), G_ini.imag()});
}


/**
 * Function to implement an n-loop flow (without Katanin substitution).
 * @param Psi   : Known state of the State at Lambda
 * @param Lambda:  Scale at which the calculation is being performed
 * @return dPsi : The derivative at Lambda, which includes the differential vertex as well as self-energy at scale Lambda
 */
template <typename Q>
auto rhs_n_loop_flow(const State<Q>& Psi, const double Lambda) -> State<Q>{

    static_assert(N_LOOPS>=1, "");

    State<Q> dPsi; // result

    Propagator S (Lambda, Psi.selfenergy, 's');
    Propagator G (Lambda, Psi.selfenergy, 'g');

    //For flow without self-energy, comment out this line
    selfEnergyOneLoopFlow(dPsi.selfenergy, Psi.vertex, S);

    //Propagator dG (Lambda, Psi.selfenergy, dPsi.selfenergy, 'k');
    //Run alternatively, for no self-energy feedback
    Propagator dG (Lambda, Psi.selfenergy, 's');

    vertexOneLoopFlow(dPsi.vertex, Psi.vertex, G, dG);

#if N_LOOPS>=2

    Vertex<Q> dGammaL = calculate_dGammaL(dPsi.vertex, Psi.vertex, G);
    Vertex<Q> dGammaR = calculate_dGammaR(dPsi.vertex, Psi.vertex, G);
    Vertex<Q> dGammaT = dGammaL + dGammaR;
    dPsi.vertex += dGammaT;
#endif

#if N_LOOPS>=3

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
    Vertex<Q> dGammatbar_C(n_spin);
#endif

    for(int i=3; i<N_LOOPS; i++){
        Vertex<Q> dGammaC_ap = calculate_dGammaC_ap(Psi.vertex, dGammaL, G);
        Vertex<Q> dGammaC_t  = calculate_dGammaC_t (Psi.vertex, dGammaL, G);
        Vertex<Q> dGammaC = dGammaC_ap + dGammaC_t;

        dGammaL = calculate_dGammaL(dGammaT, Psi.vertex, G);
        dGammaR = calculate_dGammaR(dGammaT, Psi.vertex, G);

        dGammaT = dGammaL + dGammaC + dGammaR;
        dPsi.vertex += dGammaT;
#ifdef SELF_ENERGY_FLOW_CORRECTIONS
        dGammatbar_C += dGammaC_ap;
#endif
        //if(vertexConvergedInLoops(dGammaT, dPsi.vertex))
        //    break;
    }

#ifdef SELF_ENERGY_FLOW_CORRECTIONS
    selfEnergyFlowCorrections(dPsi.selfenergy, dGammatbar_C, Psi, G);

    //TODO These are supposed to be lines 37-39 of pseudo-code. What do these refer to?
    //if(selfEnergyConverged(Psi.selfenergy, Lambda))
    //    break;
#endif

#endif

    return dPsi;

}

template <typename Q>
void selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& PsiVertex, const Propagator& S){
    // Self-energy flow
    loop(dPsiSelfEnergy, PsiVertex, S, true); // Loop for the Self-Energy calculation
}

template <typename Q>
void vertexOneLoopFlow(Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator& G, const Propagator& dG){

    // Vertex flow
    bubble_function(dPsiVertex, PsiVertex, PsiVertex, G, dG, 'a', true, '.');  // Differentiated bubble in the a-channel
    bubble_function(dPsiVertex, PsiVertex, PsiVertex, G, dG, 'p', true, '.');  // Differentiated bubble in the p-channel
    bubble_function(dPsiVertex, PsiVertex, PsiVertex, G, dG, 't', true, '.');  // Differentiated bubble in the t-channel
}

template <typename Q>
auto calculate_dGammaL(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator& G) -> Vertex<Q>{
    Vertex<Q> dGammaL(n_spin);

    bubble_function(dGammaL, dPsiVertex, PsiVertex, G, G, 'a', false, 'L');
    bubble_function(dGammaL, dPsiVertex, PsiVertex, G, G, 'p', false, 'L');
    bubble_function(dGammaL, dPsiVertex, PsiVertex, G, G, 't', false, 'L');

    return dGammaL;
}
template <typename Q>
auto calculate_dGammaR(const Vertex<Q>& dPsiVertex, const Vertex<Q>& PsiVertex, const Propagator& G) -> Vertex<Q>{
    Vertex<Q> dGammaR(n_spin);

    bubble_function(dGammaR, PsiVertex, dPsiVertex, G, G, 'a', false, 'R');
    bubble_function(dGammaR, PsiVertex, dPsiVertex, G, G, 'p', false, 'R');
    bubble_function(dGammaR, PsiVertex, dPsiVertex, G, G, 't', false, 'R');

    return dGammaR;
}

template <typename Q>
auto calculate_dGammaC_ap(const Vertex<Q>& PsiVertex, const Vertex<Q>& dGammaL, const Propagator& G) -> Vertex<Q>{
    Vertex<Q> dGamma_C_ap(n_spin);

    bubble_function(dGamma_C_ap, PsiVertex, dGammaL, G, G, 'a', false, 'L');
    bubble_function(dGamma_C_ap, PsiVertex, dGammaL, G, G, 'p', false, 'L');

    return dGamma_C_ap;
}

template <typename Q>
auto calculate_dGammaC_t (const Vertex<Q>& PsiVertex, const Vertex<Q>& dGammaL, const Propagator& G) -> Vertex<Q>{
    Vertex<Q> dGamma_C_t(n_spin);
    bubble_function(dGamma_C_t, PsiVertex, dGammaL, G, G, 't', false, 'L');
    return dGamma_C_t;
}

template <typename Q>
void selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q>& dGammatbar_C, const State<Q>& Psi, const Propagator& G){

    SelfEnergy<Q> dSigma_tbar;
    SelfEnergy<Q> dSigma_t;

    loop(dSigma_tbar, dGammatbar_C, G, false);

    Propagator extension (G.Lambda, Psi.selfenergy, dSigma_tbar, 'e');
    loop(dSigma_t, Psi.vertex, extension, false);

    dPsiSelfEnergy += dSigma_tbar + dSigma_t;

}

template <typename Q>
auto vertexConvergedInLoops(Vertex<Q>& dGamma_T, Vertex<Q>&dGamma) -> bool {
    return (dGamma_T.norm() / dGamma.norm() < converged_tol);
}
template <typename Q>
auto selfEnergyConverged(SelfEnergy<Q>& dPsiSelfEnergy, SelfEnergy<Q>& PsiSelfEnergy, Propagator& dG) -> bool {
    Propagator compare(dG.Lambda, PsiSelfEnergy, dPsiSelfEnergy, 'k');
    compare += dG*((Q)-1.);

    return (  compare.norm()/ dG.norm() < converged_tol );
}

#endif //RIGHT_HAND_SIDES_H
