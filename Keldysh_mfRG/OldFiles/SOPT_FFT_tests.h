#ifndef FPP_MFRG_SOPT_FFT_TESTS_H
#define FPP_MFRG_SOPT_FFT_TESTS_H

#include "../parameters.h"
#include "../data_structures.h"
#include "../selfenergy.h"
#include "../propagator.h"
#include "fourier_trafo.h"
#include "../ODE_solvers.h"


// Temporary vectors bfreqs, ffreqs, used in right_hand_sides.h, fourier_trafo.h, testFunctions.h, integrator.h
FrequencyGrid frequencyGrid_bos ('b', 1, Lambda_ini);
FrequencyGrid frequencyGrid_fer ('f', 1, Lambda_ini);
rvec bfreqs = frequencyGrid_bos.w;
rvec ffreqs = frequencyGrid_fer.w;

/// ------ TEST FUNCTIONS ------ ///

cvec dSOPT_FFT_K1a_rhs(const cvec& K1a, const double Lambda) { // return differentiated K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin (Lambda); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<comp> G0(Lambda, SEin,'g'); // bare propagator
    Propagator<comp> S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dK1a_1(nw1_a); // output cvec: bare differentiated K1a, component 1
    diffSOPT_FFT_K1a_R(dK1a_1, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dK1a_1;
}

cvec SOPT_FFT_K1a_rhs(const double Lambda) { // return (Lambda-dependent) K1a_1 using SOPT_FFT for testing
    SelfEnergy<comp> SEin (Lambda); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<comp> G0(Lambda, SEin,'g'); // bare propagator
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
    SelfEnergy<comp> SEin (Lambda); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<comp> G0(Lambda, SEin,'g'); // bare propagator
    Propagator<comp> S0(Lambda, SEin,'s'); // bare differentiated propagator
    cvec dSEout(nSE); // output cvec: bare differentiated SE, retarded component
    diffSOPT_FFT_SE_R(dSEout, G0, S0, glb_U, 10000, 80.); // fill ouput cvec
    return dSEout;
}

cvec SOPT_FFT_SE_rhs(const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    SelfEnergy<comp> SEin (Lambda); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<comp> G0(Lambda, SEin,'g'); // bare propagator
    cvec SEout(nSE); // output cvec: bare SE, retarded component
    SOPT_FFT_SE_R(SEout, G0, glb_U, 10000, 80.); // fill ouput cvec
    return SEout;
}

SelfEnergy<comp> SOPT_FFT_SE_rhs(const SelfEnergy<comp>& SEin, const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    Propagator<comp> G(Lambda, SEin,'g'); // full propagator
    SelfEnergy<comp> SEout (Lambda); // output self-energy
    SOPT_FFT(SEout, G, glb_U, 100000, 500., false); // fill ouput self-energy
    return SEout;
}

SelfEnergy<comp> SOPTlad_FFT_SE_rhs(const SelfEnergy<comp>& SEin, const double Lambda) { // return (Lambda-dependent) SE using SOPT_FFT for testing
    Propagator<comp> G(Lambda, SEin,'g'); // full propagator
    SelfEnergy<comp> SEout (Lambda); // output self-energy
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
    SelfEnergy<comp> SEini (Lambda_fin); // trivial self-energy
    SEini.initialize(glb_U/2., 0.); // initialize with Hartree term
    SelfEnergy<comp> SEfin (Lambda_fin); // final self-energy
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
    SelfEnergy<comp> SEin (Lambda); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<comp> S0(Lambda, SEin,'s'); // bare differentiated propagator from class
    cvec dG(nSE); // bare differentiated propagator as cvec
    for (int i=0; i<nSE; ++i) // fill output cvec
        dG[i] = S0.valsmooth(0, ffreqs[i], 0);
    return dG;
}

void test_ODE_G(const int N_ODE) { // test ODE applied to bare (differentiated) propagator
    bool write_flag = false; // whether to write result in hdf5 file
    cvec G_dir(nw1_a), G_fin(nw1_a), G_ini(nw1_a); // direct, final, initial propagator as cvec
    double Lambda_fin = 1.; double Lambda_ini = 2.; // end points of flow
    SelfEnergy<comp> SEin (Lambda_ini); // trivial self-energy
    SEin.initialize(glb_U/2., 0.); // initialize with Hartree term
    Propagator<comp> G_ini_prop(Lambda_ini, SEin, 'g'); // initial propagator from class
    Propagator<comp> G_dir_prop(Lambda_fin, SEin, 'g'); // final propagator from class
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
 * Function that prints out a .h5 file with the value of the rhs of a SOPT flow at the given Lambda for both an FFT and an fRG calculation
 * @param Lambda    : Lambda at which the derivatives are to be calculated
 */
template <typename Q>
void test_derivatives_SE(double Lambda){
    cvec blah(nw1_a);
    State<Q> sopt (Lambda);
    sopt.initialize();
    cvec rhs_SOPT_FFT_K1a = dSOPT_FFT_SE_rhs(blah, Lambda);

    sopt_state(sopt, Lambda);
    State<Q> rhs_flow = rhs_state_flow_SOPT_0(sopt, Lambda);

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


#endif //FPP_MFRG_SOPT_FFT_TESTS_H
