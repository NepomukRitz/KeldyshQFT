#ifndef KELDYSH_MFRG_RIGHT_HAND_SIDES_H
#define KELDYSH_MFRG_RIGHT_HAND_SIDES_H

#include "propagator.h"
#include "selfenergy.h"
#include "parameters.h"
#include "fourier_trafo.h"
#include "solvers.h"

using namespace std;

cvec dSOPT_K1a_rhs(const cvec& PiaEO, const double Lambda) {
    SelfEnergy<comp> SEin;
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G0(Lambda, SEin,'g');
    Propagator S0(Lambda, SEin,'s');
    cvec dPiaEO(nw1_wa);
    diffSOPTbare_FFT_SelfEnergy(dPiaEO, G0, S0, glb_U, 10000, 80.);
    return dPiaEO;
}

cvec SOPT_K1a_rhs(const double Lambda) {
    SelfEnergy<comp> SEin;
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G0(Lambda, SEin,'g');
    cvec PiaEO(nw1_wa);
    SOPTbare_FFT_SelfEnergy(PiaEO, G0, glb_U, 10000, 80.);
    return PiaEO;
}

void test_ODE_SOPT_K1a(const int N_ODE) {
    bool write_flag = false;
    cvec K1a_dir(nw1_wa), K1a_fin(nw1_wa), K1a_ini(nw1_wa); // direct, final, initial
    double Lambda_fin = 2.; double Lambda_ini = 1.;
    SelfEnergy<comp> SEin; // trivial self-energy
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G0ini(Lambda_ini, SEin, 'g'); // initial propagator
    SOPTbare_FFT_SelfEnergy(K1a_ini, G0ini, glb_U, 10000, 80.); // initial K1a
    K1a_dir = SOPT_K1a_rhs(Lambda_fin); // direct calculation of K1a
    ODE_solver_RK4(K1a_fin, Lambda_fin, K1a_ini, Lambda_ini, dSOPT_K1a_rhs, N_ODE); // final K1a
    cvec K1a_dif = K1a_dir + ( K1a_fin*(-1.) );
    cout << "Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << K1a_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("SOPT_FFT_ODE_K1a.h5",
                                  {"v", "K1a_dir_R", "K1a_dir_I", "K1a_fin_R", "K1a_fin_I", "K1a_ini_R", "K1a_ini_I"},
                                  {bfreqs, K1a_dir.real(), K1a_dir.imag(), K1a_fin.real(), K1a_fin.imag(), K1a_ini.real(), K1a_ini.imag()});
}

cvec dG_rhs(const cvec& G, const double Lambda) {
    SelfEnergy<comp> SEin;
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator S0(Lambda, SEin,'s');
    cvec dG(nSE);
    for (int i=0; i<nSE; ++i)
        dG[i] = S0.valsmooth(0, ffreqs[i]);
    return dG;
}

void test_ODE_G(const int N_ODE) {
    bool write_flag = false;
    cvec G_dir(nw1_wa), G_fin(nw1_wa), G_ini(nw1_wa); // direct, final, initial
    double Lambda_fin = 2.; double Lambda_ini = 1.;
    SelfEnergy<comp> SEin; // trivial self-energy
    for(int i=0; i<nSE; ++i) {
        SEin.setself(0, i, glb_U/2.);
        SEin.setself(1, i, 0.);
    }
    Propagator G_ini_prop(Lambda_ini, SEin, 'g'); // initial propagator
    Propagator G_dir_prop(Lambda_fin, SEin, 'g'); // initial propagator
    for (int i=0; i<nSE; ++i) {
        G_ini[i] = G_ini_prop.valsmooth(0, ffreqs[i]);
        G_dir[i] = G_dir_prop.valsmooth(0, ffreqs[i]);
    }
    ODE_solver_RK4(G_fin, Lambda_fin, G_ini, Lambda_ini, dG_rhs, N_ODE); // final G
    cvec G_dif = G_dir + ( G_fin*(-1.) );
    cout << "Using " << N_ODE << " ODE steps, the maximal difference between direct and ODE-final result is " << G_dif.max_norm() << "." << endl;
    if(write_flag) write_h5_rvecs("ODE_G.h5",
                                  {"v", "G_dir_R", "G_dir_I", "G_fin_R", "G_fin_I", "G_ini_R", "G_ini_I"},
                                  {ffreqs, G_dir.real(), G_dir.imag(), G_fin.real(), G_fin.imag(), G_ini.real(), G_ini.imag()});
}

#endif //KELDYSH_MFRG_RIGHT_HAND_SIDES_H
