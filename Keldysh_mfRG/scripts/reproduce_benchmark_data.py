import numpy as np
import os

### Run this script from the top project directory
### $ python3 ./scripts/reproduce_benchmerk_data.py

def backup_files():
    filename_main = "main.cpp"
    filename_masterparam = "parameters/master_parameters.hpp"
    filename_frequeparam = "parameters/frequency_parameters.hpp"
    filename_techniparam = "parameters/technical_parameters.hpp"
    filenames = [filename_main, filename_masterparam, filename_frequeparam, filename_techniparam]
    texts = np.array([])
    for i in range(len(filenames)):
        with open(filenames[i], 'r') as f:
            texts = np.append(texts, f.read())
    return texts

def restore_backup_files(texts):
    filename_main = "./main.cpp"
    filename_masterparam = "./parameters/master_parameters.hpp"
    filename_frequeparam = "./parameters/frequency_parameters.hpp"
    filename_techniparam = "./parameters/technical_parameters.hpp"
    filenames = [filename_main, filename_masterparam, filename_frequeparam, filename_techniparam]
    for i in range(len(filenames)):
        with open(filenames[i], 'w') as f:
            f.write(texts[i])

def write_main():
    text = """
//////////////// AUTOMATICALLY GENERATED ////////////////
    #include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "tests/test_perturbation_theory.hpp"
#include "tests/reproduce_benchmark_data.hpp"
#include "utilities/util.hpp"
#include "utilities/hdf5_routines.hpp"
#include "tests/integrand_tests/saveIntegrand.hpp"
#include "tests/test_symmetries.hpp"
#include "perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#include "tests/test_ODE.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif


auto main(int argc, char * argv[]) -> int {
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif


    /// to reproduce a nice set of benchmark data for fast spotting of bugs, go through the following steps:
    /// 0.: choose standard settings for
    ///     MF / KF
    ///     PHS / non-PHS
    ///     T=0 / T>0
    ///     SBE / asymptotic decomposition
    ///     REG 2 / 4
    ///     fixed: DEBUG_SYMMETRIES 0, MAX_DIAG_CLASS 3, CONTOUR_BASIS 0, VECTORIZED_INTEGRATION 1, KATANIN, frequency parameters, Lambda, Gamma, U, T
    /// 1.: compute State in SOPT
    /// 2.: feed into parquet solver version 2 and do 2 parquet iterations
    /// 3.: compute one-shot result of mfRG equations with 4 loops including selfenergy corrections with Fabians' NJP formula
    reproduce_benchmark_data();



#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}

    """
    filename = "./main.cpp"
    with open(filename, 'w') as f:
        f.write(text)


def write_frequency_parameters():
    text = """
#ifndef FPP_MFRG_FREQUENCY_PARAMETERS_H
#define FPP_MFRG_FREQUENCY_PARAMETERS_H
//////////////// AUTOMATICALLY GENERATED ////////////////

/// Frequency grid parameters ///

// Grid type
#define GRID 0      // 0: use Elias' grid;      1: use hybrid grid;     2: polar coordinates (currently not for K3 yet)

//#define ADAPTIVE_GRID   // if defined: use optimization routine; if undefined: just rescale the grid;

//#define ROTATEK2 // saves and interpolates K2 data on and rotated grid (corresponds to "fermionic" parametrization)
constexpr bool BOSONIC_PARAM_FOR_K3 = false; // saves and interpolates K3 data on and rotated grid (corresponds to "bosonic" parametrization)

constexpr bool INTERPOL2D_FOR_K3 = BOSONIC_PARAM_FOR_K3 and true;


// Limits of the frequency grid vectors for the different kinds of frequencies
// (i.e. bosonic transfer frequency and fermionic frequencies

// Number of bosonic and fermionic frequency points.
// Good production values for Keldysh: nBOS = nFER = 401, nBOS2 = nFER2 = 201, nBOS3 = nFER3 = 51
//#if KELDYSH_FORMALISM
constexpr int nBOS = 401;
constexpr int nFER = 401 - (KELDYSH_FORMALISM ? 0 : 1);
// Number of frequency points for K2 and K3 classes
constexpr int nBOS2 = 41;//nBOS;
constexpr int nFER2 = 41 - (KELDYSH_FORMALISM or ZERO_TEMP ? 0 : 1);//nFER;
constexpr int nBOS3 = 11; //nBOS;
constexpr int nFER3 = 11 - (KELDYSH_FORMALISM or ZERO_TEMP ? 0 : 1); //nFER;
//#else
const int COUNT = 4;
//constexpr int nBOS = COUNT * 64 * 2 + 1;
//constexpr int nFER = COUNT * 32 * 2;
//constexpr int nBOS2 = COUNT * 8 * 2 + 1;//nBOS;
//constexpr int nFER2 = COUNT * 4 * 2;//nFER;
//constexpr int nBOS3 = COUNT * 4 * 2 + 1; //nBOS;
//constexpr int nFER3 = COUNT * 2 * 2; //nFER;
const int POSINTRANGE = 64  * COUNT;
//#endif


const double Delta_factor_K1 = 10.;
const double Delta_factor_SE = 10.;
const double Delta_factor_K2_w = 20.;
const double Delta_factor_K2_v = 10.;
const double Delta_factor_K3_w = 10.;
const double Delta_factor_K3_v = 10.;



// Number of frequency points for the self energy and the susceptibility
constexpr int nSE   = nFER;
constexpr int nPROP = nFER;

// Number of frequency points for the K1 class (bosonic freq w),
// K2 (bosonic freq w, fermionic freq v) and K3 (bosonic frequency w and fermionic freqs, both v and vp)
// for the a channel
constexpr int nw1_a  = nBOS;
constexpr int nw2_a  = nBOS2;
constexpr int nv2_a  = nFER2;
constexpr int nw3_a  = nBOS3;
constexpr int nv3_a  = nFER3;
// for the p channel
constexpr int nw1_p  = nBOS;
constexpr int nw2_p  = nBOS2;
constexpr int nv2_p  = nFER2;
constexpr int nw3_p  = nBOS3;
constexpr int nv3_p = nFER3;
// for the t channel
constexpr int nw1_t  = nBOS;
constexpr int nw2_t  = nBOS2;
constexpr int nv2_t  = nFER2;
constexpr int nw3_t  = nBOS3;
constexpr int nv3_t = nFER3;
// for a general channel r
constexpr int nw1 = nBOS;
constexpr int nw2 = nBOS2;
constexpr int nv2 = nFER2;
constexpr int nw3 = nBOS3;
constexpr int nv3 = nFER3;


#endif //FPP_MFRG_FREQUENCY_PARAMETERS_H
    """
    filename = "./parameters/frequency_parameters.hpp"
    with open(filename, 'w') as f:
        f.write(text)

def write_technical_parameters():
    text = """
    #ifndef FPP_MFRG_TECHNICAL_PARAMETERS_H
#define FPP_MFRG_TECHNICAL_PARAMETERS_H
//////////////// AUTOMATICALLY GENERATED ////////////////
//#include "../data_structures.h"

/// Technical parameters ///

//If defined, the flow of the self_energy is symmetrized, closed above and below
    //#define SYMMETRIZED_SELF_ENERGY_FLOW

// Flag whether to use MPI, comment out following to not use MPI_FLAG
#define USE_MPI
#ifdef USE_MPI
constexpr bool MPI_FLAG = true;
#else
constexpr bool MPI_FLAG = false;
#endif

//Tolerance for closeness to grid points when interpolating
constexpr double inter_tol = 1e-5;

enum interpolMethod {linear=0, linear_on_aux=1, cubic=4};
constexpr interpolMethod INTERPOLATION = linear_on_aux;

//Tolerance for loop convergence
constexpr double converged_tol = 1e-7;

//Integrator type:
// 0: Riemann sum
// 1: Simpson
// 2: Simpson + additional points
// 3: adaptive Simpson
// 4: GSL
// 5: adaptive Gauss-Lobatto with Kronrod extension (preferred)
// 6: PAID with Clenshaw-Curtis rule
constexpr int INTEGRATOR_TYPE = 5;

//Integrator tolerance
inline double integrator_tol = 1e-5;

//Simpson integraton number of steps - 10 times the largest one out of nBOS and nFER
constexpr int nINT = 1501; //(nBOS*(nBOS>=nFER) + nFER*(nBOS<nFER));

// Debug mode allows to select specific Keldysh components contributing to loop and bubbles
//#define DEBUG_MODE


#endif //FPP_MFRG_TECHNICAL_PARAMETERS_H
    """
    filename = "./parameters/technical_parameters.hpp"
    with open(filename, 'w') as f:
        f.write(text)

def write_master_parameters(isKELDYSH, isZEROTEMP, isPHS, isSBE, REG):
    text = """
    #ifndef KELDYSH_MFRG_PARAMETERS_H
#define KELDYSH_MFRG_PARAMETERS_H

    //////////////// AUTOMATICALLY GENERATED ////////////////
#include <cmath>             // log function
#include <vector>            // standard vector for Keldysh indices
#include <string>
#include <array>

// For production: uncomment the following line to switch off assert()-functions
#define NDEBUG

#define DEBUG_SYMMETRIES 0 // 0 for false; 1 for true; used for test_symmetries() -> computes the mfRG equations once without use of symmetries

constexpr bool VERBOSE = false;

#define USE_ANDERSON_ACCELERATION 0

// Determines whether the 2D Hubbard model shall be studied instead of the SIAM
//#define HUBBARD

// Determines whether the Fermi-polaron problem shall be studied instead of the SIAM
//#define FERMI_POLARON_PROBLEM

#define ZERO_TEMP """ + str(isZEROTEMP) + """ // Determines whether to work in the T = 0 limit
// Defines the formalism (not defined: Matsubara formalism, defined: Keldysh formalism)
#define KELDYSH_FORMALISM """ + str(isKELDYSH) + """ // 0 for Matsubara; 1 for Keldysh formalism
#define CONTOUR_BASIS 0     // 0 for Keldysh basis; 1 for Contour basis
#define SWITCH_SUM_N_INTEGRAL 1    // if defined: sum over internal indices within integrand
#if KELDYSH_FORMALISM or not ZERO_TEMP
#define VECTORIZED_INTEGRATION 1  // perform integrals with vector-valued integrands ; 0 for False; 1 for True;
                                 // Keldysh: vectorizes over Keldysh indices
                                 // Matsubara finite T: vectorizes Matsubara sum
#else
#define VECTORIZED_INTEGRATION 0
#endif


// Determines whether particle-hole symmetry is assumed
#define PARTICLE_HOLE_SYMM """ + str(isPHS) + """

/// Production runs parameters ///

// Defines the number of diagrammatic classes that are relevant for a code:
// 1 for only K1, 2 for K1 and K2 and 3 for the full dependencies
#define MAX_DIAG_CLASS 3
#define SBE_DECOMPOSITION """ + str(isSBE) + """
#define USE_NEW_MFRG_EQS 1      // mfRG equations for SBE approximation. Only relevant when SBE_DECOMOSITION == 1.

#define ANALYTIC_TAILS 1

#define KATANIN
#define SELF_ENERGY_FLOW_CORRECTIONS 1
const int nmax_Selfenergy_iterations = 10;
const double tol_selfenergy_correction_abs = 1e-9;
const double tol_selfenergy_correction_rel = 1e-5;
const double loop_tol_abs = 1e-9;
const double loop_tol_rel = 1e-5;

// If defined, use static K1 inter-channel feedback as done by Severin Jakobs.
// Only makes sense for pure K1 calculations.
//#define STATIC_FEEDBACK

/// Physical parameters ///

constexpr double glb_mu = 0.0;                     // Chemical potential -- w.l.o.g. ALWAYS set to zero (for the SIAM)

constexpr double glb_V = 0.;                       // Bias voltage (glb_V == 0. in equilibrium)
constexpr bool EQUILIBRIUM = true;                 // If defined, use equilibrium FDT's for propagators
                                                   // (only sensible when glb_V = 0)
#define USE_FDT 0

/// Spin parameters ///

// Number of independent spin components. n_spin = 1 with SU(2) symmetry.
#if DEBUG_SYMMETRIES
constexpr int n_spin = 2;
#else
constexpr int n_spin = 1;
#endif
constexpr int n_spin_expanded = 2;

/// Parameters for Fermi-polaron problem ///

extern double glb_muc;
extern double glb_mud;
extern double glb_mc;
extern double glb_md;
extern double glb_ainv;

/// Parameters for internal structure ///

#ifdef HUBBARD
constexpr int n_in_K1 = glb_N_transfer;                         // Number of internal indices needed for K1-objects
constexpr int n_in_K2 = glb_N_transfer * glb_N_ff;              // Number of internal indices needed for K2-objects
constexpr int n_in_K3 = glb_N_transfer * glb_N_ff * glb_N_ff;   // Number of internal indices needed for K3-objects
constexpr int n_in = n_in_K3;                                   // maximal number of internal indices
#else
constexpr int n_in_K1 = 1;
constexpr int n_in_K2 = 1;
constexpr int n_in_K3 = 1;
constexpr int n_in = 1;
#endif

/// fRG parameters ///

// Regulator
// 1: sharp cutoff, 2: hybridization flow, 3: frequency regulator (as used in Vienna, Stuttgart, Tuebingen), including Fabian's idea for Keldysh
// 4: interaction cutoff
#define REG """ + str(REG) + """



// if the following is not defined, we flow with the parameter Lambda. The flowgrid just provides suggestions for stepsizes
// if the following is     defined, we flow with t via Lambda(t) <-- flowgrid;
#define REPARAMETRIZE_FLOWGRID


// ODE solvers:
// 1 -> basic Runge-Kutta 4; // WARNING: non-adaptive!
// 2 -> Bogacki–Shampine
// 3 -> Cash-Carp
// 4 -> Dormand-Prince
#define ODEsolver 3

// Limits of the fRG flow
#if REG == 4
constexpr double Lambda_ini = 0.;// 1e4;                // NOLINT(cert-err58-cpp)
constexpr double Lambda_fin = 1;// 1e-4;
#else
const double Lambda_ini = 19.8;//pow(10,  1) ;// 1e4;
const double Lambda_fin = 0.2 ;// 1e-4;
#endif
constexpr double Lambda_scale = 1./200.;             //Scale of the log substitution
constexpr double dLambda_initial = 0.5;             //Initial step size for ODE solvers with adaptive step size control

#if REG == 2
// Vector with the values of U for which we have NRG data to compare with (exclude zero!)
// Attention: these values are in units of Delta/2, not Delta -> corresponding U_fRG values are twice as large!
const std::vector<double> U_NRG {0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1., 1.2, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3., 5.};                                                    // NOLINT(cert-err58-cpp)
#else
const std::vector<double> U_NRG {};
#endif






/// Set flags used in code; DO NOT TOUCH!!!///

#ifdef HUBBARD
constexpr bool HUBBARD_MODEL = true;
#else
constexpr bool HUBBARD_MODEL = false;
#endif

#ifdef FERMI_POLARON_PROBLEM
constexpr bool FPP = true;
#else
constexpr bool FPP = false;
#endif

#if KELDYSH_FORMALISM
constexpr bool KELDYSH = true;
#else
constexpr bool KELDYSH = false;
#endif // KELDYSH_FORMALISM
#if ZERO_TEMP
constexpr bool ZERO_T = true;
#else
constexpr bool ZERO_T = false;
#endif // ZERO_TEMP

#if PARTICLE_HOLE_SYMM
constexpr bool PARTICLE_HOLE_SYMMETRY = true;
#else
constexpr bool PARTICLE_HOLE_SYMMETRY = false;
#endif

#if not KELDYSH_FORMALISM and not ZERO_TEMP
    using freqType = double;
#else
    using freqType = double;
#endif

inline std::string data_dir;


#include "frequency_parameters.hpp"
#include "technical_parameters.hpp"
#include "momentum_parameters.hpp"

#endif //KELDYSH_MFRG_PARAMETERS_H

    """
    filename = "./parameters/master_parameters.hpp"
    with open(filename, 'w') as f:
        f.write(text)

if __name__ == '__main__':
    backups = backup_files()

    write_main()
    write_frequency_parameters()
    write_technical_parameters()

    for isK in [0, 1]:
        for isT0 in [0,1]:
            for isPHS in [0,1]:
                for isSBE in [0]: #[0,1]:
                    for REG in [2]: #[2,4]
                        write_master_parameters(isK, isT0, isPHS, isSBE, REG)
                        os.system("bash ./scripts/compile_cluster.sh --ASC")
                        os.system("./Keldysh_mfRG 1 1")

    restore_backup_files(backups)
