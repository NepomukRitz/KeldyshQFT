#ifndef KELDYSH_MFRG_BUBBLE_FUNCTION_HPP
#define KELDYSH_MFRG_BUBBLE_FUNCTION_HPP

#include <cmath>                            // for using the macro M_PI as pi
#include "../symmetries/Keldysh_symmetries.hpp"  // for independent Keldysh components and utilities
#include "../correlation_functions/four_point/vertex.hpp"                         // vertex class
#include "../correlation_functions/two_point/selfenergy.hpp"                     // self-energy class
#include "../correlation_functions/two_point/propagator.hpp"                     // propagator class
#include "../integrator/integrator.hpp"          // integration routines
#include "../utilities/util.hpp"                 // measuring time, printing text output
#include "../utilities/mpi_setup.hpp"            // mpi parallelization routines
#include "../asymptotic_corrections/correction_functions.hpp"            // correction terms due to finite integration range
#include "../utilities/write_data2file.hpp"      // write vectors into hdf5 file
#include "../grids/momentum_grid.hpp"            // Momentum grid specific to the 2D Hubbard model
#include "../correlation_functions/four_point/vertex_data.hpp"
#include "bubble.hpp"
#include "precalculated_bubble.hpp"
#include "integrand.hpp"


template <typename Q,
        symmetryType symmetry_result,
        symmetryType symmetry_left,
        symmetryType symmetry_right,
        class Bubble_Object>
class BubbleFunctionCalculator{
    private:
    GeneralVertex<Q, symmetry_result>& dgamma;

    const GeneralVertex<Q, symmetry_left>& vertex1;
    const GeneralVertex<Q, symmetry_right>& vertex2;

    const Bubble_Object& Pi;
    const char channel;
    const bool diff = Pi.diff;

    const double Delta = (Pi.g.Lambda + glb_Gamma) / 2.; // hybridization (needed for proper splitting of the integration domain)

    int nw1_w = 0, nw2_w = 0, nw2_v = 0, nw3_w = 0, nw3_v = 0, nw3_v_p = 0;
    Q prefactor = 1.;

    int mpi_size = mpi_world_size(); // number of mpi processes
    int mpi_rank = mpi_world_rank(); // number of the current mpi process

    double vmin = 0, vmax = 0;
    int Nmin, Nmax; // Matsubara indices for minimal and maximal frequency. Only needed for finite-temperature Matsubara calculations!

    double tK1 = 0, tK2 = 0, tK3 = 0;

    void set_channel_specific_freq_ranges_and_prefactor();
    void find_vmin_and_vmax();

    bool missing_cross_projection(); // Needed for the Hubbard model.

    void calculate_bubble_function(K_class diag_class);
    Q get_value(int i_mpi, int i_omp, int n_omp, K_class diag_class);
    void calculate_value(Q &value, int i0, int i_in,  int ispin, int iw, double w, double v, double vp, K_class k);

    void write_out_results(const vec<Q>& Ordered_result, K_class diag_class);
    void write_out_results_K1(const vec<Q>& K1_ordered_result);
    void write_out_results_K2(const vec<Q>& K2_ordered_result);
    void write_out_results_K2b(const vec<Q>& K2b_ordered_result);
    void write_out_results_K3(const vec<Q>& K3_ordered_result);

    void set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, K_class diag_class);

    void convert_external_MPI_OMP_indices_to_physical_indices_K1(int& iK1, int& i0, int& ispin, int& iw, int& i_in, double& w,
                                                                 int i_mpi, int n_omp, int i_omp);
    void convert_external_MPI_OMP_indices_to_physical_indices_K2(int& iK2, int& i0, int& ispin, int& iw, int& iv, int& i_in,
                                                                 double& w, double& v,
                                                                 int i_mpi, int n_omp, int i_omp);
    void convert_external_MPI_OMP_indices_to_physical_indices_K2b(int& iK2, int& i0, int& ispin, int& iw, int& ivp, int& i_in,
                                                                 double& w, double& vp,
                                                                 int i_mpi, int n_omp, int i_omp);
    void convert_external_MPI_OMP_indices_to_physical_indices_K3(int& iK3, int& i0, int& ispin, int& iw, int& iv, int& ivp, int& i_in,
                                                                 double& w, double& v, double& vp,
                                                                 int i_mpi, int n_omp, int i_omp);

    int get_trafo_K1(int i0, double w);
    int get_trafo_K2(int i0, double w, double v);
    int get_trafo_K3(int i0, double w, double v, double vp);

    void get_Matsubara_integration_intervals(size_t& num_intervals, vec<vec<double>>& intervals, double w);

    Q bubble_value_prefactor();

    public:
    void perform_computation();

    BubbleFunctionCalculator(GeneralVertex<Q, symmetry_result>& dgamma_in,
                             const GeneralVertex<Q, symmetry_left>& vertex1_in,
                             const GeneralVertex<Q, symmetry_right>& vertex2_in,
                             const Bubble_Object& Pi_in,
                             const char channel_in)
                             :dgamma(dgamma_in), vertex1(vertex1_in), vertex2(vertex2_in),
                             Pi(Pi_in), channel(channel_in){
        set_channel_specific_freq_ranges_and_prefactor();
        find_vmin_and_vmax();

        // For Hubbard model computations, make sure that the internal structures of the vertices are parametrized correctly.
        if (HUBBARD_MODEL && (MAX_DIAG_CLASS > 1) && missing_cross_projection()) { // Cross projected parts are only needed for the Hubbard model in K2 and K3.
            print("Error! Needed crossprojection still has to be computed. Abort.");
            assert(false);
        }

        /// TODO(high): Figure out computations which need gamma_a_uu = gamma_a_ud - gamma_t_ud in a t-bubble,
        ///  i.e. CP_to_t(gamma_a_uu) = CP_to_t(gamma_a_ud) - CP_to_a(gamma_t_ud).
        ///  The integrand will need vertex AND vertex_initial to have access to cross-projected parts and non-crossprojected parts.
    }
};

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
Bubble_Object>::set_channel_specific_freq_ranges_and_prefactor() {
// set channel-specific frequency ranges and prefactor (1, 1, -1 for a, p, t) for sum over spins.
    switch (channel) {
        case 'a':
            nw1_w = nw1_a;
            nw2_w = nw2_a;
            nw2_v = nv2_a;
            nw3_w = nw3_a;
            nw3_v = nv3_a;
            nw3_v_p = nv3_a;
            prefactor = 1.;
            break;
        case 'p':
            nw1_w = nw1_p;
            nw2_w = nw2_p;
            nw2_v = nv2_p;
            nw3_w = nw3_p;
            nw3_v = nv3_p;
            nw3_v_p = nv3_p;
            prefactor = 1.;
            break;
        case 't':
            nw1_w = nw1_t;
            nw2_w = nw2_t;
            nw2_v = nv2_t;
            nw3_w = nw3_t;
            nw3_v = nv3_t;
            nw3_v_p = nv3_t;
            prefactor = -1.;
            break;
        default: ;
    }
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>::find_vmin_and_vmax() {
    // use std::min/std::max of selfenergy/K1 frequency grids as integration limits
    if (HUBBARD_MODEL){ // In the HM we have a larger frequency box for the SE and want to limit us to the range of bosonic vertex frequencies.
        vmin = dgamma.avertex().K1.K1_get_wlower();
        vmax = dgamma.avertex().K1.K1_get_wupper();
    }
    else{
        vmin = std::min(dgamma.avertex().K1.K1_get_wlower(), Pi.g.selfenergy.frequencies.w_lower);
        vmax = std::max(dgamma.avertex().K1.K1_get_wupper(), Pi.g.selfenergy.frequencies.w_upper);
    }

    if (MAX_DIAG_CLASS >= 2){
        // use std::min/std::max of selfenergy/K1/K2 frequency grids as integration limits
        vmin = std::min(vmin, dgamma.avertex().K2.K2_get_wlower_f());
        vmax = std::max(vmax, dgamma.avertex().K2.K2_get_wupper_f());
    }
    if (MAX_DIAG_CLASS >= 3){
        // use std::min/std::max of selfenergy/K1/K2/K3 frequency grids as integration limits
        vmin = std::min(vmin, dgamma.avertex().K3.K3_get_wlower_f());
        vmax = std::max(vmax, dgamma.avertex().K3.K3_get_wupper_f());
    }
    if ((!KELDYSH) && (!ZERO_T)) { // for finite-temperature Matsubara calculations
        // make sure that the limits for the Matsubara sum are fermionic
        Nmin = (int) (vmin/(M_PI*glb_T)-1)/2;
        Nmax = - Nmin;
        vmin = (Nmin*2+1)*(M_PI*glb_T);
        vmax = (Nmax*2+1)*(M_PI*glb_T);
    }
}

template<typename Q,
        symmetryType symmetry_result,
        symmetryType symmetry_left,
        symmetryType symmetry_right,
        class Bubble_Object>
bool
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>::missing_cross_projection() {
    switch (channel) {
        case 'a':
            if (!vertex1.pvertex().calculated_crossprojections ||
                !vertex1.tvertex().calculated_crossprojections ||
                !vertex2.pvertex().calculated_crossprojections ||
                !vertex2.tvertex().calculated_crossprojections)
            {
                return true;
            }
            break;
        case 'p':
            if (!vertex1.avertex().calculated_crossprojections ||
                !vertex1.tvertex().calculated_crossprojections ||
                !vertex2.avertex().calculated_crossprojections ||
                !vertex2.tvertex().calculated_crossprojections)
            {
                return true;
            }
            break;
        case 't':
            if (!vertex1.avertex().calculated_crossprojections ||
                !vertex1.pvertex().calculated_crossprojections ||
                !vertex2.avertex().calculated_crossprojections ||
                !vertex2.pvertex().calculated_crossprojections)
            {
                return true;
            }
            break;
        default:;
    }
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::perform_computation(){

    double t_start = get_time();
    if (MAX_DIAG_CLASS >= 0) {
        calculate_bubble_function(k1);
        tK1 = get_time() - t_start;
        print("K1", channel, " done, ");
        get_time(t_start);
    }
    if (MAX_DIAG_CLASS >= 2) {
        t_start = get_time();
        calculate_bubble_function(k2);
        tK2 = get_time() - t_start;
        print("K2", channel, " done, ");
        get_time(t_start);

#ifdef DEBUG_SYMMETRIES
        t_start = get_time();
        calculate_bubble_function(k2b);
        tK2 = get_time() - t_start;
        print("K2b", channel, " done, ");
        get_time(t_start);
#endif
    }
    if (MAX_DIAG_CLASS >= 3) {
        t_start = get_time();
        calculate_bubble_function(k3);
        tK3 = get_time() - t_start;
        print("K3", channel, " done, ");
        get_time(t_start);
    }

}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_bubble_function(const K_class diag_class){
    if (diag_class < k1 || diag_class > k3){print("Incompatible diagrammatic class! Abort."); assert(false); return;}

    int n_mpi, n_omp;
    set_external_arguments_for_parallelization(n_mpi, n_omp, diag_class);


    // initialize buffer into which each MPI process writes their results
    vec<Q> Buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);
    vertex1.initializeInterpol();
    vertex2.initializeInterpol();

    // start for-loop over external arguments, using MPI and OMP
    int iterator = 0;
    for (int i_mpi = 0; i_mpi < n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for schedule(dynamic) default(none) shared(n_omp, i_mpi, iterator, Buffer) // ,diag_class ) // needed for KCS
            for (int i_omp = 0; i_omp < n_omp; ++i_omp) {
                Buffer[iterator*n_omp + i_omp] = get_value(i_mpi, i_omp, n_omp, diag_class); // write result of integration into MPI buffer
            }
            ++iterator;
        }
    }
    // collect+combine results from different MPI processes, reorder them appropriately
    vec<Q> Result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(Buffer, Result, n_mpi, n_omp);
    vec<Q> Ordered_result = mpi_reorder_result(Result, n_mpi, n_omp);

    write_out_results(Ordered_result, diag_class);

    vertex1.set_initializedInterpol(false);
    vertex2.set_initializedInterpol(false);

}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
Q
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_value(const int i_mpi, const int i_omp, const int n_omp, const K_class diag_class){
    Q value = 0.;
    int iK1, iK2, i0, ispin, iw, iv, ivp, i_in;
    double w, v, vp;
    int trafo;
    switch (diag_class) {
        case k1:
            convert_external_MPI_OMP_indices_to_physical_indices_K1(iK1, i0, ispin,iw, i_in, w,
                                                                    i_mpi, n_omp, i_omp);
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
#else
            trafo = get_trafo_K1(i0, w);
#endif
            if (trafo == 0 || HUBBARD_MODEL) {

                calculate_value(value, i0, i_in, ispin, 0, w, 0, 0, k1);
            } // TODO: Freqency symmetries for the Hubbard model?
            break;
        case k2:
            convert_external_MPI_OMP_indices_to_physical_indices_K2(iK2, i0, ispin, iw, iv, i_in, w, v,
                                                                    i_mpi, n_omp, i_omp);
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
#else
            trafo = get_trafo_K2(i0, w, v);
#endif
            if (trafo == 0 || HUBBARD_MODEL) {calculate_value(value, i0, i_in, ispin, 0, w, v, 0, k2); }
            break;
#ifdef DEBUG_SYMMETRIES
        case k2b:
            convert_external_MPI_OMP_indices_to_physical_indices_K2b(iK2, i0, ispin, iw, ivp, i_in, w, vp,
                                                                    i_mpi, n_omp, i_omp);
            trafo = 0; // compute integrals for all frequency components
#else
            //trafo = get_trafo_K2(i0, w, v);
#endif
            if (trafo == 0) {calculate_value(value, i0, i_in, ispin, 0, w, 0, vp, k2b); }
            break;
        case k3:
            convert_external_MPI_OMP_indices_to_physical_indices_K3(iK2, i0, ispin, iw, iv, ivp, i_in, w, v, vp,
                                                                    i_mpi, n_omp, i_omp);
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
#else
            trafo = get_trafo_K3(i0, w, v, vp);
#endif
            if (trafo == 0 || HUBBARD_MODEL) {
                if (channel == 'a') {calculate_value(value, i0, i_in, ispin, iw, w, v, vp, k3); } // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
                if (channel == 'p') {calculate_value(value, i0, i_in, ispin, iv, w, v, vp, k3); } // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
                if (channel == 't') {calculate_value(value, i0, i_in, ispin, ivp,w, v, vp, k3); } // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
            }
            break;
        default:;
    }
    return value;
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_value(Q& value, const int i0, const int i_in, const int ispin, const int iw,
                                           const double w, const double v, const double vp, const K_class k){
#ifndef SWITCH_SUM_N_INTEGRAL
    for (int i2 : glb_non_zero_Keldysh_bubble) {
        int n_spin_sum = 1;                  // number of summands in spin sum (=1 in the a channel)
        if ((channel == 't' and ispin == 0) or (channel == 'a' and ispin == 1)) n_spin_sum = 3;  // in the t channel, spin sum includes three terms
        for (int i_spin=0; i_spin < n_spin_sum; ++i_spin) {
#else
            int i2 = 0;
            int i_spin = 0;
#endif
            // initialize the integrand object and perform frequency integration
            Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>
                    integrand(vertex1, vertex2, Pi, i0, i2, ispin, iw, w, v, vp, i_in, i_spin, channel, diff, k);
            if (KELDYSH) {
                value += bubble_value_prefactor() * integrator<Q>(integrand, vmin, vmax, -w / 2., w / 2., Delta);
            }
            else {
                if (ZERO_T) {
                    switch (k) {
                        case k1:
                            value += bubble_value_prefactor() *
                                     integrator_Matsubara_T0<Q, 0>(integrand, vmin, vmax, std::abs(w / 2),
                                                                   {}, Delta, false);
                            break;
                        case k2:
                            value += bubble_value_prefactor() *
                                     integrator_Matsubara_T0<Q, 3>(integrand, vmin, vmax, std::abs(w / 2),
                                                                   {v, v + w, v - w}, Delta, false);
                            break;
                        case k3:
                            value += bubble_value_prefactor() *
                                     integrator_Matsubara_T0<Q, 6>(integrand, vmin, vmax, std::abs(w / 2),
                                                                   {v, vp, w - vp, w + vp, w - v, abs(w) + abs(v)}, Delta,
                                                                   false);
                            break;
                        case k2b:
                            value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,3>(integrand, vmin, vmax, std::abs(w/2), {vp, vp+w, vp-w}, Delta, false);
                            break;
                        default:
                            break;
                    }
                }
                else {
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP) // TODO(high): Figure out type problems in matsubarasum
                    int interval_correction =  (int)(signFlipCorrection_MF(w)/(2*M_PI*glb_T) + 0.1);
                    // if interval_correction=-1, then the integrand is symmetric around v=-M_PI*glb_T
                    value += bubble_value_prefactor()*(2*M_PI) * glb_T * matsubarasum<Q>(integrand, Nmin, Nmax  + interval_correction);
#endif
                }
            }

#ifndef SWITCH_SUM_N_INTEGRAL
        }
#else
        for (int i2 : glb_non_zero_Keldysh_bubble) {
#endif
        // asymptotic corrections include spin sum
        if (not HUBBARD_MODEL) value += bubble_value_prefactor() *
                                        asymp_corrections_bubble(k, vertex1, vertex2, Pi.g,
                                                                 vmin, vmax, w, v, vp, i0, i2, i_in, channel, diff, ispin);

    }
}


template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results(const vec<Q>& Ordered_result, const K_class diag_class){
    switch (diag_class) {
        case k1:
            write_out_results_K1(Ordered_result);
            break;
        case k2:
            write_out_results_K2(Ordered_result);
            break;
#ifdef DEBUG_SYMMETRIES
        case k2b:
            write_out_results_K2b(Ordered_result);
            break;
#endif
        case k3:
            write_out_results_K3(Ordered_result);
            break;
        default: ;
    }
    dgamma.set_initializedInterpol(false);      // above initialization of the Interpolator is with the symmetry-reduced sector only (rest = zero)
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::write_out_results_K1(const vec<Q>& K1_ordered_result){
    switch (channel) {
        case 'a':
            dgamma.avertex().K1.add_vec(K1_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.avertex().enforce_freqsymmetriesK1(dgamma.avertex());
#endif
            break;
        case 'p':
            dgamma.pvertex().K1.add_vec(K1_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.pvertex().enforce_freqsymmetriesK1(dgamma.pvertex());
#endif
            break;
        case 't':
            dgamma.tvertex().K1.add_vec(K1_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.tvertex().enforce_freqsymmetriesK1(dgamma.tvertex());
#endif
            break;
        default: ;
    }
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K2(const vec<Q>& K2_ordered_result){
    assert( K2_ordered_result.size() == nK_K2*n_spin*nBOS2*nFER2*n_in);
    switch (channel) {
        case 'a':
            dgamma.avertex().K2.add_vec(K2_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.avertex().enforce_freqsymmetriesK2(dgamma.avertex());
#endif
            break;
        case 'p':
            dgamma.pvertex().K2.add_vec(K2_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.pvertex().enforce_freqsymmetriesK2(dgamma.pvertex());
#endif
            break;
        case 't':
            dgamma.tvertex().K2.add_vec(K2_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.tvertex().enforce_freqsymmetriesK2(dgamma.tvertex());
#endif
            break;
        default: ;
    }
#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and defined(USE_FDT)
    compute_components_through_FDTs(dgamma.half1(), dgamma.half1(), dgamma.half1(), channel);
#endif
}

#ifdef DEBUG_SYMMETRIES
template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K2b(const vec<Q>& K2b_ordered_result){
    assert( K2b_ordered_result.size() == nK_K2*n_spin*nBOS2*nFER2*n_in);
    switch (channel) {
        case 'a':
            dgamma.avertex().K2b.add_vec(K2b_ordered_result);
#ifndef DEBUG_SYMMETRIES
            //dgamma.avertex().enforce_freqsymmetriesK2(dgamma.avertex());
#endif
            break;
        case 'p':
            dgamma.pvertex().K2b.add_vec(K2b_ordered_result);
#ifndef DEBUG_SYMMETRIES
            //dgamma.pvertex().enforce_freqsymmetriesK2(dgamma.pvertex());
#endif
            break;
        case 't':
            dgamma.tvertex().K2b.add_vec(K2b_ordered_result);
#ifndef DEBUG_SYMMETRIES
            //dgamma.tvertex().enforce_freqsymmetriesK2(dgamma.tvertex());
#endif
            break;
        default: ;
    }
}
#endif

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K3(const vec<Q>& K3_ordered_result){
    switch (channel) {
        case 'a':
            dgamma.avertex().K3.add_vec(K3_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.avertex().enforce_freqsymmetriesK3(dgamma.avertex());
#endif
            break;
        case 'p':
            dgamma.pvertex().K3.add_vec(K3_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.pvertex().enforce_freqsymmetriesK3(dgamma.pvertex());
#endif
            break;
        case 't':
            dgamma.tvertex().K3.add_vec(K3_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
            if (not HUBBARD_MODEL) dgamma.tvertex().enforce_freqsymmetriesK3(dgamma.tvertex());
#endif
            break;
        default: ;
    }
#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and defined(USE_FDT)
    compute_components_through_FDTs(dgamma.half1(), dgamma.half1(), dgamma.half1(), channel);
#endif
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, const K_class diag_class){
    switch (diag_class) {
        case k1:
            n_mpi = nK_K1 * n_spin;        // set external arguments for MPI-parallelization (# of tasks distributed via MPI)
            n_omp = nw1_w * n_in_K1; // set external arguments for OMP-parallelization (# of tasks per MPI-task distributed via OMP)
            break;
        case k2:
        case k2b:
            n_mpi = nK_K2 * n_spin * nw2_w;
            n_omp = nw2_v * n_in_K2;
            break;
        case k3:
            n_mpi = nK_K3 * n_spin * nw3_w;
            n_omp = nw3_v * nw3_v_p * n_in_K3;
            break;
        default: ;
    }
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K1(int& iK1, int& i0, int& ispin, int& iw, int& i_in, double& w,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK1 = i_mpi * n_omp + i_omp;
    //i0 = iK1/(nw1_w*n_in_K1);                              // exterior Keldysh indices of the bubble
    //iw = iK1/(n_in_K1) - i0*nw1_w;                         // frequency index
    //i_in = iK1 - i0*nw1_w*n_in_K1 - iw*n_in_K1;            // internal index
    getMultIndex<4,int,int,int,int>(i0, ispin, iw, i_in, iK1, vertex1.avertex().K1.get_dims());

    if (channel == 'a') dgamma.avertex().K1.K1_get_freq_w(w, iw);           // frequency acc. to frequency index
    if (channel == 'p') dgamma.pvertex().K1.K1_get_freq_w(w, iw);           // frequency acc. to frequency index
    if (channel == 't') dgamma.tvertex().K1.K1_get_freq_w(w, iw);           // frequency acc. to frequency index
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2(int& iK2, int& i0, int& ispin,  int& iw, int& iv, int& i_in,
                                                                                double& w, double& v,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;
    //i0 = iK2 / (nw2_w * nw2_v * n_in_K2);
    //iw = iK2 / (nw2_v * n_in_K2) - i0 * nw2_w;
    //iv = iK2 / n_in_K2 - iw * nw2_v - i0 * nw2_w * nw2_v;
    //i_in = iK2 - iv * n_in_K2 - iw * nw2_v * n_in_K2 - i0 * nw2_w * nw2_v * n_in_K2;
    getMultIndex<5,int,int,int,int,int>(i0, ispin, iw, iv, i_in, iK2, vertex1.avertex().K2.get_dims());
    if (channel == 'a') dgamma.avertex().K2.K2_get_freqs_w(w, v, iw, iv);
    if (channel == 'p') dgamma.pvertex().K2.K2_get_freqs_w(w, v, iw, iv);
    if (channel == 't') dgamma.tvertex().K2.K2_get_freqs_w(w, v, iw, iv);
}

#ifdef DEBUG_SYMMETRIES
template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2b(int& iK2, int& i0, int& ispin,  int& iw, int& ivp, int& i_in,
                                                                                double& w, double& vp,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;
    //i0 = iK2 / (nw2_w * nw2_v * n_in_K2);
    //iw = iK2 / (nw2_v * n_in_K2) - i0 * nw2_w;
    //ivp= iK2 / n_in_K2 - iw * nw2_v - i0 * nw2_w * nw2_v;
    //i_in = iK2 - ivp * n_in_K2 - iw * nw2_v * n_in_K2 - i0 * nw2_w * nw2_v * n_in_K2;
    getMultIndex<5,int,int,int,int,int>(i0, ispin, iw, ivp, i_in, iK2, vertex1.avertex().K2b.get_dims());
    if (channel == 'a') dgamma.avertex().K2b.K2_get_freqs_w(w, vp, iw, ivp);
    if (channel == 'p') dgamma.pvertex().K2b.K2_get_freqs_w(w, vp, iw, ivp);
    if (channel == 't') dgamma.tvertex().K2b.K2_get_freqs_w(w, vp, iw, ivp);
}
#endif

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K3(int& iK3, int& i0, int& ispin,  int& iw, int& iv, int& ivp, int& i_in,
                                                                                     double& w, double& v, double& vp,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK3 = i_mpi * n_omp + i_omp;
    //i0 = iK3/(nw3_w * nw3_v * nw3_v_p * n_in_K3);
    //iw = iK3/(nw3_v * nw3_v_p * n_in_K3) - i0*nw3_w;
    //iv = iK3/(nw3_v * n_in_K3) - i0*nw3_w*nw3_v - iw*nw3_v;
    //ivp =iK3/(n_in_K3) - i0*nw3_w*nw3_v*nw3_v_p - iw*nw3_v*nw3_v_p - iv*nw3_v_p;
    //i_in = iK3 - i0*nw3_w*nw3_v*nw3_v_p*n_in_K3 - iw*nw3_v*nw3_v_p*n_in_K3 - iv*nw3_v_p*n_in_K3 - ivp*n_in_K3;
    getMultIndex<6,int,int,int,int,int,int>(i0, ispin, iw, iv, ivp, i_in, iK3, vertex1.avertex().K3.get_dims());
    if (channel == 'a') dgamma.avertex().K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, 'a');
    if (channel == 'p') dgamma.pvertex().K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, 'p');
    if (channel == 't') dgamma.tvertex().K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, 't');
}


template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
int
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::get_trafo_K1(const int i0, const double w){
    int trafo = 1;
    const double safety = 1e-10;
    int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
    switch (channel) {
        case 'a':
            trafo = TransformaK1a[i0][sign_w];
            //cout << "Ping!" << trafo << "\n";
            break;
        case 'p':
            trafo = TransformaK1p[i0][sign_w];
            break;
        case 't':
            trafo = TransformaK1t[i0][sign_w];
            break;
        default:
            print("Something went wrong in get_trafo_K1! Abort."); assert(false);
    }
    return trafo;
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
int
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_trafo_K2(const int i0, const double w, const double v){
    int trafo = 1;
    const double safety = 1e-10;
    int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
    int sign_v = sign_index<double>(v - safety); // safety to ensure that w=0 gets sign_w=-1
    switch (channel) {
        case 'a':
            trafo = TransformaK2a[i0][sign_w * 2 + sign_v];
            //cout << "Ping!" << trafo << "\n";
            break;
        case 'p':
            trafo = TransformaK2p[i0][sign_w * 2 + sign_v];
            break;
        case 't':
            trafo = TransformaK2t[i0][sign_w * 2 + sign_v];
            break;
        default:
            print("Something went wrong in get_trafo_K2! Abort."); assert(false);
    }
#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and defined(USE_FDT)
    switch (channel) {
        case 'a':
            if ((i0 == 0 or i0 == 2 or i0 == 3) and abs(w)>glb_T*26.) trafo = -1; // components can be determined via FDTs, no need to compute it via integration
            break;
        case 'p':
            if ((i0 == 0 or i0 == 1 or i0 == 3) and abs(w)>glb_T*26.) trafo = -1; // components can be determined via FDTs, no need to compute it via integration
            break;
        case 't':
            if ((i0 == 0 or i0 == 2 or i0 == 3) and abs(w)>glb_T*26.) trafo = -1; // components can be determined via FDTs, no need to compute it via integration
            break;
        default:
            break;
    }
#endif
    return trafo;
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
int
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_trafo_K3(const int i0, const double w, const double v, const double vp){
    int trafo = 1;
    const double safety = 1e-10;
    int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
    int sign_f = sign_index(v + vp - safety);
    int sign_fp= sign_index(v - vp - safety);
    switch (channel) {
        case 'a':
            trafo = TransformaK3a[i0][sign_w * 4 + sign_f * 2 + sign_fp];
            //cout << "Ping!" << trafo << "\n";
            break;
        case 'p':
            trafo = TransformaK3p[i0][sign_w * 4 + sign_f * 2 + sign_fp];
            break;
        case 't':
            trafo = TransformaK3t[i0][sign_w * 4 + sign_f * 2 + sign_fp];
            break;
        default:
            print("Something went wrong in get_trafo_K3! Abort."); assert(false);
    }
#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and defined(USE_FDT)
    if (i0 == 0 or i0 == 1) trafo = -1; // components can be determined via FDTs, no need to compute it via integration
#endif
    return trafo;
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::get_Matsubara_integration_intervals(size_t& num_intervals, vec<vec<double>>& intervals,
                                                                 const double w){
    if( -std::abs(w/2)+inter_tol < std::abs(w/2)-inter_tol){
        intervals = {{vmin, -std::abs(w/2)-inter_tol}, {-std::abs(w/2)+inter_tol, std::abs(w/2)-inter_tol}, {std::abs(w/2)+inter_tol, vmax}};
        num_intervals = 3;
    }
    else {
        intervals = {{vmin, -std::abs(w/2)-inter_tol}, {std::abs(w/2)+inter_tol, vmax}};
        num_intervals = 2;
    }
}

template<typename Q, symmetryType symmetry_result, symmetryType symmetry_left,
        symmetryType symmetry_right, class Bubble_Object>
Q
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::bubble_value_prefactor(){
    if constexpr (KELDYSH) return prefactor * (1. / (2. * M_PI * glb_i));
    else                   return prefactor * (1. / (2. * M_PI));
}



// bubble_function using the new class BubbleFunctionCalculator
template <typename Q,
        symmetryType symmetry_result,
        symmetryType symmetry_left,
        symmetryType symmetry_right,
        class Bubble_Object>
void bubble_function(GeneralVertex<Q, symmetry_result>& dgamma,
                                 const GeneralVertex<Q, symmetry_left>& vertex1,
                                 const GeneralVertex<Q, symmetry_right>& vertex2,
                                 const Bubble_Object& Pi,
                                 const char channel){
    BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>
            BubbleComputer (dgamma, vertex1, vertex2, Pi, channel);
    if (channel == 'a' or channel == 'p' or channel == 't') BubbleComputer.perform_computation();
    //else {print("Error! Incompatible channel given to bubble_function. Abort"); }

}

/// Overload for bubble_function in case no Bubble object has been initialized yet. ONLY WORKS FOR SIAM!!
template <typename Q,
        symmetryType symmetry_result,
        symmetryType symmetry_left,
        symmetryType symmetry_right>
void bubble_function(GeneralVertex<Q, symmetry_result>& dgamma,
                     const GeneralVertex<Q, symmetry_left>& vertex1,
                     const GeneralVertex<Q, symmetry_right>& vertex2,
                     const Propagator<Q>& G, const Propagator<Q>& S, const char channel, const bool diff){
    Bubble<Q> Pi(G, S, diff);
    bubble_function(dgamma, vertex1, vertex2, Pi, channel);
}

#endif //KELDYSH_MFRG_BUBBLE_FUNCTION_HPP