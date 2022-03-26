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
#include "bubble.hpp"
#include "precalculated_bubble.hpp"
#include "integrand.hpp"


template <char channel,
        typename Q,
        vertexType symmetry_result,
        vertexType symmetry_left,
        vertexType symmetry_right,
        class Bubble_Object>
class BubbleFunctionCalculator{
    private:
    GeneralVertex<Q, symmetry_result>& dgamma;

    const GeneralVertex<Q, symmetry_left> vertex1;  /// THIS IS A COPY; needed for symmetry-expansion
    const GeneralVertex<Q, symmetry_right> vertex2; /// THIS IS A COPY; needed for symmetry-expansion

    const Bubble_Object& Pi;
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

    template<K_class diag_class> void calculate_bubble_function();
    template<K_class diag_class> Q get_value(int i_mpi, int i_omp, int n_omp);
    template<K_class diag_class,int spin> void calculate_value(Q &value, int i0, int i_in, int iw, double w, double v, double vp);

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
                             const Bubble_Object& Pi_in)
                             :dgamma(dgamma_in), vertex1(vertex1_in), vertex2(vertex2_in),
                             Pi(Pi_in){
        set_channel_specific_freq_ranges_and_prefactor();
        find_vmin_and_vmax();

        // For Hubbard model computations, make sure that the internal structures of the vertices are parametrized correctly.
        if (HUBBARD_MODEL && (MAX_DIAG_CLASS > 1) && missing_cross_projection()) { // Cross projected parts are only needed for the Hubbard model in K2 and K3.
            print("Error! Needed crossprojection still has to be computed. Abort.");
            assert(false);
        }
#ifdef SWITCH_SUM_N_INTEGRAL
        vertex1.template symmetry_expand<channel,true>();
        vertex2.template symmetry_expand<channel,false>();
#endif
        /// TODO(high): Figure out computations which need gamma_a_uu = gamma_a_ud - gamma_t_ud in a t-bubble,
        ///  i.e. CP_to_t(gamma_a_uu) = CP_to_t(gamma_a_ud) - CP_to_a(gamma_t_ud).
        ///  The integrand will need vertex AND vertex_initial to have access to cross-projected parts and non-crossprojected parts.
    }
};

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>::find_vmin_and_vmax() {
    // use std::min/std::max of selfenergy/K1 frequency grids as integration limits
    if (HUBBARD_MODEL){ // In the HM we have a larger frequency box for the SE and want to limit us to the range of bosonic vertex frequencies.
        vmin = dgamma.avertex().K1.frequencies.get_wupper_b();
        vmax = dgamma.avertex().K1.frequencies.get_wupper_b();
    }
    else{
        vmin = std::min(dgamma.avertex().K1.frequencies.get_wupper_b(), Pi.g.selfenergy.Sigma.frequencies.primary_grid.w_lower);
        vmax = std::max(dgamma.avertex().K1.frequencies.get_wupper_b(), Pi.g.selfenergy.Sigma.frequencies.primary_grid.w_upper);
    }

    if constexpr(MAX_DIAG_CLASS >= 2){
        // use std::min/std::max of selfenergy/K1/K2 frequency grids as integration limits
        vmin = std::min(vmin, dgamma.avertex().K2.frequencies.get_wlower_f());
        vmax = std::max(vmax, dgamma.avertex().K2.frequencies.get_wupper_f());
    }
    if constexpr(MAX_DIAG_CLASS >= 3){
        // use std::min/std::max of selfenergy/K1/K2/K3 frequency grids as integration limits
        vmin = std::min(vmin, dgamma.avertex().K3.frequencies.get_wlower_f());
        vmax = std::max(vmax, dgamma.avertex().K3.frequencies.get_wupper_f());
    }
    if constexpr((!KELDYSH) && (!ZERO_T)) { // for finite-temperature Matsubara calculations
        // make sure that the limits for the Matsubara sum are fermionic
        Nmin = - POSINTRANGE; // (int) (vmin/(M_PI*glb_T)-1)/2;
        Nmax = - Nmin - 1;
        vmin = (Nmin*2+1)*(M_PI*glb_T);
        vmax = (Nmax*2+1)*(M_PI*glb_T);
    }
}

template<char channel, typename Q,
        vertexType symmetry_result,
        vertexType symmetry_left,
        vertexType symmetry_right,
        class Bubble_Object>
bool
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>::missing_cross_projection() {
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::perform_computation(){

    double t_start = get_time();
    if constexpr(MAX_DIAG_CLASS >= 0) {
        calculate_bubble_function<k1>();
        tK1 = get_time() - t_start;
        print("K1", channel, " done, ");
        get_time(t_start);
    }
    if constexpr(MAX_DIAG_CLASS >= 2) {
        t_start = get_time();
        calculate_bubble_function<k2>();
        tK2 = get_time() - t_start;
        print("K2", channel, " done, ");
        get_time(t_start);

#ifdef DEBUG_SYMMETRIES
        t_start = get_time();
        calculate_bubble_function<k2b>();
        tK2 = get_time() - t_start;
        print("K2b", channel, " done, ");
        get_time(t_start);
#endif
    }
    if constexpr(MAX_DIAG_CLASS >= 3) {
        t_start = get_time();
        calculate_bubble_function<k3>();
        tK3 = get_time() - t_start;
        print("K3", channel, " done, ");
        get_time(t_start);
    }

}

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
                template<K_class diag_class>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_bubble_function(){
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
#pragma omp parallel for schedule(dynamic)
            for (int i_omp = 0; i_omp < n_omp; ++i_omp) {
                Buffer[iterator*n_omp + i_omp] = get_value<diag_class>(i_mpi, i_omp, n_omp); // write result of integration into MPI buffer
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
template<K_class diag_class>
Q
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_value(const int i_mpi, const int i_omp, const int n_omp){
    Q value = 0.;
    int iK1, iK2, i0, ispin, iw, iv, ivp, i_in;
    double w, v, vp;
    int trafo;
    switch (diag_class) {
        case k1:
            convert_external_MPI_OMP_indices_to_physical_indices_K1(iK1, i0, ispin,iw, i_in, w,
                                                                    i_mpi, n_omp, i_omp);

            trafo = get_trafo_K1(i0, w);
            if (trafo == 0 || HUBBARD_MODEL) {
                if (ispin == 0 or n_spin == 1)  calculate_value<diag_class,0>(value, i0, i_in, 0, w, 0, 0);
                else                            calculate_value<diag_class,1>(value, i0, i_in, 0, w, 0, 0);
            } // TODO: Freqency symmetries for the Hubbard model?
            break;
        case k2:
            convert_external_MPI_OMP_indices_to_physical_indices_K2(iK2, i0, ispin, iw, iv, i_in, w, v,
                                                                    i_mpi, n_omp, i_omp);

            trafo = get_trafo_K2(i0, w, v);
            if (trafo == 0 || HUBBARD_MODEL) {
                if (ispin == 0 or n_spin == 1) calculate_value<diag_class,0>(value, i0, i_in, 0, w, v, 0);
                else                           calculate_value<diag_class,1>(value, i0, i_in, 0, w, v, 0);}
            break;
#ifdef DEBUG_SYMMETRIES
        case k2b:
            convert_external_MPI_OMP_indices_to_physical_indices_K2b(iK2, i0, ispin, iw, ivp, i_in, w, vp,
                                                                    i_mpi, n_omp, i_omp);
            trafo = 0; // compute integrals for all frequency components
            if (!KELDYSH and !ZERO_T and -vp + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K2b.frequencies.get_wlower_f()) {
                trafo = -1;
            }
            if (trafo == 0) {
                if (ispin == 0 or n_spin == 1) calculate_value<diag_class,0>(value, i0, i_in, 0, w, 0, vp);
                else                           calculate_value<diag_class,1>(value, i0, i_in, 0, w, 0, vp);
            }
            break;
#endif
        case k3:
            convert_external_MPI_OMP_indices_to_physical_indices_K3(iK2, i0, ispin, iw, iv, ivp, i_in, w, v, vp,
                                                                    i_mpi, n_omp, i_omp);

            trafo = get_trafo_K3(i0, w, v, vp);
            if (trafo == 0 || HUBBARD_MODEL) {
                if (ispin == 0 or n_spin == 1) calculate_value<diag_class,0>(value, i0, i_in, channel=='a' ? iw : (channel=='p' ? iv : ivp), w, v, vp);  // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
                else                           calculate_value<diag_class,1>(value, i0, i_in, channel=='a' ? iw : (channel=='p' ? iv : ivp), w, v, vp);
            }
            break;
        default:;
    }
    return value;
}

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
template<K_class k, int spin>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_value(Q& value, const int i0, const int i_in, const int iw,
                                           const double w, const double v, const double vp){

    double vmin_temp = vmin;
    double vmax_temp = vmax;
#ifndef SWITCH_SUM_N_INTEGRAL
    static_assert(n_spin == 1, "SWITCH_SUM_N_INTEGRAL not ready for DEBUG_SYMMETRIES.");
    for (int i2 : glb_non_zero_Keldysh_bubble) {
        int n_spin_sum = 1;                  // number of summands in spin sum (=1 in the a channel)
        if (channel == 't') n_spin_sum = 3;  // in the t channel, spin sum includes three terms
        for (int i_spin=0; i_spin < n_spin_sum; ++i_spin) {
#else
            int i2 = 0;
            int i_spin = 0;
#endif
            // initialize the integrand object and perform frequency integration
            Integrand<k, channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object>
                    integrand(vertex1, vertex2, Pi, i0, i2, iw, w, v, vp, i_in, i_spin, diff);
            if constexpr(ZERO_T) {
                switch (k) {
                    case k1:
                        value += bubble_value_prefactor() *
                                 integrator_Matsubara_T0<Q, 0>(integrand, vmin, vmax, std::abs(w / 2),
                                                               {}, Delta, true);

                        break;
                    case k2:
                        value += bubble_value_prefactor() *
                                 integrator_Matsubara_T0<Q, 3>(integrand, vmin, vmax, std::abs(w / 2),
                                                               {v, v + w, v - w}, Delta, true);
                        break;
                    case k3:
                        value += bubble_value_prefactor() *
                                 integrator_Matsubara_T0<Q, 6>(integrand, vmin, vmax, std::abs(w / 2),
                                                               {v, vp, w - vp, w + vp, w - v, std::abs(w) + std::abs(v)}, Delta,
                                                               true);
                        break;
                    case k2b:
                        value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,3>(integrand, vmin, vmax, std::abs(w/2), {vp, vp+w, vp-w}, Delta, true);
                        break;
                    default:
                        break;
                }
            }
            else {
                if constexpr(KELDYSH) {
                    value += bubble_value_prefactor() * integrator<Q>(integrand, vmin, vmax, -w / 2., w / 2., Delta);
                }
                else {
                    int interval_correction = signFlipCorrection_MF_int(w);
                    int W = (int) (w / (2*M_PI*glb_T) + 0.1*sgn(w));
                    vmin_temp = (-POSINTRANGE - std::abs(W/2) + interval_correction) * 2 * M_PI * glb_T;
                    vmax_temp = (POSINTRANGE-1  + std::abs(W/2)) * 2 * M_PI * glb_T;
                    // if interval_correction=-1, then the integrand is symmetric_full around v=-M_PI*glb_T

                    value = bubble_value_prefactor()*(2*M_PI) * glb_T * matsubarasum<Q>(integrand, -POSINTRANGE - std::abs(W/2) + interval_correction, POSINTRANGE-1  + std::abs(W/2));
                }
            }

#ifndef SWITCH_SUM_N_INTEGRAL
        }
#else
        for (int i2 : glb_non_zero_Keldysh_bubble) {
#endif
        // asymptotic corrections include spin sum
        if ( !HUBBARD_MODEL and !ZERO_T) {
                value += bubble_value_prefactor() * asymp_corrections_bubble<channel>(k, vertex1, vertex2, Pi.g,
                                                                                  vmin_temp, vmax_temp, w, v, vp, i0, i2, i_in, diff, spin);
        }

    }
}


template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K2(const vec<Q>& K2_ordered_result){
    //assert( K2_ordered_result.size() == nK_K2*n_spin*nBOS2*nFER2*n_in);
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
template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K2b(const vec<Q>& K2b_ordered_result){
    //assert( K2b_ordered_result.size() == nK_K2*n_spin*nBOS2*nFER2*n_in);
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, const K_class diag_class){
    const int nK_K1 = channel == 'p' ? K1p_config.dims[my_defs::K1::keldysh] : K1at_config.dims[my_defs::K1::keldysh];
    const int nK_K2 = channel == 'p' ? K2p_config.dims[my_defs::K2::keldysh] : K2at_config.dims[my_defs::K2::keldysh];
    const int nK_K3 = K3_config.dims[my_defs::K3::keldysh];
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K1(int& iK1, int& i0, int& ispin, int& iw, int& i_in, double& w,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK1 = i_mpi * n_omp + i_omp;

    my_defs::K1::index_type idx;
    getMultIndex<rank_K1>(idx, iK1, vertex1.get_rvertex(channel).K1.get_dims());
    i0       = (int) idx[my_defs::K1::keldysh];
    ispin    = (int) idx[my_defs::K1::spin];
    iw       = (int) idx[my_defs::K1::omega];
    i_in     = (int) idx[my_defs::K1::internal];
   //getMultIndex<4,int,int,int,int>(ispin, iw, i0, i_in, iK1, vertex1.avertex().K1.get_dims());

    if (channel == 'a') dgamma.avertex().K1.frequencies.get_freqs_w(w, iw);           // frequency acc. to frequency index
    if (channel == 'p') dgamma.pvertex().K1.frequencies.get_freqs_w(w, iw);           // frequency acc. to frequency index
    if (channel == 't') dgamma.tvertex().K1.frequencies.get_freqs_w(w, iw);           // frequency acc. to frequency index
}

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2(int& iK2, int& i0, int& ispin,  int& iw, int& iv, int& i_in,
                                                                                double& w, double& v,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;

    my_defs::K2::index_type idx;
    getMultIndex<rank_K2>(idx, iK2, vertex1.get_rvertex(channel).K2.get_dims());
    i0       = (int) idx[my_defs::K2::keldysh];
    ispin    = (int) idx[my_defs::K2::spin];
    iw       = (int) idx[my_defs::K2::omega];
    iv       = (int) idx[my_defs::K2::nu];
    i_in     = (int) idx[my_defs::K2::internal];
    //getMultIndex<5,int,int,int,int,int>(ispin, iw, iv, i0, i_in, iK2, vertex1.avertex().K2.get_dims());
    if (channel == 'a') dgamma.avertex().K2.frequencies.get_freqs_w(w, v, iw, iv);
    if (channel == 'p') dgamma.pvertex().K2.frequencies.get_freqs_w(w, v, iw, iv);
    if (channel == 't') dgamma.tvertex().K2.frequencies.get_freqs_w(w, v, iw, iv);
}

#ifdef DEBUG_SYMMETRIES
template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2b(int& iK2, int& i0, int& ispin,  int& iw, int& ivp, int& i_in,
                                                                                double& w, double& vp,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;

    my_defs::K2::index_type idx;
    getMultIndex<rank_K2>(idx, iK2, vertex1.get_rvertex(channel).K2b.get_dims());
    i0       = (int) idx[my_defs::K2b::keldysh];
    ispin    = (int) idx[my_defs::K2b::spin];
    iw       = (int) idx[my_defs::K2b::omega];
    ivp      = (int) idx[my_defs::K2b::nup];
    i_in     = (int) idx[my_defs::K2b::internal];
    //getMultIndex<5,int,int,int,int,int>(ispin, iw, ivp, i0, i_in, iK2, vertex1.avertex().K2b.get_dims());
    if (channel == 'a') dgamma.avertex().K2b.frequencies.get_freqs_w(w, vp, iw, ivp);
    if (channel == 'p') dgamma.pvertex().K2b.frequencies.get_freqs_w(w, vp, iw, ivp);
    if (channel == 't') dgamma.tvertex().K2b.frequencies.get_freqs_w(w, vp, iw, ivp);
}
#endif

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K3(int& iK3, int& i0, int& ispin,  int& iw, int& iv, int& ivp, int& i_in,
                                                                                     double& w, double& v, double& vp,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK3 = i_mpi * n_omp + i_omp;

    my_defs::K3::index_type idx;
    getMultIndex<rank_K3>(idx, iK3, vertex1.get_rvertex(channel).K3.get_dims());
    i0       = (int) idx[my_defs::K3::keldysh];
    ispin    = (int) idx[my_defs::K3::spin];
    iw       = (int) idx[my_defs::K3::omega];
    iv       = (int) idx[my_defs::K3::nu];
    ivp      = (int) idx[my_defs::K3::nup];
    i_in     = (int) idx[my_defs::K3::internal];
    //getMultIndex<6,int,int,int,int,int,int>(ispin, iw, iv, ivp, i0, i_in, iK3, vertex1.avertex().K3.get_dims());
    if (channel == 'a') dgamma.avertex().K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp);
    if (channel == 'p') dgamma.pvertex().K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp);
    if (channel == 't') dgamma.tvertex().K3.frequencies.get_freqs_w(w, v, vp, iw, iv, ivp);
}


template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
int
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::get_trafo_K1(const int i0, const double w){
    int trafo = 1;
#ifdef DEBUG_SYMMETRIES
    trafo = 0; // compute integrals for all frequency components
#else
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
#if CONTOUR_BASIS == 1 and defined(ZERO_TEMP)
    if (is_zero_due_to_FDTs<k1>(i0, w, 0, 0, channel)) trafo = -1; // components zero according to FDTs
#endif // CONTOUR_BASIS
#endif
    return trafo;
}

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
int
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_trafo_K2(const int i0, const double w, const double v){
    int trafo = 1;
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
            if (!KELDYSH and !ZERO_T and -v + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K2.frequencies.get_wlower_f()) {
                trafo = -1;
            }
#else
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
    if constexpr(!KELDYSH and !ZERO_T and -v + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K2.frequencies.get_wlower_f()) {
        trafo = 0;
    }
#if defined(USE_FDT)
    if (EQUILIBRIUM and ! HUBBARD_MODEL) {
#if CONTOUR_BASIS != 1
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
#else
#if defined(ZERO_TEMP) and defined(USE_FDT)
    if (is_zero_due_to_FDTs<k2>(i0, w, v, 0, channel)) trafo = -1; // components zero according to FDTs
#endif //ZERO_TEMP
#endif // CONTOUR_BASIS
    }
#endif // EQUILIBRIUM...

#endif // DEBUG_SYMMETRIES
    return trafo;
}

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
int
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_trafo_K3(const int i0, const double w, const double v, const double vp){
    int trafo = 1;
#ifdef DEBUG_SYMMETRIES
    trafo = 0; // compute integrals for all frequency components

    if (!KELDYSH and !ZERO_T and (-v + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K3.frequencies.get_wlower_f() or -vp + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K3.frequencies.get_wlower_f())) {
        trafo = -1;
    }
#else
    const double safety = 1e-5;
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

    if constexpr (!KELDYSH and !ZERO_T and (-v + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K3.frequencies.get_wlower_f() or -vp + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K3.frequencies.get_wlower_f())) {
        trafo = -1;
        //std::cout << "omitted frequencies: " << v << "\t" << vp << std::endl;
        //std::cout << "with limits " << vertex1.avertex().K3.frequencies.get_wlower_f() << std::endl;
    }
#if defined(USE_FDT)
    if (EQUILIBRIUM and ! HUBBARD_MODEL) {
#if CONTOUR_BASIS != 1
        if (i0 == 0 or i0 == 1) trafo = -1; // components can be determined via FDTs, no need to compute it via integration
#else
//#ifdef ZERO_TEMP

        if (is_zero_due_to_FDTs<k3>(i0, w, v, vp, channel)) trafo = -1; // components zero according to FDTs

//#endif //ZERO_TEMP
#endif // CONTOUR_BASIS
    }
#endif //EQUILIBRIUM...

#endif // DEBUG_SYMMETRIES
    return trafo;
}

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
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

template<char channel, typename Q, vertexType symmetry_result, vertexType symmetry_left,
        vertexType symmetry_right, class Bubble_Object>
Q
BubbleFunctionCalculator<channel, Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::bubble_value_prefactor(){
    if constexpr (KELDYSH) return prefactor * (1. / (2. * M_PI * glb_i));
    else                   return prefactor * (1. / (2. * M_PI));
}



// bubble_function using the new class BubbleFunctionCalculator
template <typename Q,
        vertexType symmetry_result,
        vertexType symmetry_left,
        vertexType symmetry_right,
        class Bubble_Object>
void bubble_function(GeneralVertex<Q, symmetry_result>& dgamma,
                                 const GeneralVertex<Q, symmetry_left>& vertex1,
                                 const GeneralVertex<Q, symmetry_right>& vertex2,
                                 const Bubble_Object& Pi,
                                 const char channel){

    if (channel == 'a') {
        BubbleFunctionCalculator<'a', Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi);
        BubbleComputer.perform_computation();
    }
    else if (channel == 'p') {
        BubbleFunctionCalculator<'p', Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi);
        BubbleComputer.perform_computation();
    }
    else if (channel == 't') {
        BubbleFunctionCalculator<'t', Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi);
        BubbleComputer.perform_computation();
    }
    //else {print("Error! Incompatible channel given to bubble_function. Abort"); }

}

/// Overload for bubble_function in case no Bubble object has been initialized yet. ONLY WORKS FOR SIAM!!
template <typename Q,
        vertexType symmetry_result,
        vertexType symmetry_left,
        vertexType symmetry_right>
void bubble_function(GeneralVertex<Q, symmetry_result>& dgamma,
                     const GeneralVertex<Q, symmetry_left>& vertex1,
                     const GeneralVertex<Q, symmetry_right>& vertex2,
                     const Propagator<Q>& G, const Propagator<Q>& S, const char channel, const bool diff){
    Bubble<Q> Pi(G, S, diff);
    bubble_function(dgamma, vertex1, vertex2, Pi, channel);
}

#endif //KELDYSH_MFRG_BUBBLE_FUNCTION_HPP
