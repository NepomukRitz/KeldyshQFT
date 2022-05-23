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
        typename vertexType_result,
        typename vertexType_left,
        typename vertexType_right,
        class Bubble_Object>
class BubbleFunctionCalculator{
    using value_type = vec<Q>;
    private:
    vertexType_result& dgamma;

    const vertexType_left vertex1;  /// THIS IS A COPY; needed for symmetry-expansion
    const vertexType_right vertex2; /// THIS IS A COPY; needed for symmetry-expansion

    const Bubble_Object& Pi;
    const bool diff = Pi.diff;

    const double Delta = (Pi.g.Lambda + glb_Gamma) / 2.; // hybridization (needed for proper splitting of the integration domain)

    int nw1_w = 0, nw2_w = 0, nw2_v = 0, nw3_w = 0, nw3_v = 0, nw3_v_p = 0;
    Q prefactor = 1.;

    int mpi_size = mpi_world_size(); // number of mpi processes
    int mpi_rank = mpi_world_rank(); // number of the current mpi process
    std::array<std::size_t,my_defs::K1::rank> dimsK1;
    std::array<std::size_t,my_defs::K2::rank> dimsK2;
    std::array<std::size_t,my_defs::K3::rank> dimsK3;
    std::vector<int> indepKeldyshComponents_K1;
    std::vector<int> indepKeldyshComponents_K2;
    std::vector<int> indepKeldyshComponents_K3;

    double vmin = 0, vmax = 0;
    int Nmin, Nmax; // Matsubara indices for minimal and maximal frequency. Only needed for finite-temperature Matsubara calculations!

    double tK1 = 0, tK2 = 0, tK3 = 0;

    void check_presence_of_symmetry_related_contributions();
    void set_channel_specific_freq_ranges_and_prefactor();
    void find_vmin_and_vmax();

    bool missing_cross_projection(); // Needed for the Hubbard model.

    template<K_class diag_class> void calculate_bubble_function();
    template<K_class diag_class> value_type get_value(int i_mpi, int i_omp, int n_omp, int n_vectorization);
    template<K_class diag_class,int spin> void calculate_value(value_type &value, int i0, int i_in, int iw, double w, double v, double vp);

    void write_out_results(const vec<Q>& Ordered_result, K_class diag_class);
    void write_out_results_K1(const vec<Q>& K1_ordered_result);
    void write_out_results_K2(const vec<Q>& K2_ordered_result);
    void write_out_results_K2b(const vec<Q>& K2b_ordered_result);
    void write_out_results_K3(const vec<Q>& K3_ordered_result);

    void set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, int& n_vectorization, K_class diag_class);

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

    const bool store_integrand_for_PT = false; // set to true if the integrand for the fully retarded up-down component at zero frequency shall be saved

    public:
    void perform_computation();

    BubbleFunctionCalculator(vertexType_result& dgamma_in,
                             const vertexType_left& vertex1_in,
                             const vertexType_right& vertex2_in,
                             const Bubble_Object& Pi_in)
                             :dgamma(dgamma_in), vertex1(vertex1_in), vertex2(vertex2_in),
                             Pi(Pi_in){
#if not  DEBUG_SYMMETRIES
        //check_presence_of_symmetry_related_contributions();
#endif
        set_channel_specific_freq_ranges_and_prefactor();
        find_vmin_and_vmax();

        // For Hubbard model computations, make sure that the internal structures of the vertices are parametrized correctly.
        if (HUBBARD_MODEL && (MAX_DIAG_CLASS > 1) && missing_cross_projection()) { // Cross projected parts are only needed for the Hubbard model in K2 and K3.
            utils::print("Error! Needed crossprojection still has to be computed. Abort.");
            assert(false);
        }
#if SWITCH_SUM_N_INTEGRAL
        vertex1.template symmetry_expand<channel,true>();
        vertex2.template symmetry_expand<channel,false>();
#endif
        /// TODO(high): Figure out computations which need gamma_a_uu = gamma_a_ud - gamma_t_ud in a t-bubble,
        ///  i.e. CP_to_t(gamma_a_uu) = CP_to_t(gamma_a_ud) - CP_to_a(gamma_t_ud).
        ///  The integrand will need vertex AND vertex_initial to have access to cross-projected parts and non-crossprojected parts.


        dimsK1 = vertex1.get_rvertex(channel).K1.get_dims();
        dimsK2 = MAX_DIAG_CLASS > 1 ? vertex1.get_rvertex(channel).K2.get_dims() : std::array<std::size_t, my_defs::K2::rank>();
        dimsK3 = MAX_DIAG_CLASS > 2 ? vertex1.get_rvertex(channel).K3.get_dims() : std::array<std::size_t, my_defs::K3::rank>();
        if constexpr (VECTORIZED_INTEGRATION == 1) {
            /// vectorization over Keldysh indices
            dimsK1[my_defs::K1::keldysh] = 1;
            dimsK2[my_defs::K2::keldysh] = 1;
            dimsK3[my_defs::K3::keldysh] = 1;
        }

        indepKeldyshComponents_K1 = channel == 'a' ? non_zero_Keldysh_K1a : (channel == 'p'? non_zero_Keldysh_K1p : non_zero_Keldysh_K1t);
        indepKeldyshComponents_K2 = channel == 'a' ? non_zero_Keldysh_K2a : (channel == 'p'? non_zero_Keldysh_K2p : non_zero_Keldysh_K2t);
        indepKeldyshComponents_K3 = non_zero_Keldysh_K3;

    }
};

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
Bubble_Object>::set_channel_specific_freq_ranges_and_prefactor() {
// set channel-specific frequency ranges and prefactor (1, 1, -1 for a, p, t) for sum over spins.

    nw1_w = dgamma.get_rvertex(channel).K1.get_dims()[my_defs::K1::omega];
    nw2_w = dgamma.get_rvertex(channel).K2.get_dims()[my_defs::K2::omega];
    nw2_v = dgamma.get_rvertex(channel).K2.get_dims()[my_defs::K2::nu];
    nw3_w = dgamma.get_rvertex(channel).K3.get_dims()[my_defs::K3::omega];
    nw3_v = dgamma.get_rvertex(channel).K3.get_dims()[my_defs::K3::nu];
    nw3_v_p=dgamma.get_rvertex(channel).K3.get_dims()[my_defs::K3::nup];

    prefactor = channel == 't' ? -1. : 1.;

}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right, Bubble_Object>::find_vmin_and_vmax() {
    // use std::min/std::max of selfenergy/K1 frequency grids as integration limits
    if (HUBBARD_MODEL){ // In the HM we have a larger frequency box for the SE and want to limit us to the range of bosonic vertex frequencies.
        vmin = dgamma.avertex().K1.frequencies.get_wupper_b();
        vmax = dgamma.avertex().K1.frequencies.get_wupper_b();
    }
    else{
        vmin =-Delta * 10.; // std::min(dgamma.avertex().K1.frequencies.get_wupper_b(), Pi.g.selfenergy.Sigma.frequencies.primary_grid.w_lower);
        vmax = Delta * 10.; // std::max(dgamma.avertex().K1.frequencies.get_wupper_b(), Pi.g.selfenergy.Sigma.frequencies.primary_grid.w_upper);
    }

    if constexpr(MAX_DIAG_CLASS >= 2){
        // use std::min/std::max of selfenergy/K1/K2 frequency grids as integration limits
        //vmin = std::min(vmin, dgamma.avertex().K2.frequencies.get_wlower_f());
        //vmax = std::max(vmax, dgamma.avertex().K2.frequencies.get_wupper_f());
    }
    if constexpr(MAX_DIAG_CLASS >= 3){
        // use std::min/std::max of selfenergy/K1/K2/K3 frequency grids as integration limits
        //vmin = std::min(vmin, dgamma.avertex().K3.frequencies.get_wlower_f());
        //vmax = std::max(vmax, dgamma.avertex().K3.frequencies.get_wupper_f());
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
        typename vertexType_result,
        typename vertexType_left,
        typename vertexType_right,
        class Bubble_Object>
bool
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right, Bubble_Object>::missing_cross_projection() {
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

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
                Bubble_Object>::perform_computation(){

    double t_start = utils::get_time();
    if constexpr(MAX_DIAG_CLASS >= 0) {
        calculate_bubble_function<k1>();
        tK1 = utils::get_time() - t_start;
        //utils::print("K1", channel, " done, ");
        //utils::get_time(t_start);
    }
    if constexpr(MAX_DIAG_CLASS >= 2) {
        t_start = utils::get_time();
        calculate_bubble_function<k2>();
        tK2 = utils::get_time() - t_start;
        //utils::print("K2", channel, " done, ");
        //utils::get_time(t_start);

#if DEBUG_SYMMETRIES
        t_start = utils::get_time();
        calculate_bubble_function<k2b>();
        tK2 = utils::get_time() - t_start;
        //utils::print("K2b", channel, " done, ");
        //utils::get_time(t_start);
#endif
    }
    if constexpr(MAX_DIAG_CLASS >= 3) {
        t_start = utils::get_time();
        calculate_bubble_function<k3>();
        tK3 = utils::get_time() - t_start;
        //utils::print("K3", channel, " done, ");
        //utils::get_time(t_start);
    }

}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
                template<K_class diag_class>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::calculate_bubble_function(){
    if (diag_class < k1 || diag_class > k3){utils::print("Incompatible diagrammatic class! Abort."); assert(false); return;}

    int n_mpi, n_omp, n_vectorization;
    set_external_arguments_for_parallelization(n_mpi, n_omp, n_vectorization, diag_class);


    // initialize buffer into which each MPI process writes their results
    vec<Q> Buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp * n_vectorization);
    vertex1.initializeInterpol();
    vertex2.initializeInterpol();

    // start for-loop over external arguments, using MPI and OMP
    int iterator = 0;
    for (int i_mpi = 0; i_mpi < n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for schedule(dynamic)
            for (int i_omp = 0; i_omp < n_omp; ++i_omp) {
                value_type value = get_value<diag_class>(i_mpi, i_omp, n_omp, n_vectorization);
                for (int k = 0; k < n_vectorization; k++) {
                    Buffer[(iterator * n_omp + i_omp) * n_vectorization + k] = value[k]; // write result of integration into MPI buffer
                }
            }
            ++iterator;
        }
    }
    // collect+combine results from different MPI processes, reorder them appropriately
    vec<Q> Result = mpi_initialize_result<Q> (n_mpi, n_omp * n_vectorization);
    mpi_collect(Buffer, Result, n_mpi, n_omp * n_vectorization);
    vec<Q> Ordered_result = mpi_reorder_result(Result, n_mpi, n_omp * n_vectorization);

    write_out_results(Ordered_result, diag_class);

    vertex1.set_initializedInterpol(false);
    vertex2.set_initializedInterpol(false);

}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
template<K_class diag_class>
auto
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::get_value(const int i_mpi, const int i_omp, const int n_omp, const int n_vectorization) -> value_type{
    value_type value(n_vectorization);
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
#if DEBUG_SYMMETRIES
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

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
template<K_class k, int spin>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::calculate_value(value_type& value, const int i0, const int i_in, const int iw,
                                           const double w, const double v, const double vp){
#if VECTORIZED_INTEGRATION==1
    using Integrand_class = Integrand<k, channel, spin, Q, vertexType_left, vertexType_right, Bubble_Object,Eigen::Matrix<Q,4,4>>;
#else
    using Integrand_class = Integrand<k, channel, spin, Q, vertexType_left, vertexType_right, Bubble_Object,Q>;
#endif
    using integrand_valtype = std::result_of_t<Integrand_class(double)>;
    integrand_valtype integration_result = myzero<integrand_valtype>();

#if not SWITCH_SUM_N_INTEGRAL
    static_assert(n_spin == 1, "SWITCH_SUM_N_INTEGRAL not ready for DEBUG_SYMMETRIES.");
    for (int i2 : glb_non_zero_Keldysh_bubble) {
        int n_spin_sum = 1;                  // number of summands in spin sum (=1 in the a channel)
        if ((channel == 't' and spin == 0) or (channel == 'a' and spin == 1)) n_spin_sum = 3;  // in the t channel, spin sum includes three terms
        for (int i_spin=0; i_spin < n_spin_sum; ++i_spin) {
#else
    int i2 = 0;
    int i_spin = 0;
#endif
            // initialize the integrand object and perform frequency integration
            Integrand_class integrand(vertex1, vertex2, Pi, i0, i2, iw, w, v, vp, i_in, i_spin, diff);

            if (store_integrand_for_PT){ // for PT-Calculations: Store the integrand for the fully retarded up-down component at zero frequency:
#if SWITCH_SUM_N_INTEGRAL
                assert(false); // Cannot do it this way if the integration is vectorized over the internal Keldysh sum
#endif
                if (channel == 'a' and i0 == 7 and i_spin == 0 and w == 0 and v == 0 and vp == 0) integrand.save_integrand();
            }


            if constexpr(ZERO_T) {
                switch (k) {
                    case k1:
                        integration_result += bubble_value_prefactor() *
                                 integrator_Matsubara_T0(integrand, vmin, vmax, std::abs(w / 2),
                                                               {0.}, Delta, true);

                break;
            case k2:
                integration_result += bubble_value_prefactor() * integrator_Matsubara_T0(integrand, vmin, vmax, std::abs(w / 2), {0., v, v + w, v - w}, Delta, true);
                break;
            case k3:
                integration_result += bubble_value_prefactor() * integrator_Matsubara_T0(integrand, vmin, vmax, std::abs(w / 2), {0, v, vp, w - vp, w + vp, w - v, w + v}, Delta, true);
                break;
            case k2b:
                integration_result += bubble_value_prefactor() * integrator_Matsubara_T0(integrand, vmin, vmax, std::abs(w / 2), {0, vp, vp + w, vp - w}, Delta, true);
                break;
            default:
                break;
        }
    } else {
        if constexpr(KELDYSH) {
            integration_result += bubble_value_prefactor() * integrator(integrand, vmin, vmax, -w / 2., w / 2., Delta, true);
        } else {
            int interval_correction = signFlipCorrection_MF_int(w);
            int W = (int) (w / (2 * M_PI * glb_T) + 0.1 * sgn(w));
            double vmin_temp = (-POSINTRANGE - std::abs(W / 2) + interval_correction) * 2 * M_PI * glb_T;
            double vmax_temp = (POSINTRANGE - 1 + std::abs(W / 2)) * 2 * M_PI * glb_T;
            // if interval_correction=-1, then the integrand is symmetric_full around v=-M_PI*glb_T

            integration_result = bubble_value_prefactor() * (2 * M_PI) * glb_T *
                    matsubarasum<Q>(integrand, -POSINTRANGE - std::abs(W / 2) + interval_correction, POSINTRANGE - 1 + std::abs(W / 2));

            integration_result +=
                    bubble_value_prefactor() * asymp_corrections_bubble<channel>(k, vertex1, vertex2, Pi.g,
                                                                                 vmin_temp, vmax_temp, w, v, vp, i0, i2,
                                                                                 i_in, diff, spin);
        }
    }

#if not SWITCH_SUM_N_INTEGRAL
    }
#else
    //for (int i2: glb_non_zero_Keldysh_bubble)
    {
#endif
        // asymptotic corrections include spin sum
        if constexpr (!HUBBARD_MODEL and (ZERO_T or KELDYSH) and false) {
            integration_result +=
                    bubble_value_prefactor() * asymp_corrections_bubble<channel>(k, vertex1, vertex2, Pi.g,
                                                                                 vmin, vmax, w, v, vp, i0, i2,
                                                                                 i_in, diff, spin);
        }

    }

    /// write integration_result into value
    if constexpr(VECTORIZED_INTEGRATION == 1) {
        // for vector-/matrix-valued result:
        if constexpr(DEBUG_SYMMETRIES) {
            // if DEBUG_SYMMETRIES is true, we compute and store ALL components
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    value[rotate_Keldysh_matrix<channel, true>(i * 4 + j)] = integration_result(i, j);
                }
            }
        }
        else {
            // if DEBUG_SYMMETRIES is false, we compute and store symmetry-reduced components (given in indepKeldyshComponents_Ki)
            std::vector<int> &indepKeldyshComponents = (k == k1 ? indepKeldyshComponents_K1 : (k == k3 ? indepKeldyshComponents_K3 : indepKeldyshComponents_K2));
            const int size = indepKeldyshComponents.size();
            for (int i = 0; i < size; i++) {
                int left, right;
                get_i0_left_right<channel>(indepKeldyshComponents[i], left, right);
                const Q val_temp = integration_result(left, right);
                value[i] = val_temp;
            }
        }
    } else {
        // for scalar result:
        value[0] = integration_result;
    }
}


template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::write_out_results(const vec<Q>& Ordered_result, const K_class diag_class){
    switch (diag_class) {
        case k1:
            write_out_results_K1(Ordered_result);
            break;
        case k2:
            write_out_results_K2(Ordered_result);
            break;
#if DEBUG_SYMMETRIES
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

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
                Bubble_Object>::write_out_results_K1(const vec<Q>& K1_ordered_result){
    dgamma.get_rvertex(channel).K1.add_vec(K1_ordered_result);
    if constexpr(not DEBUG_SYMMETRIES) {
        dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
        if (not HUBBARD_MODEL) dgamma.get_rvertex(channel).enforce_freqsymmetriesK1(dgamma.get_rvertex(channel));
    }

}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::write_out_results_K2(const vec<Q>& K2_ordered_result){
    //assert( K2_ordered_result.size() == nK_K2*n_spin*nBOS2*nFER2*n_in);
    dgamma.get_rvertex(channel).K2.add_vec(K2_ordered_result);
    if constexpr(not DEBUG_SYMMETRIES) {
        dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
        if (not HUBBARD_MODEL) dgamma.get_rvertex(channel).enforce_freqsymmetriesK2(dgamma.get_rvertex(channel));
    }
#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and USE_FDT
    compute_components_through_FDTs(dgamma.half1(), dgamma.half1(), dgamma.half1(), channel);
#endif
}

#if DEBUG_SYMMETRIES
template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::write_out_results_K2b(const vec<Q>& K2b_ordered_result){
    //assert( K2b_ordered_result.size() == nK_K2*n_spin*nBOS2*nFER2*n_in);
    dgamma.get_rvertex(channel).K2b.add_vec(K2b_ordered_result);
}
#endif

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::write_out_results_K3(const vec<Q>& K3_ordered_result){
    dgamma.get_rvertex(channel).K3.add_vec(K3_ordered_result);
    if constexpr(not DEBUG_SYMMETRIES) {
        dgamma.initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
        if (not HUBBARD_MODEL) dgamma.get_rvertex(channel).enforce_freqsymmetriesK3(dgamma.get_rvertex(channel));
    }
#if defined(EQUILIBRIUM) and not defined(HUBBARD_MODEL) and USE_FDT
    compute_components_through_FDTs(dgamma.half1(), dgamma.half1(), dgamma.half1(), channel);
#endif
}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, int& n_vectorization, const K_class diag_class){

    const int nK_K1 = dimsK1[my_defs::K1::keldysh];
    const int nK_K2 = dimsK2[my_defs::K2::keldysh];
    const int nK_K3 = dimsK3[my_defs::K3::keldysh];

    switch (diag_class) {
        case k1:
            n_mpi = nK_K1 * n_spin;        // set external arguments for MPI-parallelization (# of tasks distributed via MPI)
            n_omp = nw1_w * n_in_K1; // set external arguments for OMP-parallelization (# of tasks per MPI-task distributed via OMP)
            n_vectorization = VECTORIZED_INTEGRATION == 1 ? vertex1.get_rvertex(channel).K1.get_dims()[my_defs::K1::keldysh] : 1;
            break;
        case k2:
        case k2b:
            n_mpi = nK_K2 * n_spin * nw2_w;
            n_omp = nw2_v * n_in_K2;
            n_vectorization = VECTORIZED_INTEGRATION == 1 ? vertex1.get_rvertex(channel).K2.get_dims()[my_defs::K2::keldysh] : 1;
            break;
        case k3:
            n_mpi = nK_K3 * n_spin * nw3_w;
            n_omp = nw3_v * nw3_v_p * n_in_K3;
            n_vectorization = VECTORIZED_INTEGRATION == 1 ? vertex1.get_rvertex(channel).K3.get_dims()[my_defs::K3::keldysh] : 1;
            break;
        default: ;
    }
}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K1(int& iK1, int& i0, int& ispin, int& iw, int& i_in, double& w,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK1 = i_mpi * n_omp + i_omp;

    my_defs::K1::index_type idx;
    getMultIndex<rank_K1>(idx, iK1, dimsK1);
    i0       = (int) idx[my_defs::K1::keldysh];
    ispin    = (int) idx[my_defs::K1::spin];
    iw       = (int) idx[my_defs::K1::omega];
    i_in     = (int) idx[my_defs::K1::internal];
   //getMultIndex<4,int,int,int,int>(ispin, iw, i0, i_in, iK1, vertex1.avertex().K1.get_dims());

    if (channel == 'a') dgamma.avertex().K1.frequencies.get_freqs_w(w, iw);           // frequency acc. to frequency index
    if (channel == 'p') dgamma.pvertex().K1.frequencies.get_freqs_w(w, iw);           // frequency acc. to frequency index
    if (channel == 't') dgamma.tvertex().K1.frequencies.get_freqs_w(w, iw);           // frequency acc. to frequency index
}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2(int& iK2, int& i0, int& ispin,  int& iw, int& iv, int& i_in,
                                                                                double& w, double& v,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;

    my_defs::K2::index_type idx;
    getMultIndex<rank_K2>(idx, iK2, dimsK2);
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

#if DEBUG_SYMMETRIES
template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2b(int& iK2, int& i0, int& ispin,  int& iw, int& ivp, int& i_in,
                                                                                double& w, double& vp,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;

    my_defs::K2::index_type idx;
    getMultIndex<rank_K2>(idx, iK2, dimsK2);
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

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K3(int& iK3, int& i0, int& ispin,  int& iw, int& iv, int& ivp, int& i_in,
                                                                                     double& w, double& v, double& vp,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK3 = i_mpi * n_omp + i_omp;

    my_defs::K3::index_type idx;
    getMultIndex<rank_K3>(idx, iK3, dimsK3);
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


template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
int
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
                Bubble_Object>::get_trafo_K1(const int i0, const double w){
    int trafo = 1;

    if constexpr(DEBUG_SYMMETRIES) {
        trafo = 0; // compute integrals for all frequency components
    }
    else {
        if constexpr(VECTORIZED_INTEGRATION) {
            // Make sure that the frequency point does not belong to the symmetry-reduced sector for all relevant Keldysh components
            // otherwise we have to compute that point
            const double safety = 1e-10;
            int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
            const int number_of_indep_Components = dgamma.get_rvertex(channel).K1.get_dims()[my_defs::K1::keldysh];
            for (int i0_temp = 0; i0_temp < number_of_indep_Components; i0_temp++) {
                if (dgamma.get_rvertex(channel).freq_transformations.K1[i0_temp][sign_w] == 0) trafo = 0;
            }

        }
        else {
            const double safety = 1e-10;
            int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
            trafo = dgamma.get_rvertex(channel).freq_transformations.K1[i0][sign_w];

#if CONTOUR_BASIS == 1 and defined(ZERO_TEMP) and USE_FDT
            if (is_zero_due_to_FDTs<k1>(i0, w, 0, 0, channel)) trafo = -1; // components zero according to FDTs
#endif // CONTOUR_BASIS
        } // VECTORIZED_INTEGRATION

    } // DEBUG_SYMMETRIES
    return trafo;
}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
int
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::get_trafo_K2(const int i0, const double w, const double v){
    int trafo = 1;
    if constexpr (DEBUG_SYMMETRIES) {
        trafo = 0; // compute integrals for all frequency components
        if (!KELDYSH and !ZERO_T and -v + signFlipCorrection_MF(w) * 0.5 < vertex1.avertex().K2.frequencies.get_wlower_f()) {
            trafo = -1;
        }
    }
    else {
        if constexpr (VECTORIZED_INTEGRATION) {
            // Make sure that the frequency point does not belong to the symmetry-reduced sector for any relevant Keldysh component
            // otherwise we have to compute that point via quadrature
            const double safety = 1e-10;
            int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
            int sign_v = sign_index<double>(v - safety); // safety to ensure that w=0 gets sign_w=-1
            const int number_of_indep_Components = dgamma.get_rvertex(channel).K2.get_dims()[my_defs::K2::keldysh];

            for (int i0_temp = 0; i0_temp < number_of_indep_Components; i0_temp++) {
                if (dgamma.get_rvertex(channel).freq_transformations.K2[i0_temp][sign_w * 2 + sign_v] == 0) trafo = 0;
            }

        }
        else{
            const double safety = 1e-10;
            int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
            int sign_v = sign_index<double>(v - safety); // safety to ensure that w=0 gets sign_w=-1
            trafo = dgamma.get_rvertex(channel).freq_transformations.K2[i0][sign_w*2 + sign_v];

        if (!KELDYSH and !ZERO_T and -v + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K2.frequencies.get_wlower_f()) {
            trafo = 0;
        }
    #if USE_FDT
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
    #if defined(ZERO_TEMP) and USE_FDT
        if (is_zero_due_to_FDTs<k2>(i0, w, v, 0, channel)) trafo = -1; // components zero according to FDTs
    #endif //ZERO_TEMP
    #endif // CONTOUR_BASIS
        }
    #endif // EQUILIBRIUM...
        } // VECTORIZED_INTEGRATION
    } // DEBUG_SYMMETRIES
    return trafo;
}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
int
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
        Bubble_Object>::get_trafo_K3(const int i0, const double w, const double v, const double vp){
    int trafo = 1;
    if constexpr(DEBUG_SYMMETRIES) {
        trafo = 0; // compute integrals for all frequency components

        if (!KELDYSH and !ZERO_T and
            (-v + signFlipCorrection_MF(w) * 0.5 < vertex1.avertex().K3.frequencies.get_wlower_f() or
             -vp + signFlipCorrection_MF(w) * 0.5 < vertex1.avertex().K3.frequencies.get_wlower_f())) {
            trafo = -1;
        }
    }
    else {
        if constexpr(VECTORIZED_INTEGRATION) {
            // Make sure that the frequency point does not belong to the symmetry-reduced sector for all relevant Keldysh components
            // otherwise we have to compute that point
            const double safety = 1e-10;
            int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
            int sign_f = sign_index(v + vp - safety);
            int sign_fp = sign_index(v - vp - safety);
            const int number_of_indep_Components = dgamma.get_rvertex(channel).K3.get_dims()[my_defs::K3::keldysh];

            for (int i0_temp = 0; i0_temp < number_of_indep_Components; i0_temp++) {
                if (dgamma.get_rvertex(channel).freq_transformations.K3[i0_temp][sign_w * 4 + sign_f * 2 + sign_fp] ==
                    0)
                    trafo = 0;
            }

        }
        else {
            const double safety = 1e-5;
            int sign_w = sign_index<double>(w - safety); // safety to ensure that w=0 gets sign_w=-1
            int sign_f = sign_index(v + vp - safety);
            int sign_fp= sign_index(v - vp - safety);
            trafo = dgamma.get_rvertex(channel).freq_transformations.K3[i0][sign_w * 4 + sign_f * 2 + sign_fp];


            if (!KELDYSH and !ZERO_T and (-v + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K3.frequencies.get_wlower_f() or -vp + signFlipCorrection_MF(w)*0.5 < vertex1.avertex().K3.frequencies.get_wlower_f())) {
                trafo = -1;
                //std::cout << "omitted frequencies: " << v << "\t" << vp << std::endl;
                //std::cout << "with limits " << vertex1.avertex().K3.frequencies.get_wlower_f() << std::endl;
            }
#if USE_FDT
            if (EQUILIBRIUM and ! HUBBARD_MODEL) {
#if CONTOUR_BASIS != 1
                if (i0 == 0 or i0 == 1) trafo = -1; // components can be determined via FDTs, no need to compute it via integration
#else
#ifdef ZERO_TEMP
                if (is_zero_due_to_FDTs<k3>(i0, w, v, vp, channel)) trafo = -1; // components zero according to FDTs
#endif //ZERO_TEMP
#endif // CONTOUR_BASIS
            }
#endif //EQUILIBRIUM...

        } // VECTORIZED_INTEGRATION
    } // DEBUG_SYMMETRIES
    return trafo;
}

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
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

template<char channel, typename Q, typename vertexType_result, typename vertexType_left,
        typename vertexType_right, class Bubble_Object>
Q
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right,
                Bubble_Object>::bubble_value_prefactor(){
    if constexpr (KELDYSH) return prefactor * (1. / (2. * M_PI * glb_i));
    else                   return prefactor * (1. / (2. * M_PI));
}

template<char channel, typename Q,
        typename vertexType_result,
        typename vertexType_left,
        typename vertexType_right,
        class Bubble_Object>
void
BubbleFunctionCalculator<channel, Q, vertexType_result, vertexType_left, vertexType_right, Bubble_Object>::check_presence_of_symmetry_related_contributions() {
    bool vertex1_is_bare = false;
    bool vertex2_is_bare = false;
    if ((vertex1.avertex().max_norm() < 1e-18)
        && (vertex1.pvertex().max_norm() < 1e-18)
        && (vertex1.tvertex().max_norm() < 1e-18)) {vertex1_is_bare = true;}
    if ((vertex2.avertex().max_norm() < 1e-18)
        && (vertex2.pvertex().max_norm() < 1e-18)
        && (vertex2.tvertex().max_norm() < 1e-18)) {vertex2_is_bare = true;}

    if (channel == 'a'){ // There must be a non-vanishing contribution in the t-channel for both vertices
        if (not vertex1_is_bare) assert(vertex1.tvertex().max_norm() > 1e-18);
        if (not vertex2_is_bare) assert(vertex2.tvertex().max_norm() > 1e-18);
    }
    if (channel == 't'){ // There must be a non-vanishing contribution in the a-channel for both vertices
        if (not vertex1_is_bare) assert(vertex1.avertex().max_norm() > 1e-18);
        if (not vertex2_is_bare) assert(vertex2.avertex().max_norm() > 1e-18);
    }
}


// bubble_function using the new class BubbleFunctionCalculator
template <
        typename vertexType_result,
        typename vertexType_left,
        typename vertexType_right,
        class Bubble_Object>
void bubble_function(vertexType_result& dgamma,
                                 const vertexType_left& vertex1,
                                 const vertexType_right& vertex2,
                                 const Bubble_Object& Pi,
                                 const char channel){
    using Q = typename vertexType_result::base_type;
    if (channel == 'a') {
        BubbleFunctionCalculator<'a', Q, vertexType_result, vertexType_left, vertexType_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi);
        BubbleComputer.perform_computation();
    }
    else if (channel == 'p') {
        BubbleFunctionCalculator<'p', Q, vertexType_result, vertexType_left, vertexType_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi);
        BubbleComputer.perform_computation();
    }
    else if (channel == 't') {
        BubbleFunctionCalculator<'t', Q, vertexType_result, vertexType_left, vertexType_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi);
        BubbleComputer.perform_computation();
    }
    //else {utils::print("Error! Incompatible channel given to bubble_function. Abort"); }

}

/// Overload for bubble_function in case no Bubble object has been initialized yet. ONLY WORKS FOR SIAM!!
template <typename Q,
        typename vertexType_result,
        typename vertexType_left,
        typename vertexType_right>
void bubble_function(vertexType_result& dgamma,
                     const vertexType_left& vertex1,
                     const vertexType_right& vertex2,
                     const Propagator<Q>& G, const Propagator<Q>& S, const char channel, const bool diff){
    Bubble<Q> Pi(G, S, diff);
    bubble_function(dgamma, vertex1, vertex2, Pi, channel);
}

#endif //KELDYSH_MFRG_BUBBLE_FUNCTION_HPP
