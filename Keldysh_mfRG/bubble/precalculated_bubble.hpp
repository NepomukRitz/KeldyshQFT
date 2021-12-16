#ifndef KELDYSH_MFRG_PRECALCULATED_BUBBLE_HPP
#define KELDYSH_MFRG_PRECALCULATED_BUBBLE_HPP

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

template <typename Q>
class PrecalculatedBubble{
    const Bubble<Q> Helper_Bubble;

    void compute_FermionicBubble();
    void compute_FermionicBubble_HUBBARD();
    void compute_FermionicBubble_SIAM();
    void perform_internal_sum(int iK, int iv1, int iv2);

    int get_iK_bubble(int iK_actual) const;
    int get_iK_actual(int iK_bubble) const;

    // Hubbard model specific functions
    void perform_internal_sum_2D_Hubbard(int iK, int iv1, int iv2,
                                         Minimal_2D_FFT_Machine& Swave_Bubble_Calculator);
    void compute_internal_bubble(int iK, double v1, double v2,
                                 Minimal_2D_FFT_Machine& Swave_Bubble_Calculator, vec<comp>& values_of_bubble) const;
    void set_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                         int iK, double v1, double v2,
                         vec<Q>& first_propagator, vec<Q>& second_propagator) const;
    void set_Matsubara_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                                   double v1, double v2,
                                   vec<Q>& first_propagator, vec<Q>& second_propagator) const;
    void set_Keldysh_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                                 int iK, double v1, double v2,
                                 vec<Q>& first_propagator, vec<Q>& second_propagator) const;

public:
    const Propagator<Q>& g; // Access needed when computing the full bubble
    const Propagator<Q>& s;
    const bool diff;
    FrequencyGrid fermionic_grid;
    vec<Q> FermionicBubble = vec<Q> (glb_number_of_Keldysh_components_bubble * nFER * nFER * n_in); // 9 non-zero Keldysh components

    PrecalculatedBubble(const Propagator<Q>& G_in, const Propagator<Q>& S_in,
                        const bool diff_in)
                       :g(G_in), s(S_in), diff(diff_in),
                       Helper_Bubble(G_in, S_in, diff),
                       fermionic_grid(G_in.selfenergy.frequencies){
        if (diff) {print("Precalculating a differentiated bubble...", true);}
        else {print("Precalculating a regular bubble...", true);}
        compute_FermionicBubble();
        print("...done.", true);
    }
    auto value(int iK, double w, double vpp, int i_in, char channel) const -> Q;
    auto value_on_fermionic_grid(int iK_bubble, double v1, double v2, int i_in) const -> Q;
    int composite_index(int iK_bubble, int iv1, int iv2, int i_in) const;
};

template <typename Q> auto PrecalculatedBubble<Q>::value(int iK, double w, double vpp,
                                                         int i_in, char channel) const -> Q {
    if (KELDYSH && ((iK == 0) || (iK == 1) || (iK == 2) || (iK == 4) || (iK == 5) || (iK == 8) || (iK == 10))){
        return 0.;
    } // Catch trivial Keldysh indices
    Q Pival;
    switch (channel)
    {
        case 'a':
            Pival = value_on_fermionic_grid(get_iK_bubble(iK), vpp - w / 2., vpp + w / 2., i_in);    //vppa-1/2wa, vppa+1/2wa for the a-channel
            break;
        case 'p':
            Pival = value_on_fermionic_grid(get_iK_bubble(iK), w / 2. + vpp, w / 2. - vpp, i_in);    //wp/2+vppp, wp/2-vppp for the p-channel
            break;
        case 't':
            Pival = value_on_fermionic_grid(get_iK_bubble(iK), vpp - w / 2., vpp + w / 2., i_in);    //vppt-1/2wt, vppt+1/2wt for the t-channel
            break;
        default:;
    }
    return Pival;
}

template <typename Q> auto PrecalculatedBubble<Q>::value_on_fermionic_grid(const int iK_bubble, const double v1, const double v2, const int i_in) const -> Q{
    if (    std::abs(v1) + inter_tol < fermionic_grid.w_upper
            && std::abs(v2) + inter_tol < fermionic_grid.w_upper) {

        Q result = interpolate_lin2D<Q>(v1, v2, fermionic_grid, fermionic_grid,
                                    [&](int i, int j) -> Q {return FermionicBubble[composite_index(iK_bubble, i, j, i_in)];});
        return result;

    }
    else {
        //std::cout << "Out of interpolation tolerance! \n";
        return 0.;
    }
}

template <typename Q> void PrecalculatedBubble<Q>::compute_FermionicBubble(){
    if (HUBBARD_MODEL) compute_FermionicBubble_HUBBARD();
    else compute_FermionicBubble_SIAM();
}

template <typename Q> void PrecalculatedBubble<Q>::compute_FermionicBubble_HUBBARD(){
    double starting_time = get_time();
    std::vector<Minimal_2D_FFT_Machine> FFT_Machinery(omp_get_max_threads());
    double end_time = get_time();
    double time_diff = (end_time - starting_time); // time given in seconds
    //std::cout << "Time for FFT initialization = " << time_diff << " s." << "\n";

    for (int iK_bubble = 0; iK_bubble < glb_number_of_Keldysh_components_bubble; ++iK_bubble) {
        int iK = get_iK_actual(iK_bubble);
        //std::cout << "Now calculating iK = " << iK << "\n";
#pragma omp parallel for schedule(dynamic) default(none) shared(FFT_Machinery, iK)
        for (int iv = 0; iv < nFER * nFER; ++iv) {
            const int iv1 = iv / nFER; // integer division, always rounds down
            const int iv2 = iv - iv1 * nFER;
            //if (iv2 == 0) {std::cout << "Now calculating iK = " << iK << ", iv1 = " << iv1 << "\n";}
            perform_internal_sum_2D_Hubbard(iK, iv1, iv2, FFT_Machinery[omp_get_thread_num()]);
        }
    }
}

template <typename Q> void PrecalculatedBubble<Q>::compute_FermionicBubble_SIAM(){
    for (int iK_bubble = 0; iK_bubble < glb_number_of_Keldysh_components_bubble; ++iK_bubble) {
        int iK = get_iK_actual(iK_bubble);
        for (int iv1 = 0; iv1 < nFER; ++iv1) {
            for (int iv2 = 0; iv2 < nFER; ++iv2) {
                perform_internal_sum(iK, iv1, iv2);
            }
        }
    }
}

template <typename Q> void PrecalculatedBubble<Q>::perform_internal_sum(const int iK, const int iv1, const int iv2){
    double v1 = fermionic_grid.get_ws(iv1);
    double v2 = fermionic_grid.get_ws(iv2);
    for (int i_in = 0; i_in < n_in; ++i_in) {
        FermionicBubble[composite_index(get_iK_bubble(iK), iv1, iv2, i_in)] =
                Helper_Bubble.value(iK, v1, v2, i_in);
    }
}

//TODO(high): Write a flexible method for this in utilities which can handle all types of data (bubbles, vertices, self-energies etc.)
template <typename Q> int PrecalculatedBubble<Q>::composite_index(const int iK_bubble, const int iv1, const int iv2, const int i_in) const{
    return iK_bubble*nFER*nFER*n_in + iv1*nFER*n_in + iv2*n_in + i_in;
}

/* Map the 16 Keldysh indices to the 9 non-trivial ones. */
template<typename Q>
int PrecalculatedBubble<Q>::get_iK_bubble(const int iK_actual) const {
    int iK_bubble = 0;
    if (KELDYSH){
        switch (iK_actual) {
            case  3: iK_bubble = 0; break;
            case  6: iK_bubble = 1; break;
            case  7: iK_bubble = 2; break;
            case  9: iK_bubble = 3; break;
            case 11: iK_bubble = 4; break;
            case 12: iK_bubble = 5; break;
            case 13: iK_bubble = 6; break;
            case 14: iK_bubble = 7; break;
            case 15: iK_bubble = 8; break;
            default:
                print("ERROR! Trivial Keldysh index not filtered out yet! Abort."); assert(false);
        }
    }
    return iK_bubble;
}

/* Map the 9 non-trivial Keldysh indices for the bubble to the 16 original ones. */
template<typename Q>
int PrecalculatedBubble<Q>::get_iK_actual(const int iK_bubble) const {
    int iK_actual = 0;
    if (KELDYSH){
        switch (iK_bubble) {
            case 0: iK_actual =  3; break;
            case 1: iK_actual =  6; break;
            case 2: iK_actual =  7; break;
            case 3: iK_actual =  9; break;
            case 4: iK_actual = 11; break;
            case 5: iK_actual = 12; break;
            case 6: iK_actual = 13; break;
            case 7: iK_actual = 14; break;
            case 8: iK_actual = 15; break;
            default:
                print("ERROR! Number of nine non-trivial Keldysh-indices exceeded! Abort."); assert(false);
        }
    }
    return iK_actual;
}

// Hubbard model specific functions
template<typename Q>
void PrecalculatedBubble<Q>::perform_internal_sum_2D_Hubbard(const int iK, const int iv1, const int iv2,
                                                             Minimal_2D_FFT_Machine& Swave_Bubble_Calculator) {
    const double v1 = fermionic_grid.get_ws(iv1);
    const double v2 = fermionic_grid.get_ws(iv2);

    vec<comp> values_of_bubble (n_in);
    compute_internal_bubble(iK, v1, v2, Swave_Bubble_Calculator, values_of_bubble);

    for (int i_in = 0; i_in < n_in; ++i_in) {
        FermionicBubble[composite_index(get_iK_bubble(iK), iv1, iv2, i_in)] = values_of_bubble[i_in];
    }
}

template<typename Q>
void PrecalculatedBubble<Q>::compute_internal_bubble(const int iK, const double v1, const double v2,
                                                     Minimal_2D_FFT_Machine& Swave_Bubble_Calculator,
                                                     vec<comp>& values_of_bubble) const {
    vec<comp> first_propagator (glb_N_transfer); // input for FFT
    vec<comp> second_propagator (glb_N_transfer); // input for FFT

    if (diff){
        set_propagators(g, s, iK, v1, v2, first_propagator, second_propagator);
        values_of_bubble = Swave_Bubble_Calculator.compute_swave_bubble(first_propagator, second_propagator);

        set_propagators(s, g, iK, v1, v2, first_propagator, second_propagator);
        values_of_bubble += Swave_Bubble_Calculator.compute_swave_bubble(first_propagator, second_propagator);
    }
    else{
        set_propagators(g, g, iK, v1, v2, first_propagator, second_propagator);
        values_of_bubble = Swave_Bubble_Calculator.compute_swave_bubble(first_propagator, second_propagator);
    }
}

template<typename Q>
void PrecalculatedBubble<Q>::set_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                                             const int iK, const double v1, const double v2,
                                             vec<Q>& first_propagator, vec<Q>& second_propagator) const {
    if (KELDYSH){
        set_Keldysh_propagators(g1, g2, iK, v1, v2, first_propagator, second_propagator);
    }
    else{
        set_Matsubara_propagators(g1, g2, v1, v2, first_propagator, second_propagator);
    }
}

template<typename Q>
void PrecalculatedBubble<Q>::set_Matsubara_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                                                       const double v1, const double v2,
                                                       vec<Q>& first_propagator, vec<Q>& second_propagator) const {
    for (int i_in = 0; i_in < glb_N_transfer; ++i_in) { //TODO: Careful! This only works for s-wave. Otherwise n_in > glb_N_transfer!
        first_propagator[i_in]  = g1.valsmooth(0, v1, i_in);
        second_propagator[i_in] = g2.valsmooth(0, v2, i_in);
    }
}

template<typename Q>
void
PrecalculatedBubble<Q>::set_Keldysh_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                                                const int iK, const double v1, const double v2,
                                                vec<Q>& first_propagator, vec<Q>& second_propagator) const {
    for (int i_in = 0; i_in < glb_N_transfer; ++i_in) { //TODO: Careful! This only works for s-wave. Otherwise n_in > glb_N_transfer!
        switch (iK) {
            case 3: //AA
                first_propagator[i_in]  = myconj(g1.valsmooth(0, v1, i_in));
                second_propagator[i_in] = myconj(g2.valsmooth(0, v2, i_in));
                break;
            case 6: //AR
                first_propagator[i_in]  = myconj(g1.valsmooth(0, v1, i_in));
                second_propagator[i_in] = g2.valsmooth(0, v2, i_in);
                break;
            case 7: //AK
                first_propagator[i_in]  = myconj(g1.valsmooth(0, v1, i_in));
                second_propagator[i_in] = g2.valsmooth(1, v2, i_in);
                break;
            case 9: //RA
                first_propagator[i_in]  = g1.valsmooth(0, v1, i_in);
                second_propagator[i_in] = myconj(g2.valsmooth(0, v2, i_in));
                break;
            case 11://KA
                first_propagator[i_in]  = g1.valsmooth(1, v1, i_in);
                second_propagator[i_in] = myconj(g2.valsmooth(0, v2, i_in));
                break;
            case 12://RR
                first_propagator[i_in]  = g1.valsmooth(0, v1, i_in);
                second_propagator[i_in] = g2.valsmooth(0, v2, i_in);
                break;
            case 13://RK
                first_propagator[i_in]  = g1.valsmooth(0, v1, i_in);
                second_propagator[i_in] = g2.valsmooth(1, v2, i_in);
                break;
            case 14://KR
                first_propagator[i_in]  = g1.valsmooth(1, v1, i_in);
                second_propagator[i_in] = g2.valsmooth(0, v2, i_in);
                break;
            case 15://KK
                first_propagator[i_in]  = g1.valsmooth(1, v1, i_in);
                second_propagator[i_in] = g2.valsmooth(1, v2, i_in);
                break;
            default:
                first_propagator[i_in]  = 0.;
                second_propagator[i_in] = 0.;
        }
    }
}



#endif //KELDYSH_MFRG_PRECALCULATED_BUBBLE_HPP
