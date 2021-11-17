/**
 * Classes / functions for computing bubbles.
 *
 * class Bubble            : Object containing two propagators, with call operator
 * class Integrand_Ki      : Object with call operator providing the integrand for bubble frequency integration
 *                           in diagrammatic class Ki:
 *                           Computing Vertex * Bubble * Vertex, and performing the internal Keldysh summation
 * class Integrand_Ki_diff : Same as Integrand_Ki, for differentiated bubble
 * bubble_function()       : Computing the bubble frequency integral, using the integrand classes and the integrator
 *                           from "integrator.h". Using MPI+OMP parallelization in the external arguments.
 */

#ifndef KELDYSH_MFRG_BUBBLES_H
#define KELDYSH_MFRG_BUBBLES_H

#include <cmath>                            // for using the macro M_PI as pi
#include "symmetries/Keldysh_symmetries.h"  // for independent Keldysh components and utilities
#include "vertex.h"                         // vertex class
#include "selfenergy.h"                     // self-energy class
#include "propagator.h"                     // propagator class
#include "integrator/integrator.h"          // integration routines
#include "utilities/util.h"                 // measuring time, printing text output
#include "utilities/mpi_setup.h"            // mpi parallelization routines
#include "correctionFunctions.h"            // correction terms due to finite integration range
#include "utilities/write_data2file.h"      // write vectors into hdf5 file
#include "grids/momentum_grid.h"            // Momentum grid specific to the 2D Hubbard model
#include "vertex_data.h"


/// Class combining two propagators, either GG or GS+SG
template <typename Q>
class Bubble{
public:
    const Propagator<Q>& g; // Access needed when computing the full bubble
    const Propagator<Q>& s;
    const bool diff;
    /**
     * Constructor:
     * @param propagatorG : first propagator (always a standard one)
     * @param propagatorS : second propagator (standard or single-scale/differentiated, depending on "diff_in")
     * @param diff_in      : whether to compute standard (false) or differentiated (true) bubble
     */
    Bubble(const Propagator<Q>& propagatorG, const Propagator<Q>& propagatorS, const bool diff_in)
        :g(propagatorG), s(propagatorS), diff(diff_in) {};

    /**
     * Call operator:
     * @param iK    : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param v1    : frequency of first propagator
     * @param v2    : frequency of second propagator
     * @param i_in  : internal structure index
     * @return comp : value of the bubble evaluated at (iK, v1, v2)
     */
    auto value(int iK, double v1, double v2, int i_in) const -> Q{
        Q ans;
        if(diff){
            if (KELDYSH){
                switch (iK) {
                    case 3: //AA
                        ans = myconj(g.valsmooth(0, v1, i_in)) * myconj(s.valsmooth(0, v2, i_in)) + myconj(s.valsmooth(0, v1, i_in)) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 6: //AR
                        ans = myconj(g.valsmooth(0, v1, i_in)) * s.valsmooth(0, v2, i_in) + myconj(s.valsmooth(0, v1, i_in)) * g.valsmooth(0, v2, i_in);
                        break;
                    case 7: //AK
                        ans = myconj(g.valsmooth(0, v1, i_in)) * s.valsmooth(1, v2, i_in) + myconj(s.valsmooth(0, v1, i_in)) * g.valsmooth(1, v2, i_in);
                        break;
                    case 9: //RA
                        ans = g.valsmooth(0, v1, i_in) * myconj(s.valsmooth(0, v2, i_in)) + s.valsmooth(0, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 11://KA
                        ans = g.valsmooth(1, v1, i_in) * myconj(s.valsmooth(0, v2, i_in)) + s.valsmooth(1, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 12://RR
                        ans = g.valsmooth(0, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
                        break;
                    case 13://RK
                        ans = g.valsmooth(0, v1, i_in) * s.valsmooth(1, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(1, v2, i_in);
                        break;
                    case 14://KR
                        ans = g.valsmooth(1, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(1, v1, i_in) *  g.valsmooth(0, v2, i_in);
                        break;
                    case 15://KK
                        ans = g.valsmooth(1, v1, i_in) * s.valsmooth(1, v2, i_in) + s.valsmooth(1, v1, i_in) * g.valsmooth(1, v2, i_in);
                        break;
                    default:
                        return 0.;
                }
            }
            else{ // Matsubara
                ans = g.valsmooth(0, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
            }
        }
        else {
            if (KELDYSH){
                switch (iK){ // labelling propagators from top (t: left) to bottom (t: right); a,t: G(v+w/2)G(v-w/2), p: G(w/2-v)G(w/2+v)
                    case 3: //AA
                        ans = myconj(g.valsmooth(0, v1, i_in)) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 6: //AR
                        ans = myconj(g.valsmooth(0, v1, i_in)) * g.valsmooth(0, v2, i_in);
                        break;
                    case 7: //AK
                        ans = myconj(g.valsmooth(0, v1, i_in)) * g.valsmooth(1, v2, i_in);
                        break;
                    case 9: //RA
                        ans = g.valsmooth(0, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 11://KA
                        ans = g.valsmooth(1, v1, i_in) * myconj(g.valsmooth(0, v2, i_in));
                        break;
                    case 12://RR
                        ans = g.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
                        break;
                    case 13://RK
                        ans = g.valsmooth(0, v1, i_in) * g.valsmooth(1, v2, i_in);
                        break;
                    case 14://KR
                        ans =  g.valsmooth(1, v1, i_in) *  g.valsmooth(0, v2, i_in);
                        break;
                    case 15://KK
                        ans =  g.valsmooth(1, v1, i_in) *  g.valsmooth(1, v2, i_in);
                        break;
                    default:
                        return 0.;
                }
            }
            else{ // Matsubara
                ans = g.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
            }
        }
        if (!KELDYSH) {
            if (PARTICLE_HOLE_SYMMETRY) {
                ans *= -1.;     // -1=glb_i^2; needed for particle-hole symmetry in Matsubara (we only save the imaginary part of self-energy and propagators)
            }
        }

        assert(isfinite(ans) == true);
        return ans;
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
        double wa_1, wa_2, wp_1, wp_2, wt_1, wt_2;
        if (KELDYSH || ZERO_T){
            wa_1 = wa_2 = wp_1 = wp_2 = wt_1 = wt_2 = w / 2.;
        }
        else{ // bosonic frequencies have to be rounded to an integer value for finite-temperature Matsubara calculation
            wa_1 = wp_2 = wt_1 = floor2bfreq(w / 2.);
            wa_2 = wp_1 = wt_2 = ceil2bfreq(w / 2.);
        } // TODO(medium): Put this first part into an extra function?
        Q Pival;
        switch (channel) {
            case 'a':
                Pival = value(iK, vpp - wa_1, vpp + wa_2, i_in);    //vppa-1/2wa, vppa+1/2wa for the a-channel
                break;
            case 'p':
                Pival = value(iK, wp_1 + vpp, wp_2 - vpp, i_in);    //wp/2+vppp, wp/2-vppp for the p-channel
                break;
            case 't':
                Pival = value(iK, vpp - wt_1, vpp + wt_2, i_in);    //vppt-1/2wt, vppt+1/2wt for the t-channel
                break;
            default:;
        }
        assert(isfinite(Pival) == true);
        return Pival;
    }
};

template <typename Q>
class PrecalculateBubble{
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

    PrecalculateBubble(const Propagator<Q>& G_in, const Propagator<Q>& S_in,
                       const bool diff_in)
                       :g(G_in), s(S_in), diff(diff_in),
                       Helper_Bubble(G_in, S_in, diff),
                       fermionic_grid('f', 1, G_in.Lambda){
        if (diff) {print("Precalculating a differentiated bubble...", true);}
        else {print("Precalculating a regular bubble...", true);}
        compute_FermionicBubble();
        print("...done.", true);
    }
    auto value(int iK, double w, double vpp, int i_in, char channel) const -> Q;
    auto value_on_fermionic_grid(int iK_bubble, double v1, double v2, int i_in) const -> Q;
    int composite_index(int iK_bubble, int iv1, int iv2, int i_in) const;
};

template <typename Q> auto PrecalculateBubble<Q>::value(int iK, double w, double vpp,
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

template <typename Q> auto PrecalculateBubble<Q>::value_on_fermionic_grid(const int iK_bubble, const double v1, const double v2, const int i_in) const -> Q{
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


template <typename Q> void PrecalculateBubble<Q>::compute_FermionicBubble(){
    if (HUBBARD_MODEL) compute_FermionicBubble_HUBBARD();
    else compute_FermionicBubble_SIAM();
}

template <typename Q> void PrecalculateBubble<Q>::compute_FermionicBubble_HUBBARD(){
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

template <typename Q> void PrecalculateBubble<Q>::compute_FermionicBubble_SIAM(){
    for (int iK_bubble = 0; iK_bubble < glb_number_of_Keldysh_components_bubble; ++iK_bubble) {
        int iK = get_iK_actual(iK_bubble);
        for (int iv1 = 0; iv1 < nFER; ++iv1) {
            for (int iv2 = 0; iv2 < nFER; ++iv2) {
                perform_internal_sum(iK, iv1, iv2);
            }
        }
    }
}

template <typename Q> void PrecalculateBubble<Q>::perform_internal_sum(const int iK, const int iv1, const int iv2){
    double v1 = fermionic_grid.get_ws(iv1);
    double v2 = fermionic_grid.get_ws(iv2);
    for (int i_in = 0; i_in < n_in; ++i_in) {
        FermionicBubble[composite_index(get_iK_bubble(iK), iv1, iv2, i_in)] =
                Helper_Bubble.value(iK, v1, v2, i_in);
    }
}

//TODO(high): Write a flexible method for this in utilities which can handle all types of data (bubbles, vertices, self-energies etc.)
template <typename Q> int PrecalculateBubble<Q>::composite_index(const int iK_bubble, const int iv1, const int iv2, const int i_in) const{
    return iK_bubble*nFER*nFER*n_in + iv1*nFER*n_in + iv2*n_in + i_in;
}

/* Map the 16 Keldysh indices to the 9 non-trivial ones. */
template<typename Q>
int PrecalculateBubble<Q>::get_iK_bubble(const int iK_actual) const {
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
int PrecalculateBubble<Q>::get_iK_actual(const int iK_bubble) const {
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
void PrecalculateBubble<Q>::perform_internal_sum_2D_Hubbard(const int iK, const int iv1, const int iv2,
                                                            Minimal_2D_FFT_Machine& Swave_Bubble_Calculator) {
    double v1 = fermionic_grid.get_ws(iv1);
    double v2 = fermionic_grid.get_ws(iv2);

    vec<comp> values_of_bubble (glb_N_transfer);
    compute_internal_bubble(iK, v1, v2, Swave_Bubble_Calculator, values_of_bubble);

    for (int i_in = 0; i_in < n_in; ++i_in) {
        FermionicBubble[composite_index(get_iK_bubble(iK), iv1, iv2, i_in)] = values_of_bubble[i_in];
    }
}

template<typename Q>
void PrecalculateBubble<Q>::compute_internal_bubble(const int iK, const double v1, const double v2,
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
void PrecalculateBubble<Q>::set_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
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
void PrecalculateBubble<Q>::set_Matsubara_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
                                                      const double v1, const double v2,
                                                      vec<Q>& first_propagator, vec<Q>& second_propagator) const {
    for (int i_in = 0; i_in < glb_N_transfer; ++i_in) { //TODO: Careful! This only works for s-wave. Otherwise n_in > glb_N_transfer!
        first_propagator[i_in]  = g1.valsmooth(0, v1, i_in);
        second_propagator[i_in] = g2.valsmooth(0, v2, i_in);
    }
}

template<typename Q>
void
PrecalculateBubble<Q>::set_Keldysh_propagators(const Propagator<Q>& g1, const Propagator<Q>& g2,
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
// End of Hubbard model specific functions


//Class created for debugging of the Bubbles
template <typename Q>
class IntegrandBubble{
    const Propagator<Q>& g1;
    const Propagator<Q>& g2;
    bool diff;
    double w;
    int iK;
    char channel;

public:
    /**
     * Constructor for the IntegrandBubble
     * @param g1_in     : Propagator object for the lower (right) leg
     * @param g2_in     : Propagator object for the upper (left) leg
     * @param diff_in   : Boolean defining whether bubble is differentiated or not
     * @param w_in      : Transfer frequency at which the bubble should be integrated
     * @param iK_in     : Keldysh index to be taken
     * @param channel_in: Char indicating the channel in which the bubble should be calculated and which determines the frequency transformations
     */
    IntegrandBubble(const Propagator<Q>& g1_in, const Propagator<Q>& g2_in, bool diff_in, double w_in, int iK_in, char channel_in)
            : g1(g1_in), g2(g2_in), diff(diff_in), w(w_in), iK(iK_in), channel(channel_in) {};

    /**
     * Call operator
     * @param vpp : v'', the frequency over which is being integrated
     * @return The value g1(v1)*g2(v2), where v1 and v2 are calculated according to the channel. The components of the
     * propagators taken depend on the Keldysh component
     */
    auto operator() (double vpp) const -> Q {
        Q ans;
        double v1, v2;
        Bubble<Q> Pi(g1, g2, diff);

        switch(channel){
            case 'p':
                v1 = w/2.+vpp;
                v2 = w/2.-vpp;
                break;
            case 'a': case 't':
                v1 = vpp-w/2.;
                v2 = vpp+w/2.;
                break;
            default:
                v1 = 0.;
                v2 = 0.;
                print("Error in IntegrandBubble! Abort.");
                assert(false);
        }
        //Make reference to the Bubble object of the actual code, making this into a useful test of code correctnes and compliance
        return Pi.value(iK, v1, v2, 0)/(2.*M_PI*glb_i);
    }
};

/// Refactoring of the classes Integrand_K1, Integrand_K2, Integrand_K3 into one single class
template <typename Q,
        template <typename> class symmetry_left,
        template <typename> class symmetry_right,
        class Bubble_Object>
class Integrand {
private:
    const GeneralVertex<Q, symmetry_left>& vertex1;
    const GeneralVertex<Q, symmetry_right>& vertex2;
    const Bubble_Object& Pi;
    int i0 = 0;
    const int i2;
    const int i_in;
    const char channel;
    const int iw=0;
    const double w, v = 0., vp = 0.;
    const bool diff;
    const int spin;

    K_class diag_class;

    Q res_l_V_initial, res_r_V_initial, res_l_Vhat_initial, res_r_Vhat_initial; // To be precomputed for K1

    void set_Keldysh_index_i0(int i0_in);
    void precompute_vertices();

    bool case_always_has_to_be_zero() const;

    void compute_vertices(double vpp, Q& res_l_V, Q& res_r_V, Q& res_l_Vhat, Q& res_r_Vhat) const;

public:
    /**
     * Constructor for K1-class:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0 or 1) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     * @param diff_in    : determines whether to compute differentiated or non-differentiated bubble
     */


    /**
     * Constructor for K3-class:
     * Same as for K2 plus additionally
     * @param vp_in      : external fermionic frequency \nu'
     */
    Integrand(const GeneralVertex<Q, symmetry_left>& vertex1_in,
              const GeneralVertex<Q, symmetry_right>& vertex2_in,
              const Bubble_Object& Pi_in,
              int i0_in, int i2_in, const int iw_in, const double w_in, const double v_in, const double vp_in, const int i_in_in,
              const char ch_in, const bool diff_in, const int spin_in, const K_class k_in)
              :vertex1(vertex1_in), vertex2(vertex2_in), Pi(Pi_in),
              i2(i2_in), iw(iw_in), w(w_in), v(v_in), vp(vp_in), i_in(i_in_in), channel(ch_in), diff(diff_in), spin(spin_in), diag_class(k_in){
        set_Keldysh_index_i0(i0_in);
    }

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q;

    void save_integrand() const;
    void save_integrand(const rvec& freqs, const std::string& filename_prefix) const;
    void get_integrand_vals(const rvec& freqs, rvec& integrand_re, rvec& integrand_im, rvec& Pival_re, rvec& Pival_im)  const;
};

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::set_Keldysh_index_i0(const int i0_in) {
    if (KELDYSH){
#ifndef DEBUG_SYMMETRIES
        switch (diag_class) {
            case k1: // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
                switch (channel) {
                    case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
                    case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
                    case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
                    default: ;
                }
                break;
            case k2: // converting index i0_in (0,...,4) into actual Keldysh index i0 (0,...,15)
                switch (channel) {
                    case 'a': i0 = non_zero_Keldysh_K2a[i0_in]; break;
                    case 'p': i0 = non_zero_Keldysh_K2p[i0_in]; break;
                    case 't': i0 = non_zero_Keldysh_K2t[i0_in]; break;
                    default: ;
                }
                break;
            case k3:
                i0 = non_zero_Keldysh_K3[i0_in]; // converting index i0_in (0,...,5) into actual Keldysh index i0 (0,...,15)
                break;
            default: ;
        }
#else
    i0 = i0_in;
#endif
    }
    else{
        i0 = 0;
    }
}

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::precompute_vertices(){
    // For K1 class, left and right vertices do not depend on integration frequency
    // -> precompute them to save time
#ifndef DEBUG_SYMMETRIES
#ifdef KELDYSH_FORMALISM
    std::vector<int> indices = indices_sum(i0, i2, channel);

    VertexInput input_l (indices[0], w, 0., 0., i_in, 0, channel);
    VertexInput input_r (indices[1], w, 0., 0., i_in, 0, channel);
#else
    VertexInput input_l (0, w, 0., 0., i_in, 0, channel);
    VertexInput &input_r = input_l;
#endif
    res_l_V_initial = vertex1[spin].left_same_bare(input_l);
    res_r_V_initial = vertex2[spin].right_same_bare(input_r);
    if (channel == 't') {
        input_l.spin = 1;
        input_r.spin = 1;
        res_l_Vhat_initial = vertex1[spin].left_same_bare(input_l);
        res_r_Vhat_initial = vertex2[spin].right_same_bare(input_r);
    }
#else // DEBUG_SYMMETRIES

#ifdef KELDYSH_FORMALISM
    std::vector<int> indices = indices_sum(i0, i2, channel);

    VertexInput input_l (indices[0], w, 0., 0., i_in, spin, channel);
    VertexInput input_r (indices[1], w, 0., 0., i_in, spin, channel);
#else
    VertexInput input_l (0, w, 0., 0., i_in, spin, channel);
    VertexInput &input_r = input_l;
#endif
    res_l_V_initial = vertex1.left_same_bare(input_l);
    res_r_V_initial = vertex2.right_same_bare(input_r);
    if (channel == 't' and spin == 0) {
        input_l.spin = 1 - input_l.spin; // flip spin 0 <-> 1
        input_r.spin = 1 - input_r.spin; // flip spin 0 <-> 1
        res_l_Vhat_initial = vertex1.left_same_bare(input_l);
        res_r_Vhat_initial = vertex2.right_same_bare(input_r);
    }
    else if (channel == 'a' and spin == 1) {
        input_l.spin = 1 - input_l.spin; // flip spin 0 <-> 1
        input_r.spin = 1 - input_r.spin; // flip spin 0 <-> 1
        res_l_Vhat_initial = vertex1.left_same_bare(input_l);
        res_r_Vhat_initial = vertex2.right_same_bare(input_r);
    }
#endif // DEBUG_SYMMETRIES
}

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
auto Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::operator()(double vpp) const -> Q {
    if (case_always_has_to_be_zero()) {return 0.;}
    Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
    compute_vertices(vpp, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat);

    Q Pival = Pi.value(i2, w, vpp, i_in, channel);
    Q result;
    if ((channel == 't' and spin == 0)
#ifdef DEBUG_SYMMETRIES
    or (channel == 'a' and spin == 1)
#endif
    )
        result = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
    else
        result = res_l_V * Pival * res_r_V;

#ifdef DEBUG_SYMMETRIES
        if (channel == 'p' and spin == 0) {
        result = (res_l_V * Pival * res_r_V +  res_l_Vhat * Pival * res_r_Vhat) * 0.5;
    }
    else if (channel == 'p' and spin == 1) {
        result = (res_l_V * Pival * res_r_Vhat + res_l_Vhat * Pival * res_r_V) * 0.5;
        }
#endif

    return result;
}

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
bool Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::case_always_has_to_be_zero() const {
    bool zero_result = false;
    if (KELDYSH && (MAX_DIAG_CLASS <= 1)){
        if (!diff) {
            switch (channel) {
                case 'a':
                    // only nonzero combinations of \int dvpp Gamma_0 Pi(vpp) Gamma_0
                    if (i0 == 1 && (i2 != 11 && i2 != 13)) {zero_result = true;}
                    if (i0 == 3 &&
                        (i2 != 6 && i2 != 7 && i2 != 9 && i2 != 11 && i2 != 13 && i2 != 14 && i2 != 15))
                        {zero_result = true;}
                    break;
                case 'p':
                    // only nonzero combinations of \int dvpp Gamma_0 Pi(vpp) Gamma_0
                    if (i0 == 1 && (i2 != 7 && i2 != 11)) {zero_result = true;}
                    if (i0 == 5 &&
                        (i2 != 3 && i2 != 7 && i2 != 11 && i2 != 12 && i2 != 13 && i2 != 14 && i2 != 15))
                        {zero_result = true;}
                    break;
                case 't':
                    // only nonzero combinations of \int dvpp Gamma_0 Pi(vpp) Gamma_0
                    if (i0 == 1 && (i2 != 11 && i2 != 13)) {zero_result = true;}
                    if (i0 == 3 &&
                        (i2 != 6 && i2 != 7 && i2 != 9 && i2 != 11 && i2 != 13 && i2 != 14 && i2 != 15))
                        {zero_result = true;}
                    break;
                default:;
            }
        }
    }
    return zero_result;
}

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::compute_vertices(const double vpp,
                                                                                  Q& res_l_V, Q& res_r_V,
                                                                                  Q& res_l_Vhat, Q& res_r_Vhat) const{
    if (MAX_DIAG_CLASS <= 1){
        res_l_V = res_l_V_initial;
        res_r_V = res_r_V_initial;
        res_l_Vhat = res_l_Vhat_initial;
        res_r_Vhat = res_r_Vhat_initial;
    }
    else{
        std::vector<int> indices = indices_sum(i0, i2, channel);
        VertexInput input_l (indices[0], w, v, vpp,  i_in, spin, channel, diag_class, iw);
        VertexInput input_r (indices[1], w, vpp, vp, i_in, spin, channel, diag_class, iw);

        if (diag_class == k1 or diag_class == k2b)
            res_l_V = vertex1.left_same_bare(input_l);
        else
            res_l_V = vertex1.left_diff_bare(input_l);

        if (diag_class == k3 or diag_class == k2b)
            res_r_V = vertex2.right_diff_bare(input_r);
        else
            res_r_V = vertex2.right_same_bare(input_r);

        if ((channel == 't' and spin == 0)
#ifdef DEBUG_SYMMETRIES
        or (channel == 'a' and spin == 1) or (channel == 'p')
#endif
        ){
            if(channel == 't') assert(0 == input_l.spin);
            if(channel == 'a') assert(1 == input_l.spin);
            input_l.spin = 1 - input_l.spin;
            input_r.spin = 1 - input_r.spin;
            if (diag_class == k1 or diag_class == k2b)
                res_l_Vhat = vertex1.left_same_bare(input_l);
            else
                res_l_Vhat = vertex1.left_diff_bare(input_l);

            if (diag_class == k3 or diag_class == k2b)
                res_r_Vhat = vertex2.right_diff_bare(input_r);
            else
                res_r_Vhat = vertex2.right_same_bare(input_r);
        }
    }
}


template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::get_integrand_vals(const rvec& freqs, rvec& integrand_re, rvec& integrand_im, rvec& Pival_re, rvec& Pival_im) const {
    int npoints = freqs.size();
    for (int i=0; i<npoints; ++i) {

        double vpp = freqs[i];


        Q Pival, integrand_value;
        if (not KELDYSH and std::abs(std::abs(w/2)-std::abs(vpp)) < 1e-6) {
            Pival = std::numeric_limits<Q>::infinity();
            integrand_value = std::numeric_limits<Q>::infinity();
        }
        else {
            Pival = Pi.value(i2, w, vpp, i_in, channel);
            integrand_value = (*this)(vpp);
        }
        if constexpr (PARTICLE_HOLE_SYMMETRY && (!KELDYSH)){
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
            Pival_re[i] = Pival;
            Pival_im[i] = 0.;
        }
        else{
            integrand_re[i] = myreal(integrand_value);
            integrand_im[i] = myimag(integrand_value);
            Pival_re[i] = myreal(Pival);
            Pival_im[i] = myimag(Pival);
        }
    }


}

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::save_integrand() const {
    /// Define standard frequency points on which to evaluate the integrand
    int npoints = nBOS;
    if (diag_class == k2) {npoints = 1000;}
    else if (diag_class == k3) {npoints = 100;}

    rvec freqs (npoints);

    for (int i=0; i<npoints; ++i) {
        double wl, wu;
        switch (diag_class) {
            case k1:
                wl = vertex1[spin].avertex().K1_get_wlower() * 2.;
                wu = vertex1[spin].avertex().K1_get_wupper() * 2.;
                break;
            case k2:
                wl = vertex1[spin].avertex().K2_get_wlower_f();
                wu = vertex1[spin].avertex().K2_get_wupper_f();
                break;
            case k3:
                wl = vertex1[spin].avertex().K3_get_wlower_f();
                wu = vertex1[spin].avertex().K3_get_wupper_f();
                break;
            default:;
        }
        double vpp = wl + i * (wu - wl) / (npoints - 1);
        if (diag_class == k1) { vpp = vertex1[spin].avertex().K1_get_freq_w(vpp, i); }
        freqs[i] = vpp;
    }

    save_integrand(freqs, "");

}

template<typename Q, template <typename> class symmetry_left, template <typename> class symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::save_integrand(const rvec& freqs, const std::string& filename_prefix) const {
    /// evaluate the integrand on frequency points in freqs
    int npoints = freqs.size();

    rvec integrand_re (npoints);
    rvec integrand_im (npoints);
    rvec Pival_re (npoints);
    rvec Pival_im (npoints);

    get_integrand_vals(freqs, integrand_re, integrand_im, Pival_re, Pival_im);

    std::string filename = data_dir+filename_prefix+"integrand_K" + std::to_string(diag_class);
    filename += channel;
    filename += "_i0=" + std::to_string(i0)
                + "_i2=" + std::to_string(i2)
                + "_w=" + std::to_string(w);
    if (diag_class == k2) {filename += "_v=" + std::to_string(v);}
    else if (diag_class == k3) {filename += "_vp=" + std::to_string(vp);}
    filename += + ".h5";
    write_h5_rvecs(filename,
                   {"v", "integrand_re", "integrand_im", "Pival_re", "Pival_im"},
                   {freqs, integrand_re, integrand_im, Pival_re, Pival_im});
}




template <typename Q,
        template <typename> class symmetry_result,
        template <typename> class symmetry_left,
        template <typename> class symmetry_right,
        class Bubble_Object>
class BubbleFunctionCalculator{
    private:
    GeneralVertex<Q, symmetry_result>& dgamma;

    const GeneralVertex<Q, symmetry_left>& vertex1;
    const GeneralVertex<Q, symmetry_right>& vertex2;

    const Bubble_Object& Pi;
    const char channel;
    const bool diff = Pi.diff;
    const int spin;

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

    void calculate_value_K1(Q& value, int i0, int i_in, double w);
    void calculate_value_K2(Q& value, int i0, int i_in, double w, double v);
    void calculate_value_K2b(Q& value, int i0, int i_in, double w, double v);
    void calculate_value_K3(Q& value, int i0, int i_in, const int iw, double w, double v, double vp);

    void write_out_results(const vec<Q>& Ordered_result, K_class diag_class);
    void write_out_results_K1(const vec<Q>& K1_ordered_result);
    void write_out_results_K2(const vec<Q>& K2_ordered_result);
    void write_out_results_K2b(const vec<Q>& K2b_ordered_result);
    void write_out_results_K3(const vec<Q>& K3_ordered_result);

    void set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, K_class diag_class);

    void convert_external_MPI_OMP_indices_to_physical_indices_K1(int& iK1, int& i0, int& iw, int& i_in, double& w,
                                                                 int i_mpi, int n_omp, int i_omp);
    void convert_external_MPI_OMP_indices_to_physical_indices_K2(int& iK2, int& i0, int& iw, int& iv, int& i_in,
                                                                 double& w, double& v,
                                                                 int i_mpi, int n_omp, int i_omp);
    void convert_external_MPI_OMP_indices_to_physical_indices_K2b(int& iK2, int& i0, int& iw, int& ivp, int& i_in,
                                                                 double& w, double& vp,
                                                                 int i_mpi, int n_omp, int i_omp);
    void convert_external_MPI_OMP_indices_to_physical_indices_K3(int& iK3, int& i0, int& iw, int& iv, int& ivp, int& i_in,
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
                             const char channel_in, int spin_in)
                             :dgamma(dgamma_in), vertex1(vertex1_in), vertex2(vertex2_in),
                             Pi(Pi_in), channel(channel_in), spin(spin_in){
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

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
Bubble_Object>::set_channel_specific_freq_ranges_and_prefactor() {
// set channel-specific frequency ranges and prefactor (1, 0.5, -1 for a, p, t) for sum over spins.
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

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>::find_vmin_and_vmax() {
    // use std::min/std::max of selfenergy/K1 frequency grids as integration limits
    vmin = std::min(dgamma[spin].avertex().K1.K1_get_wlower(), Pi.g.selfenergy.frequencies.w_lower);
    vmax = std::max(dgamma[spin].avertex().K1.K1_get_wupper(), Pi.g.selfenergy.frequencies.w_upper);
    if (MAX_DIAG_CLASS >= 2){
        // use std::min/std::max of selfenergy/K1/K2 frequency grids as integration limits
        vmin = std::min(vmin, dgamma[spin].avertex().K2.K2_get_wlower_f());
        vmax = std::max(vmax, dgamma[spin].avertex().K2.K2_get_wupper_f());
    }
    if (MAX_DIAG_CLASS >= 3){
        // use std::min/std::max of selfenergy/K1/K2/K3 frequency grids as integration limits
        vmin = std::min(vmin, dgamma[spin].avertex().K3.K3_get_wlower_f());
        vmax = std::max(vmax, dgamma[spin].avertex().K3.K3_get_wupper_f());
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
        template <typename> class symmetry_result,
        template <typename> class symmetry_left,
        template <typename> class symmetry_right,
        class Bubble_Object>
bool
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>::missing_cross_projection() {
    switch (channel) {
        case 'a':
            if (!vertex1[spin].pvertex().calculated_crossprojections ||
                !vertex1[spin].tvertex().calculated_crossprojections ||
                !vertex2[spin].pvertex().calculated_crossprojections ||
                !vertex2[spin].tvertex().calculated_crossprojections)
            {
                return true;
            }
            break;
        case 'p':
            if (!vertex1[spin].avertex().calculated_crossprojections ||
                !vertex1[spin].tvertex().calculated_crossprojections ||
                !vertex2[spin].avertex().calculated_crossprojections ||
                !vertex2[spin].tvertex().calculated_crossprojections)
            {
                return true;
            }
            break;
        case 't':
            if (!vertex1[spin].avertex().calculated_crossprojections ||
                !vertex1[spin].pvertex().calculated_crossprojections ||
                !vertex2[spin].avertex().calculated_crossprojections ||
                !vertex2[spin].pvertex().calculated_crossprojections)
            {
                return true;
            }
            break;
        default:;
    }
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::perform_computation(){

    double t_start = get_time();
    if (MAX_DIAG_CLASS >= 0) {
        calculate_bubble_function(k1);
        tK1 = get_time() - t_start;
        //print("K1", channel, " done, ");
        //get_time(t_start);
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

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
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

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
Q
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::get_value(const int i_mpi, const int i_omp, const int n_omp, const K_class diag_class){
    Q value = 0.;
    int iK1, iK2, i0, iw, iv, ivp, i_in;
    double w, v, vp;
    int trafo;
    switch (diag_class) {
        case k1:
            convert_external_MPI_OMP_indices_to_physical_indices_K1(iK1, i0, iw, i_in, w,
                                                                    i_mpi, n_omp, i_omp);
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
#else
            trafo = get_trafo_K1(i0, w);
#endif
            if (trafo == 0) {calculate_value_K1(value, i0, i_in, w); }
            break;
        case k2:
            convert_external_MPI_OMP_indices_to_physical_indices_K2(iK2, i0, iw, iv, i_in, w, v,
                                                                    i_mpi, n_omp, i_omp);
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
#else
            trafo = get_trafo_K2(i0, w, v);
#endif
            if (trafo == 0) {calculate_value_K2(value, i0, i_in, w, v); }
            break;
#ifdef DEBUG_SYMMETRIES
        case k2b:
            convert_external_MPI_OMP_indices_to_physical_indices_K2b(iK2, i0, iw, ivp, i_in, w, vp,
                                                                    i_mpi, n_omp, i_omp);
            trafo = 0; // compute integrals for all frequency components
#else
            //trafo = get_trafo_K2(i0, w, v);
#endif
            if (trafo == 0) {calculate_value_K2b(value, i0, i_in, w, vp); }
            break;
        case k3:
            convert_external_MPI_OMP_indices_to_physical_indices_K3(iK2, i0, iw, iv, ivp, i_in, w, v, vp,
                                                                    i_mpi, n_omp, i_omp);
#ifdef DEBUG_SYMMETRIES
            trafo = 0; // compute integrals for all frequency components
#else
            trafo = get_trafo_K3(i0, w, v, vp);
#endif
            if (trafo == 0) {
                if (channel == 'a') {calculate_value_K3(value, i0, i_in, iw, w, v, vp); } // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
                if (channel == 'p') {calculate_value_K3(value, i0, i_in, iv, w, v, vp); } // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
                if (channel == 't') {calculate_value_K3(value, i0, i_in, ivp,w, v, vp); } // for 2D interpolation of K3 we need to know the index of the constant bosonic frequency w_r (r = channel of the bubble)
            }
            break;
        default:;
    }
    return value;
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::calculate_value_K1(Q& value, const int i0, const int i_in, const double w){
    for (int i2 : glb_non_zero_Keldysh_bubble) {
        Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>
                integrand_K1(vertex1, vertex2, Pi, i0, i2, 0, w, 0., 0., i_in, channel, diff, spin, k1);
        if (KELDYSH){
            value += bubble_value_prefactor() * integrator<Q>(integrand_K1, vmin, vmax, -w / 2., w / 2., Delta);
        }
        else{
            if (ZERO_T){
                //if (std::abs(w) < 1e-2) integrand_K1.save_integrand();
                value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,0>(integrand_K1, vmin, vmax, std::abs(w/2), {}, Delta, false);
            }
            else{
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP) // TODO(high): Figure out type problems in matsubarasum
                //int interval_correction =  (int)(signFlipCorrection_MF(w)/(2*M_PI*glb_T) + 0.1); // if interval_correction=-1, then the integrand is symmetric around v=-M_PI*glb_T
                value += bubble_value_prefactor()*(2*M_PI) * glb_T * matsubarasum<Q>(integrand_K1, Nmin, Nmax);
#endif
            }
        }
        value += bubble_value_prefactor() *
                asymp_corrections_bubble(k1, vertex1, vertex2, Pi.g, vmin, vmax,
                                         w, 0., 0., i0, i2, i_in, channel, diff, spin);

    }
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_value_K2(Q& value, const int i0, const int i_in, const double w, const double v){
    if (false and vertex2[spin].Ir()) {value = 0.;} // right part of multi-loop contribution does not contribute to K2 class
                                                    // don't do this; this is a bug for the Parquet solver
    else {
        for (int i2 : glb_non_zero_Keldysh_bubble) { // TODO(medium): Add a sum over form factors here, as well. (the integrand class will need to take another index.)
            Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>
                    integrand_K2(vertex1, vertex2, Pi, i0, i2, 0, w, v, 0., i_in, channel, diff, spin, k2);
            if (KELDYSH){
                value += bubble_value_prefactor() * integrator<Q>(integrand_K2, vmin, vmax, -w / 2., w / 2., Delta);
            }
            else{
                if (ZERO_T){
                    //if (std::abs(w) < 1e-2 or std::abs(w+0.54) < 1e-2) integrand_K2.save_integrand();
                    value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,3>(integrand_K2, vmin, vmax, std::abs(w/2), {v, v+w, v-w}, Delta, false);
                    //value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,0>(integrand_K2, vmin, vmax, std::abs(w/2), {}, Delta); // TODO(high): Remove?!
                }
                else{
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP) // TODO(high): Figure out type problems in matsubarasum
                    //int interval_correction =  (int)(signFlipCorrection_MF(w)/(2*M_PI*glb_T) + 0.1);
                    // if interval_correction=-1, then the integrand is symmetric around v=-M_PI*glb_T
                    value += bubble_value_prefactor()*(2*M_PI) * glb_T * matsubarasum<Q>(integrand_K2, Nmin, Nmax);
#endif
                }
            }
            value += bubble_value_prefactor() *
                     asymp_corrections_bubble(k2, vertex1, vertex2, Pi.g,
                                              vmin, vmax, w, v, 0., i0, i2, i_in, channel, diff, spin);

        }
    }
}

#ifdef DEBUG_SYMMETRIES
template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_value_K2b(Q& value, const int i0, const int i_in, const double w, const double vp){
    if (false and vertex1[spin].Ir()) {value = 0.;} // left part of multi-loop contribution does not contribute to K2 class
                                                    // don't do this; this is a bug for the Parquet solver
    else {
        for (int i2 : glb_non_zero_Keldysh_bubble) { // TODO(medium): Add a sum over form factors here, as well. (the integrand class will need to take another index.)
            Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>
                    integrand_K2b(vertex1, vertex2, Pi, i0, i2, 0, w, 0., vp, i_in, channel, diff, spin, k2b);
            if (KELDYSH){
                value += bubble_value_prefactor() * integrator<Q>(integrand_K2b, vmin, vmax, -w / 2., w / 2., Delta);
            }
            else{
                if (ZERO_T){
                    //if (std::abs(w) < 1e-2 or std::abs(w+0.54) < 1e-2) integrand_K2.save_integrand();
                    value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,3>(integrand_K2b, vmin, vmax, std::abs(w/2), {vp, vp+w, vp-w}, Delta, false);
                    //value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,0>(integrand_K2, vmin, vmax, std::abs(w/2), {}, Delta); // TODO(high): Remove?!
                }
                else{
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP) // TODO(high): Figure out type problems in matsubarasum
                    //int interval_correction =  (int)( signFlipCorrection_MF(w)/(2*M_PI*glb_T) + 0.1);
                    // if interval_correction=-1, then the integrand is symmetric around v=-M_PI*glb_T
                    value += bubble_value_prefactor()*(2*M_PI) * glb_T * matsubarasum<Q>(integrand_K2b, Nmin, Nmax);
#endif
                }
            }
            value += bubble_value_prefactor() *
                     asymp_corrections_bubble(k2b, vertex1, vertex2, Pi.g,
                                              vmin, vmax, w, vp, 0., i0, i2, i_in, channel, diff, spin);

        }
    }
}
#endif

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::calculate_value_K3(Q& value, const int i0, const int i_in,
                                                const int iw, const double w, const double v, const double vp){
    for (int i2 : glb_non_zero_Keldysh_bubble) { // TODO(medium): Add two sums over form factors here, as well. (the integrand class will need to take another index.)
        // initialize the integrand object and perform frequency integration
        Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>
                integrand_K3(vertex1, vertex2, Pi, i0, i2, iw, w, v, vp, i_in, channel, diff, spin, k3);
        if (KELDYSH){
            value += bubble_value_prefactor() * integrator<Q>(integrand_K3, vmin, vmax, -w / 2., w / 2., Delta);
        }
        else{
            if (ZERO_T){
                value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,6>(integrand_K3, vmin, vmax, std::abs(w/2), {v, vp, w-vp, w+vp, w-v, abs(w)+abs(v)}, Delta, false);
                //value += bubble_value_prefactor() * integrator_Matsubara_T0<Q,0>(integrand_K3, vmin, vmax, std::abs(w/2), {}, Delta); // TODO(high): Remove?!
            }
            else{
#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP) // TODO(high): Figure out type problems in matsubarasum
                //int interval_correction =  (int)(signFlipCorrection_MF(w)/(2*M_PI*glb_T) + 0.1);
                // if interval_correction=-1, then the integrand is symmetric around v=-M_PI*glb_T
                value += bubble_value_prefactor()*(2*M_PI) * glb_T * matsubarasum<Q>(integrand_K3, Nmin, Nmax);
#endif
            }
        }
        value += bubble_value_prefactor() *
                asymp_corrections_bubble(k3, vertex1, vertex2, Pi.g,
                                         vmin, vmax, w, v, vp, i0, i2, i_in, channel, diff, spin);

    }
}


template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results(const vec<Q>& Ordered_result, const K_class diag_class){
    dgamma[spin].half1().initializeInterpol();     // initialize Interpolator with the symmetry-reduced sector of the vertex to retrieve all remaining entries
                                                /// TODO: does cubic interpolation overshoot at the edges of the symmetry-reduced sector?
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
    dgamma[spin].half1().set_initializedInterpol(false);      // above initialization of the Interpolator is with the symmetry-reduced sector only (rest = zero)
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::write_out_results_K1(const vec<Q>& K1_ordered_result){
    switch (channel) {
        case 'a':
            dgamma[spin].avertex().K1.add_vec(K1_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].avertex().enforce_freqsymmetriesK1(dgamma[spin].avertex());
#endif
            break;
        case 'p':
            dgamma[spin].pvertex().K1.add_vec(K1_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].pvertex().enforce_freqsymmetriesK1(dgamma[spin].pvertex());
#endif
            break;
        case 't':
            dgamma[spin].tvertex().K1.add_vec(K1_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].tvertex().enforce_freqsymmetriesK1(dgamma[spin].tvertex());
#endif
            break;
        default: ;
    }
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K2(const vec<Q>& K2_ordered_result){
    assert( K2_ordered_result.size() == nK_K2*nBOS2*nFER2*n_in);
    switch (channel) {
        case 'a':
            dgamma[spin].avertex().K2.add_vec(K2_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].avertex().enforce_freqsymmetriesK2(dgamma[spin].avertex());
#endif
            break;
        case 'p':
            dgamma[spin].pvertex().K2.add_vec(K2_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].pvertex().enforce_freqsymmetriesK2(dgamma[spin].pvertex());
#endif
            break;
        case 't':
            dgamma[spin].tvertex().K2.add_vec(K2_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].tvertex().enforce_freqsymmetriesK2(dgamma[spin].tvertex());
#endif
            break;
        default: ;
    }
}

#ifdef DEBUG_SYMMETRIES
template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K2b(const vec<Q>& K2b_ordered_result){
    assert( K2b_ordered_result.size() == nK_K2*nBOS2*nFER2*n_in);
    switch (channel) {
        case 'a':
            dgamma[spin].avertex().K2b.add_vec(K2b_ordered_result);
#ifndef DEBUG_SYMMETRIES
            //dgamma[spin].avertex().enforce_freqsymmetriesK2(dgamma[spin].avertex());
#endif
            break;
        case 'p':
            dgamma[spin].pvertex().K2b.add_vec(K2b_ordered_result);
#ifndef DEBUG_SYMMETRIES
            //dgamma[spin].pvertex().enforce_freqsymmetriesK2(dgamma[spin].pvertex());
#endif
            break;
        case 't':
            dgamma[spin].tvertex().K2b.add_vec(K2b_ordered_result);
#ifndef DEBUG_SYMMETRIES
            //dgamma[spin].tvertex().enforce_freqsymmetriesK2(dgamma[spin].tvertex());
#endif
            break;
        default: ;
    }
}
#endif

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::write_out_results_K3(const vec<Q>& K3_ordered_result){
    switch (channel) {
        case 'a':
            dgamma[spin].avertex().K3.add_vec(K3_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].avertex().enforce_freqsymmetriesK3(dgamma[spin].avertex());
#endif
            break;
        case 'p':
            dgamma[spin].pvertex().K3.add_vec(K3_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].pvertex().enforce_freqsymmetriesK3(dgamma[spin].pvertex());
#endif
            break;
        case 't':
            dgamma[spin].tvertex().K3.add_vec(K3_ordered_result);
#ifndef DEBUG_SYMMETRIES
            dgamma[spin].tvertex().enforce_freqsymmetriesK3(dgamma[spin].tvertex());
#endif
            break;
        default: ;
    }
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::set_external_arguments_for_parallelization(int& n_mpi, int& n_omp, const K_class diag_class){
    switch (diag_class) {
        case k1:
            n_mpi = nK_K1;        // set external arguments for MPI-parallelization (# of tasks distributed via MPI)
            n_omp = nw1_w * n_in_K1; // set external arguments for OMP-parallelization (# of tasks per MPI-task distributed via OMP)
            break;
        case k2:
        case k2b:
            n_mpi = nK_K2 * nw2_w;
            n_omp = nw2_v * n_in_K2;
            break;
        case k3:
            n_mpi = nK_K3 * nw3_w;
            n_omp = nw3_v * nw3_v_p * n_in_K3;
            break;
        default: ;
    }
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K1(int& iK1, int& i0, int& iw, int& i_in, double& w,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK1 = i_mpi * n_omp + i_omp;
    i0 = iK1/(nw1_w*n_in_K1);                              // exterior Keldysh indices of the bubble
    iw = iK1/(n_in_K1) - i0*nw1_w;                         // frequency index
    i_in = iK1 - i0*nw1_w*n_in_K1 - iw*n_in_K1;            // internal index
    if (channel == 'a') dgamma[spin].avertex().K1.K1_get_freq_w(w, iw);           // frequency acc. to frequency index
    if (channel == 'p') dgamma[spin].pvertex().K1.K1_get_freq_w(w, iw);           // frequency acc. to frequency index
    if (channel == 't') dgamma[spin].tvertex().K1.K1_get_freq_w(w, iw);           // frequency acc. to frequency index
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2(int& iK2, int& i0, int& iw, int& iv, int& i_in,
                                                                                double& w, double& v,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;
    i0 = iK2 / (nw2_w * nw2_v * n_in_K2);
    iw = iK2 / (nw2_v * n_in_K2) - i0 * nw2_w;
    iv = iK2 / n_in_K2 - iw * nw2_v - i0 * nw2_w * nw2_v;
    i_in = iK2 - iv * n_in_K2 - iw * nw2_v * n_in_K2 - i0 * nw2_w * nw2_v * n_in_K2;
    if (channel == 'a') dgamma[spin].avertex().K2.K2_get_freqs_w(w, v, iw, iv);
    if (channel == 'p') dgamma[spin].pvertex().K2.K2_get_freqs_w(w, v, iw, iv);
    if (channel == 't') dgamma[spin].tvertex().K2.K2_get_freqs_w(w, v, iw, iv);
}

#ifdef DEBUG_SYMMETRIES
template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K2b(int& iK2, int& i0, int& iw, int& ivp, int& i_in,
                                                                                double& w, double& vp,
                                                                                const int i_mpi, const int n_omp, const int i_omp){
    iK2 = i_mpi * n_omp + i_omp;
    i0 = iK2 / (nw2_w * nw2_v * n_in_K2);
    iw = iK2 / (nw2_v * n_in_K2) - i0 * nw2_w;
    ivp= iK2 / n_in_K2 - iw * nw2_v - i0 * nw2_w * nw2_v;
    i_in = iK2 - ivp * n_in_K2 - iw * nw2_v * n_in_K2 - i0 * nw2_w * nw2_v * n_in_K2;
    if (channel == 'a') dgamma[spin].avertex().K2b.K2_get_freqs_w(w, vp, iw, ivp);
    if (channel == 'p') dgamma[spin].pvertex().K2b.K2_get_freqs_w(w, vp, iw, ivp);
    if (channel == 't') dgamma[spin].tvertex().K2b.K2_get_freqs_w(w, vp, iw, ivp);
}
#endif

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
void
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
        Bubble_Object>::convert_external_MPI_OMP_indices_to_physical_indices_K3(int& iK3, int& i0, int& iw, int& iv, int& ivp, int& i_in,
                                                                                     double& w, double& v, double& vp,
                                                                                     const int i_mpi, const int n_omp, const int i_omp){
    iK3 = i_mpi * n_omp + i_omp;
    i0 = iK3/(nw3_w * nw3_v * nw3_v_p * n_in_K3);
    iw = iK3/(nw3_v * nw3_v_p * n_in_K3) - i0*nw3_w;
    iv = iK3/(nw3_v * n_in_K3) - i0*nw3_w*nw3_v - iw*nw3_v;
    ivp =iK3/(n_in_K3) - i0*nw3_w*nw3_v*nw3_v_p - iw*nw3_v*nw3_v_p - iv*nw3_v_p;
    i_in = iK3 - i0*nw3_w*nw3_v*nw3_v_p*n_in_K3 - iw*nw3_v*nw3_v_p*n_in_K3 - iv*nw3_v_p*n_in_K3 - ivp*n_in_K3;
    if (channel == 'a') dgamma[spin].avertex().K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, 'a');
    if (channel == 'p') dgamma[spin].pvertex().K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, 'p');
    if (channel == 't') dgamma[spin].tvertex().K3.K3_get_freqs_w(w, v, vp, iw, iv, ivp, 't');
}


template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
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

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
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
    return trafo;
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
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
    return trafo;
}

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
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

template<typename Q, template <typename> class symmetry_result, template <typename> class symmetry_left,
        template <typename> class symmetry_right, class Bubble_Object>
Q
BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right,
                Bubble_Object>::bubble_value_prefactor(){
    if constexpr (KELDYSH) return prefactor * (1. / (2. * M_PI * glb_i));
    else                   return prefactor * (1. / (2. * M_PI));
}


// bubble_function using the new class BubbleFunctionCalculator
template <typename Q,
        template <typename> class symmetry_result,
        template <typename> class symmetry_left,
        template <typename> class symmetry_right,
        class Bubble_Object>
void bubble_function(GeneralVertex<Q, symmetry_result>& dgamma,
                                 const GeneralVertex<Q, symmetry_left>& vertex1,
                                 const GeneralVertex<Q, symmetry_right>& vertex2,
                                 const Bubble_Object& Pi,
                                 const char channel){
    for (int i_spin = 0; i_spin < n_spin; i_spin++) {
        BubbleFunctionCalculator<Q, symmetry_result, symmetry_left, symmetry_right, Bubble_Object>
                BubbleComputer (dgamma, vertex1, vertex2, Pi, channel, i_spin);
        if (channel == 'a' or channel == 'p' or channel == 't') BubbleComputer.perform_computation();
        //else {print("Error! Incompatible channel given to bubble_function. Abort"); }
    }
}

/// Overload for bubble_function in case no Bubble object has been initialized yet. ONLY WORKS FOR SIAM!!
template <typename Q,
        template <typename> class symmetry_result,
        template <typename> class symmetry_left,
        template <typename> class symmetry_right>
void bubble_function(GeneralVertex<Q, symmetry_result>& dgamma,
                     const GeneralVertex<Q, symmetry_left>& vertex1,
                     const GeneralVertex<Q, symmetry_right>& vertex2,
                     const Propagator<Q>& G, const Propagator<Q>& S, const char channel, const bool diff){
    Bubble<Q> Pi(G, S, diff);
    bubble_function(dgamma, vertex1, vertex2, Pi, channel);
}

#endif //KELDYSH_MFRG_BUBBLES_H
