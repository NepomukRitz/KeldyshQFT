#ifndef KELDYSH_MFRG_INTEGRAND_HPP
#define KELDYSH_MFRG_INTEGRAND_HPP

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

/// TODO: write move constructor (necessary for move in PAID integrator)
/// TODO: write vectorized variant of integrand

/// Possible (unit-)tests
/// for operations: move, call
/// [MISSING] copy integrand, perform move and compare with copy
/// [IMPLEMENTED in test_PT_state()] compute bubble and compare with known result
/// check combination of vertices and propagators and internal sum over indices
/// [IMPLEMENTED in integrand_tets/...] load vertices and selfenergies from file and output integrands to hdf5 file
/// [MISSING] output integrand without and with internal sum, then compare the output (are they identical?)


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
        symmetryType symmetry_left,
        symmetryType symmetry_right,
        class Bubble_Object>
class Integrand {
private:
    const GeneralVertex<Q, symmetry_left>& vertex1;
    const GeneralVertex<Q, symmetry_right>& vertex2;
    const Bubble_Object& Pi;
    int i0 = 0;
    const int i2;
    const int i_in;
    const int i_spin;
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
              int i0_in, int i2_in, const int spin_in, const int iw_in, const double w_in, const double v_in, const double vp_in, const int i_in_in,
              const int i_spin_in, const char ch_in, const bool diff_in, const K_class k_in)
              :vertex1(vertex1_in), vertex2(vertex2_in), Pi(Pi_in),
              i2(i2_in), spin(spin_in), iw(iw_in), w(w_in), v(v_in), vp(vp_in), i_in(i_in_in), i_spin(i_spin_in), channel(ch_in), diff(diff_in), diag_class(k_in){
        set_Keldysh_index_i0(i0_in);
        if (MAX_DIAG_CLASS < 2) {precompute_vertices();}
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

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
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

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::precompute_vertices(){
    // For K1 class, left and right vertices do not depend on integration frequency
    // -> precompute them to save time
#ifndef DEBUG_SYMMETRIES
#ifdef KELDYSH_FORMALISM
    std::vector<int> indices = indices_sum(i0, i2, channel);

    VertexInput input_l (indices[0], spin, w, 0., 0., i_in, channel);
    VertexInput input_r (indices[1], spin, w, 0., 0., i_in, channel);
#else
    VertexInput input_l (0, spin, w, 0., 0., i_in, channel);
    VertexInput &input_r = input_l;
#endif
    res_l_V_initial = vertex1.left_same_bare(input_l);
    res_r_V_initial = vertex2.right_same_bare(input_r);
    if (channel == 't') {
        input_l.spin = 1;
        input_r.spin = 1;
        res_l_Vhat_initial = vertex1.left_same_bare(input_l);
        res_r_Vhat_initial = vertex2.right_same_bare(input_r);
    }
#else // DEBUG_SYMMETRIES

#ifdef KELDYSH_FORMALISM
    std::vector<int> indices = indices_sum(i0, i2, channel);

    VertexInput input_l (indices, spin, w, 0., 0., i_in, channel);
    VertexInput input_r (indices[1], spin, w, 0., 0., i_in, channel);
#else
    VertexInput input_l (0, spin, w, 0., 0., i_in, channel);
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

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
auto Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::operator()(double vpp) const -> Q {
    Q result;
#ifndef SWITCH_SUM_N_INTEGRAL

    if (case_always_has_to_be_zero()) {return 0.;}
    Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
    compute_vertices(vpp, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat);

    Q Pival = Pi.value(i2, w, vpp, i_in, channel);

    if (channel != 't')  // no spin sum in a and p channel
        result = res_l_V * Pival * res_r_V;
    else {
        // in t channel, spin sum has 3 terms:
        // result = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
        //        = 2. * res_l_V * Pival * res_r_V + res_l_V * Pival * res_r_Vhat + res_l_Vhat * Pival * res_r_V;
        switch (i_spin) {
            case 0:
                result = 2. * res_l_V * Pival * res_r_V;
                break;
            case 1:
                result = res_l_V * Pival * res_r_Vhat;
                break;
            case 2:
                result = res_l_Vhat * Pival * res_r_V;
                break;
            default:;
        }
    }
#ifdef DEBUG_SYMMETRIES
    if (spin == 0 and channel ==  'p') {
        result = (res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat) * 0.5;
    }
    if (spin == 1) {
        if (channel == 't')  // no spin sum in t channel
            result = res_l_V * Pival * res_r_V;
        else if (channel == 'p') {
            result = (res_l_V * Pival * res_r_Vhat + res_l_Vhat * Pival * res_r_V) * 0.5;
        }
        else { // channel == 'a'
            // in a channel, spin sum has 3 terms:
            // result = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
            //        = 2. * res_l_V * Pival * res_r_V + res_l_V * Pival * res_r_Vhat + res_l_Vhat * Pival * res_r_V;
            switch (i_spin) {
                case 0:
                    result = 2. * res_l_V * Pival * res_r_V;
                    break;
                case 1:
                    result = res_l_V * Pival * res_r_Vhat;
                    break;
                case 2:
                    result = res_l_Vhat * Pival * res_r_V;
                    break;
                default:;
            }
        }
    }
#endif // DEBUG_SYMMETRIES


#else

    /// TODO: fix conflict with sum over internal spins
    //assert(false);
    VertexInput input_external (i0, spin, w, v, vp, i_in, channel, diag_class, iw);

    if (channel == 'a')
        result = sum_over_internal<'a',Q>([&](const VertexInput& input_l) -> Q {if (diag_class == k1 or diag_class == k2b) return vertex1.left_same_bare(input_l) ; else return vertex1.left_diff_bare(input_l);},
                      Pi,
                      [&](const VertexInput& input_r) -> Q {if (diag_class == k3 or diag_class == k2b) return vertex2.right_diff_bare(input_r); else return vertex2.right_same_bare(input_r);},
                      input_external, vpp
    );
    else if (channel == 'p') {
        result = sum_over_internal<'p',Q>([&](const VertexInput& input_l) -> Q {if (diag_class == k1 or diag_class == k2b) return vertex1.left_same_bare(input_l) ; else return vertex1.left_diff_bare(input_l);},
                               Pi,
                               [&](const VertexInput& input_r) -> Q {if (diag_class == k3 or diag_class == k2b) return vertex2.right_diff_bare(input_r); else return vertex2.right_same_bare(input_r);},
                               input_external, vpp
        );
    }
    else {
        result = sum_over_internal<'t',Q>([&](const VertexInput& input_l) -> Q {if (diag_class == k1 or diag_class == k2b) return vertex1.left_same_bare(input_l) ; else return vertex1.left_diff_bare(input_l);},
                               Pi,
                               [&](const VertexInput& input_r) -> Q {if (diag_class == k3 or diag_class == k2b) return vertex2.right_diff_bare(input_r); else return vertex2.right_same_bare(input_r);},
                               input_external, vpp
        );
    }
#endif

    return result;
}

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
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

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
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
        VertexInput input_l (indices[0], spin, w, v, vpp,  i_in, channel, diag_class, iw);
        VertexInput input_r (indices[1], spin, w, vpp, vp, i_in, channel, diag_class, iw);

        if (i_spin == 0) { // first summand in all channels is res_l_V * Pival * res_r_V
            if (diag_class == k1 or diag_class == k2b)
                res_l_V = vertex1.left_same_bare(input_l);
            else
                res_l_V = vertex1.left_diff_bare(input_l);

            if (diag_class == k3 or diag_class == k2b)
                res_r_V = vertex2.right_diff_bare(input_r);
            else
                res_r_V = vertex2.right_same_bare(input_r);
#ifdef DEBUG_SYMMETRIES
            if (channel == 'p') {  // res_l_Vhat * Pival * res_r_Vhat
                // compute res_l_Vhat
                input_l.spin = 1 - spin;
                if (diag_class == k1 or diag_class == k2b)
                    res_l_Vhat = vertex1.left_same_bare(input_l);
                else
                    res_l_Vhat = vertex1.left_diff_bare(input_l);
                // compute res_r_Vhat
                input_r.spin = 1 - spin;
                if (diag_class == k3 or diag_class == k2b)
                    res_r_Vhat = vertex2.right_diff_bare(input_r);
                else
                    res_r_Vhat = vertex2.right_same_bare(input_r);
            }
#endif
        }
        else { // relevant for t-channel (spin component 0) and a-channel (spin component 1): there i_spin = 0, 1, 2

            if (channel == 't'
#ifdef DEBUG_SYMMETRIES
                or channel == 'a'
#endif
            ) {
                // channel = t, i_spin = 1
                if (i_spin == 1) {  // res_l_V * Pival * res_r_Vhat
                    // compute res_l_V
                    if (diag_class == k1 or diag_class == k2b)
                        res_l_V = vertex1.left_same_bare(input_l);
                    else
                        res_l_V = vertex1.left_diff_bare(input_l);
                    // compute res_r_Vhat
                    input_r.spin = 1 - spin;
                    if (diag_class == k3 or diag_class == k2b)
                        res_r_Vhat = vertex2.right_diff_bare(input_r);
                    else
                        res_r_Vhat = vertex2.right_same_bare(input_r);
                }
                // channel = t, i_spin = 2
                else {              // res_l_Vhat * Pival * res_r_V
                    assert(i_spin == 2);
                    // compute res_r_V
                    if (diag_class == k3 or diag_class == k2b)
                        res_r_V = vertex2.right_diff_bare(input_r);
                    else
                        res_r_V = vertex2.right_same_bare(input_r);
                    // compute res_l_Vhat
                    input_l.spin = 1 - spin;
                    if (diag_class == k1 or diag_class == k2b)
                        res_l_Vhat = vertex1.left_same_bare(input_l);
                    else
                        res_l_Vhat = vertex1.left_diff_bare(input_l);
                }
            }
        }
    }
}

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
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
            integrand_re[i] = myreal(integrand_value);
            integrand_im[i] = 0.;
            Pival_re[i] = myreal(Pival);
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

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::save_integrand() const {
    /// Define standard frequency points on which to evaluate the integrand
    int npoints = 10000;//nBOS;
    if (diag_class == k2) {npoints = 1000;}
    else if (diag_class == k3) {npoints = 100;}

    rvec freqs (npoints);

    for (int i=0; i<npoints; ++i) {
        double wl, wu;
        switch (diag_class) {
            case k1:
                wl = vertex1.avertex().K1.K1_get_wlower();
                wu = vertex1.avertex().K1.K1_get_wupper();
                wl *= 2.;
                wu *= 2.;
                break;
            case k2:
                wl = vertex1.avertex().K2.K2_get_wlower_f();
                wu = vertex1.avertex().K2.K2_get_wupper_f();
                break;
            case k3:
                wl = vertex1.avertex().K3.K3_get_wlower_f();
                wu = vertex1.avertex().K3.K3_get_wupper_f();
                break;
            default:;
        }
        double vpp = wl + i * (wu - wl) / (npoints - 1);
        if (diag_class == k1 and not HUBBARD_MODEL) {vertex1.avertex().K1.K1_get_freq_w(vpp, i);}
        freqs[i] = vpp;
    }

    std::string filename_prefix = "";
    if (HUBBARD_MODEL) filename_prefix = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/integrands/";
    save_integrand(freqs, filename_prefix);

}

template<typename Q, symmetryType symmetry_left, symmetryType symmetry_right, class Bubble_Object>
void Integrand<Q, symmetry_left, symmetry_right, Bubble_Object>::save_integrand(const rvec& freqs, const std::string& filename_prefix) const {
    /// evaluate the integrand on frequency points in freqs
    int npoints = freqs.size();

    rvec integrand_re (npoints);
    rvec integrand_im (npoints);
    rvec Pival_re (npoints);
    rvec Pival_im (npoints);

    get_integrand_vals(freqs, integrand_re, integrand_im, Pival_re, Pival_im);

    std::string filename = "";
    if (not HUBBARD_MODEL) filename += data_dir;
    filename += filename_prefix+"integrand_K" + std::to_string(diag_class);
    filename += channel;
    filename += "_i0=" + std::to_string(i0)
                + "_i2=" + std::to_string(i2)
                + "_w=" + std::to_string(w);
    if (diag_class == k2) {filename += "_v=" + std::to_string(v);}
    else if (diag_class == k3) {filename += "_vp=" + std::to_string(vp);}
    filename += "_i_in=" + std::to_string(i_in);
    filename += + ".h5";
    write_h5_rvecs(filename,
                   {"v", "integrand_re", "integrand_im", "Pival_re", "Pival_im"},
                   {freqs, integrand_re, integrand_im, Pival_re, Pival_im});
}


#endif //KELDYSH_MFRG_INTEGRAND_HPP
