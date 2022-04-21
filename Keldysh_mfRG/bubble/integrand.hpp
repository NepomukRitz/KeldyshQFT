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
#include "bubble.hpp"
#include "precalculated_bubble.hpp"
#include "../multidimensional/multiarray.hpp"

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
    const char channel;

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
template <K_class diag_class, char channel, int spin, typename Q,
        vertexType symmetry_left,
        vertexType symmetry_right,
        class Bubble_Object,
        typename return_type = Q>
class Integrand {
private:
    const GeneralVertex<Q, symmetry_left>& vertex1;
    const GeneralVertex<Q, symmetry_right>& vertex2;
    const Bubble_Object& Pi;
    int i0_symmred;
    int i0 = 0;
    int i0_left;
    int i0_right;
    const int i2;
    const int i_in;
    const int i_spin;
    const int iw=0;
    const double w, v = 0., vp = 0.;
    const bool diff;

#if KELDYSH_FORMALISM
    using buffer_type_vertex_l = Eigen::Matrix<Q,myRowsAtCompileTime<return_type>(),4>;
    using buffer_type_vertex_r = Eigen::Matrix<Q,4,myColsAtCompileTime<return_type>()>; // buffer_type_vertex_l;
#else
    using buffer_type_vertex_l = Q;
    using buffer_type_vertex_r = Q;
#endif

    //K_class diag_class;

    Q res_l_V_initial, res_r_V_initial, res_l_Vhat_initial, res_r_Vhat_initial; // To be precomputed for K1

    void set_Keldysh_index_i0(int i0_in);
    void set_Keldysh_index_i0_left_right(int i0_in);
    void precompute_vertices();

    bool case_always_has_to_be_zero() const;

    void compute_vertices(double vpp, Q& res_l_V, Q& res_r_V, Q& res_l_Vhat, Q& res_r_Vhat) const;

    template<int ispin> void load_vertex_keldyshComponents_left_scalar (buffer_type_vertex_l& values_vertex, const VertexInput& input) const;
    template<int ispin> void load_vertex_keldyshComponents_right_scalar(buffer_type_vertex_r& values_vertex, const VertexInput& input) const;
    template<int ispin> void load_vertex_keldyshComponents_left_vectorized (buffer_type_vertex_l& values_vertex, const VertexInput& input) const;
    template<int ispin> void load_vertex_keldyshComponents_right_vectorized(buffer_type_vertex_r& values_vertex, const VertexInput& input) const;

    Q sum_over_internal_scalar(const VertexInput& input_external, double vpp) const;
    return_type sum_over_internal_vectorized(const VertexInput& input_external, double vpp) const;
        public:
    /**
     * Constructor for asymptotic class Ki:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index specifying the (external) Keldysh component of integrand object; i0_in = [0, .. nK_Ki]
     *                     where nK_Ki is the number of symmetry-reduced Keldysh components in the result
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     * @param diff_in    : determines whether to compute differentiated or non-differentiated bubble
     */
    Integrand(const GeneralVertex<Q, symmetry_left>& vertex1_in,
              const GeneralVertex<Q, symmetry_right>& vertex2_in,
              const Bubble_Object& Pi_in,
              int i0_in, int i2_in, const int iw_in, const double w_in, const double v_in, const double vp_in, const int i_in_in,
              const int i_spin_in, const bool diff_in)
              :vertex1(vertex1_in), vertex2(vertex2_in), Pi(Pi_in), i0_symmred(i0_in),
              i2(i2_in), iw(iw_in), w(w_in), v(v_in), vp(vp_in), i_in(i_in_in), i_spin(i_spin_in), diff(diff_in){
        set_Keldysh_index_i0(i0_in);
        if (MAX_DIAG_CLASS < 2) {precompute_vertices();}
    }

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> return_type;

    void save_integrand() const;
    void save_integrand(const rvec& freqs, const std::string& filename_prefix) const;
    void get_integrand_vals(const rvec& freqs, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& integrand_vals, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& Pivals, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& vertex_vals1, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& vertex_vals2)  const;
};

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::set_Keldysh_index_i0(const int i0_in) {
    if constexpr(KELDYSH){
        if constexpr(not DEBUG_SYMMETRIES and not VECTORIZED_INTEGRATION) {
            switch (diag_class) {
                case k1: // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
                    switch (channel) {
                        case 'a':
                            i0 = non_zero_Keldysh_K1a[i0_in];
                            break;
                        case 'p':
                            i0 = non_zero_Keldysh_K1p[i0_in];
                            break;
                        case 't':
                            i0 = non_zero_Keldysh_K1t[i0_in];
                            break;
                        default:;
                    }
                    break;
                case k2: // converting index i0_in (0,...,4) into actual Keldysh index i0 (0,...,15)
                    switch (channel) {
                        case 'a':
                            i0 = non_zero_Keldysh_K2a[i0_in];
                            break;
                        case 'p':
                            i0 = non_zero_Keldysh_K2p[i0_in];
                            break;
                        case 't':
                            i0 = non_zero_Keldysh_K2t[i0_in];
                            break;
                        default:;
                    }
                    break;
                case k3:
                    i0 = non_zero_Keldysh_K3[i0_in]; // converting index i0_in (0,...,5) into actual Keldysh index i0 (0,...,15)
                    break;
                default:;
            }
        }
        else {
            i0 = i0_in;
        }
        set_Keldysh_index_i0_left_right(i0);
    }
    else{
        i0 = 0;
        i0_left = 0;
        i0_right= 0;
    }
}

// i0_in in {0,...,15}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::set_Keldysh_index_i0_left_right(const int i0_in) {
    my_index_t left, right;
    get_i0_left_right<channel>(i0_in, left, right);
    i0_left  = left * 4; //
    i0_right = right * 4; //
}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::precompute_vertices() {
    // For K1 class, left and right vertices do not depend on integration frequency
    // -> precompute them to save time
    if constexpr(not DEBUG_SYMMETRIES) {
#if KELDYSH_FORMALISM
        std::vector<int> indices = indices_sum(i0, i2, channel);

        VertexInput input_l(indices[0], spin, w, 0., 0., i_in, channel);
        VertexInput input_r(indices[1], spin, w, 0., 0., i_in, channel);
#else
        VertexInput input_l (0, spin, w, 0., 0., i_in, channel);
        VertexInput &input_r = input_l;
#endif
        res_l_V_initial = vertex1.template left_same_bare<channel>(input_l);
        res_r_V_initial = vertex2.template right_same_bare<channel>(input_r);
        if (channel == 't') {
            input_l.spin = 1;
            input_r.spin = 1;
            res_l_Vhat_initial = vertex1.template left_same_bare<channel>(input_l);
            res_r_Vhat_initial = vertex2.template right_same_bare<channel>(input_r);
        }
    }
    else {// DEBUG_SYMMETRIES

    #if KELDYSH_FORMALISM
        std::vector<int> indices = indices_sum(i0, i2, channel);
        VertexInput input_l(indices[0], spin, w, 0., 0., i_in, channel);
        VertexInput input_r(indices[1], spin, w, 0., 0., i_in, channel);
    #else
        VertexInput input_l (0, spin, w, 0., 0., i_in, channel);
        VertexInput &input_r = input_l;
    #endif
        res_l_V_initial = vertex1.template left_same_bare<channel>(input_l);
        res_r_V_initial = vertex2.template right_same_bare<channel>(input_r);
        if (channel == 't' and spin == 0) {
            input_l.spin = 1 - input_l.spin; // flip spin 0 <-> 1
            input_r.spin = 1 - input_r.spin; // flip spin 0 <-> 1
            res_l_Vhat_initial = vertex1.template left_same_bare<channel>(input_l);
            res_r_Vhat_initial = vertex2.template right_same_bare<channel>(input_r);
        } else if (channel == 'a' and spin == 1) {
            input_l.spin = 1 - input_l.spin; // flip spin 0 <-> 1
            input_r.spin = 1 - input_r.spin; // flip spin 0 <-> 1
            res_l_Vhat_initial = vertex1.template left_same_bare<channel>(input_l);
            res_r_Vhat_initial = vertex2.template right_same_bare<channel>(input_r);
        } else if (channel == 'p' and spin == 1) {
            input_l.spin = 1 - input_l.spin; // flip spin 0 <-> 1
            input_r.spin = 1 - input_r.spin; // flip spin 0 <-> 1
            res_l_Vhat_initial = vertex1.template left_same_bare<channel>(input_l);
            res_r_Vhat_initial = vertex2.template right_same_bare<channel>(input_r);
        }
    } // DEBUG_SYMMETRIES
}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
auto Integrand<diag_class, channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::operator()(double vpp) const -> return_type {
    return_type result;
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
#if DEBUG_SYMMETRIES
    if (spin == 0 and channel ==  'p') {
        result = (res_l_V * Pival * res_r_V);
    }
    if (spin == 1) {
        if (channel == 't')  // no spin sum in t channel
            result = res_l_V * Pival * res_r_V;
        else if (channel == 'p') {
            result = (res_l_V * Pival * res_r_Vhat);
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
    VertexInput input_external (i0, 0, w, v, vp, i_in, channel, diag_class, iw);

    //result = sum_over_internal_scalar(input_external, vpp);
    if constexpr(KELDYSH)
    {
        result = sum_over_internal_vectorized(input_external, vpp);
    }
    else {
        result = sum_over_internal_scalar(input_external, vpp);
    }
    /// comment in to test vectorized access
    //Q result2 = sum_over_internal_scalar(input_external, vpp);
    //assert(std::abs(result-result2) < 1e-10);

#endif

    return result;
}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
bool Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::case_always_has_to_be_zero() const {
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

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::compute_vertices(const double vpp,
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
                res_l_V = vertex1.template left_same_bare<channel>(input_l);
            else
                res_l_V = vertex1.template left_diff_bare<channel>(input_l);

            if (diag_class == k3 or diag_class == k2b)
                res_r_V = vertex2.template right_diff_bare<channel>(input_r);
            else
                res_r_V = vertex2.template right_same_bare<channel>(input_r);
#if DEBUG_SYMMETRIES
            if (channel == 'p') {  // res_l_Vhat * Pival * res_r_Vhat
                // compute res_l_Vhat
                input_l.spin = 1 - spin;
                if (diag_class == k1 or diag_class == k2b)
                    res_l_Vhat = vertex1.template left_same_bare<channel>(input_l);
                else
                    res_l_Vhat = vertex1.template left_diff_bare<channel>(input_l);
                // compute res_r_Vhat
                input_r.spin = 1 - spin;
                if (diag_class == k3 or diag_class == k2b)
                    res_r_Vhat = vertex2.template right_diff_bare<channel>(input_r);
                else
                    res_r_Vhat = vertex2.template right_same_bare<channel>(input_r);
            }
#endif
        }
        else { // relevant for t-channel (spin component 0) and a-channel (spin component 1): there i_spin = 0, 1, 2

            if constexpr(channel == 't'
#if DEBUG_SYMMETRIES
                or channel == 'a'
#endif
            ) {
                // channel = t, i_spin = 1
                if (i_spin == 1) {  // res_l_V * Pival * res_r_Vhat
                    // compute res_l_V
                    if (diag_class == k1 or diag_class == k2b)
                        res_l_V = vertex1.template left_same_bare<channel>(input_l);
                    else
                        res_l_V = vertex1.template left_diff_bare<channel>(input_l);
                    // compute res_r_Vhat
                    input_r.spin = 1 - spin;
                    if (diag_class == k3 or diag_class == k2b)
                        res_r_Vhat = vertex2.template right_diff_bare<channel>(input_r);
                    else
                        res_r_Vhat = vertex2.template right_same_bare<channel>(input_r);
                }
                // channel = t, i_spin = 2
                else {              // res_l_Vhat * Pival * res_r_V
                    assert(i_spin == 2);
                    // compute res_r_V
                    if (diag_class == k3 or diag_class == k2b)
                        res_r_V = vertex2.template right_diff_bare<channel>(input_r);
                    else
                        res_r_V = vertex2.template right_same_bare<channel>(input_r);
                    // compute res_l_Vhat
                    input_l.spin = 1 - spin;
                    if (diag_class == k1 or diag_class == k2b)
                        res_l_Vhat = vertex1.template left_same_bare<channel>(input_l);
                    else
                        res_l_Vhat = vertex1.template left_diff_bare<channel>(input_l);
                }
            }
        }
    }
}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
template<int ispin>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::load_vertex_keldyshComponents_left_scalar(buffer_type_vertex_l& values_vertex, const VertexInput& input) const {
    //size_t len_1 = values_vertex.length()[0];
    const auto vertexvalue_scalar = [&] (const VertexInput& input_l) {if constexpr(diag_class == k1 or diag_class == k2b) return vertex1.template left_same_bare_symmetry_expanded<ispin,channel,Q>(input_l) ; else return vertex1.template left_diff_bare_symmetry_expanded<ispin,channel,Q>(input_l);};
    //const auto vertexvalue_scalar = [&] (const VertexInput& input_l) {if constexpr(diag_class == k1 or diag_class == k2b) return vertex1.template left_same_bare<channel>(input_l) ; else return vertex1.template left_diff_bare<channel>(input_l);};
    //assert(len_1 == glb_number_of_Keldysh_components_bubble);

    if constexpr(not KELDYSH) {
        values_vertex =  vertexvalue_scalar(input);
    }
    else {



        //Q v11, v12, v21, v22;
        //const std::array<size_t,2> dims_vtemp = {2,2};
        // v_temp contains values for different Keldysh indices:
        // i.e. in a-channel v_temp(i,j) = v_{1',i |j,2}
        //      in p-channel v_temp(i,j) = v_{1',2'|j,i}
        //      in t-channel v_temp(i,j) = v_{i ,2'|j,2}
        //multidimensional::multiarray<Q,2> v_temp(dims_vtemp);
        std::vector<int> alpha = alphas(input.iK);

        // fill v_temp:
        const std::array<size_t,4> dims_K = {2,2, 2,2};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                //if constexpr(channel == 'a') {input_tmp.iK = getFlatIndex<4,int,int,int,int>(alpha[0]-1, i, j, alpha[3]-1, dims_K);}
                //else if constexpr(channel == 'p') {input_tmp.iK = getFlatIndex<4,int,int,int,int>(alpha[0]-1, alpha[1]-1, j, i, dims_K);}
                //else if constexpr(channel == 't') {input_tmp.iK = getFlatIndex<4,int,int,int,int>(i, alpha[1]-1, j, alpha[3]-1, dims_K);}
                //else assert(false);
                VertexInput input_tmp = input;
                input_tmp.iK = input.iK + i*2+j;

                values_vertex[i*2+j] = vertexvalue_scalar(input_tmp);
            }
        }
        /*
        // fill values_vertex:
        if constexpr(channel == 'a' or channel == 't') {
            values_vertex(1, spin_idx) =  v_temp(0,0);
            values_vertex(0, spin_idx) = values_vertex(2, spin_idx) = v_temp(1,0);
            values_vertex(5, spin_idx) = values_vertex(7, spin_idx) = v_temp(0,1);
            values_vertex(3, spin_idx) = values_vertex(4, spin_idx) = values_vertex(6, spin_idx) = values_vertex(8, spin_idx) = v_temp(1,1);
        }
        else if constexpr(channel == 'p') {
            values_vertex(0, spin_idx) =  v_temp(0,0);
            values_vertex(1, spin_idx) = values_vertex(2, spin_idx) = v_temp(1,0);
            values_vertex(3, spin_idx) = values_vertex(4, spin_idx) = v_temp(0,1);
            values_vertex(5, spin_idx) = values_vertex(6, spin_idx) = values_vertex(7, spin_idx) = values_vertex(8, spin_idx) = v_temp(1,1);
        }
        else assert(false);
        */


    }

}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
template<int ispin>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::load_vertex_keldyshComponents_right_scalar(buffer_type_vertex_r& values_vertex, const VertexInput& input) const {
    //size_t len_1 = values_vertex.length()[0];
    const auto vertexvalue_scalar = [&](const VertexInput& input_r) {if constexpr(diag_class == k3 or diag_class == k2b) return vertex2.template right_diff_bare_symmetry_expanded<ispin,channel,Q>(input_r); else return vertex2.template right_same_bare_symmetry_expanded<ispin,channel,Q>(input_r);};
    //const auto vertexvalue_scalar = [&](const VertexInput& input_r) {if constexpr(diag_class == k3 or diag_class == k2b) return vertex2.template right_diff_bare<channel>(input_r); else return vertex2.template right_same_bare<channel>(input_r);};
    //assert(len_1 == glb_number_of_Keldysh_components_bubble);

    if constexpr(not KELDYSH) {
        values_vertex =  vertexvalue_scalar(input);
    }
    else {
        //Q v11, v12, v21, v22;
        //const std::array<size_t, 2> dims_vtemp = {2, 2};
        // v_temp contains values for different Keldysh indices:
        // i.e. in a-channel v_temp(i,j) = v_{i ,2'|1,j}
        //      in p-channel v_temp(i,j) = v_{j ,i |1,2}
        //      in t-channel v_temp(i,j) = v_{1',i |1,j}
        //multidimensional::multiarray<Q, 2> v_temp(dims_vtemp);
        std::vector<int> alpha = alphas(input.iK);

        // fill v_temp:
        const std::array<size_t, 4> dims_K = {2, 2, 2, 2};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                //if constexpr(channel == 'a') { input_tmp.iK = getFlatIndex<4,int,int,int,int>(i, alpha[1] - 1, alpha[2] - 1, j, dims_K); }
                //else if constexpr(channel == 'p') { input_tmp.iK = getFlatIndex<4,int,int,int,int>(j, i, alpha[2] - 1, alpha[3] - 1, dims_K); }
                //else if constexpr(channel == 't') { input_tmp.iK = getFlatIndex<4,int,int,int,int>(alpha[0] - 1, i, alpha[2] - 1, j, dims_K); }
                //else
                //    assert(false);

                VertexInput input_tmp = input;
                input_tmp.iK = input.iK + i*2+j;

                values_vertex[i*2+j] = vertexvalue_scalar(input_tmp);
            }
        }
        /*
        // fill values_vertex:
        if constexpr(channel == 'a' or channel == 't') {
            values_vertex[3] = v_temp(0, 0);
            values_vertex[5] = values_vertex[6] = v_temp(0, 1);
            values_vertex[0] = values_vertex[4] = v_temp(1, 0);
            values_vertex[1] = values_vertex[2] = values_vertex[7] = values_vertex[8] = v_temp(
                    1, 1);
        } else if constexpr(channel == 'p') {
            values_vertex[5] = v_temp(0, 0);
            values_vertex[3] = values_vertex[6] = v_temp(1, 0);
            values_vertex[1] = values_vertex[7] = v_temp(0, 1);
            values_vertex[0] = values_vertex[2] = values_vertex[4] = values_vertex[8] = v_temp(1, 1);
        } else
            assert(false);
            */
    }
}


template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
template<int ispin>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::load_vertex_keldyshComponents_left_vectorized(buffer_type_vertex_l& values_vertex, const VertexInput& input) const {
    //size_t len_1 = values_vertex.length()[0];
    using result_type_fetch = Eigen::Matrix<Q, 4 * myRowsAtCompileTime<return_type>(),1>;
    const auto vertexvalue_vector = [&] (const VertexInput& input_l) {if constexpr(diag_class == k1 or diag_class == k2b) return vertex1.template left_same_bare_symmetry_expanded<ispin,channel,result_type_fetch>(input_l) ; else return vertex1.template left_diff_bare_symmetry_expanded<ispin,channel,result_type_fetch>(input_l);};
    //assert(len_1 == glb_number_of_Keldysh_components_bubble);

    if constexpr(not KELDYSH) {
        //static_assert(KELDYSH, "Currently vectorized integrand only available for Keldysh.");
        values_vertex = vertexvalue_vector(input);
    }
    else {
        Eigen::Matrix<Q,Eigen::Dynamic, Eigen::Dynamic> result_temp = vertexvalue_vector(input);
        result_temp.resize(4, myRowsAtCompileTime<return_type>());
        values_vertex = result_temp.transpose();
    }

}


template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
template<int ispin>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::load_vertex_keldyshComponents_right_vectorized(buffer_type_vertex_r& values_vertex, const VertexInput& input) const {
    //size_t len_1 = values_vertex.length()[0];
    using result_type_fetch = Eigen::Matrix<Q, 4 * myColsAtCompileTime<return_type>(), 1>;
    const auto vertexvalue_vector = [&](const VertexInput& input_r) {if constexpr(diag_class == k3 or diag_class == k2b) return vertex2.template right_diff_bare_symmetry_expanded<ispin,channel,result_type_fetch>(input_r); else return vertex2.template right_same_bare_symmetry_expanded<ispin,channel,result_type_fetch>(input_r);};
    //assert(len_1 == glb_number_of_Keldysh_components_bubble);

    if constexpr(not KELDYSH) {
        //static_assert(KELDYSH, "Currently vectorized integrand only available for Keldysh.");
        values_vertex = vertexvalue_vector(input);
    }
    else {

        Eigen::Matrix<Q,Eigen::Dynamic, Eigen::Dynamic> result_temp = vertexvalue_vector(input);
        result_temp.resize(4, myColsAtCompileTime<return_type>());
        values_vertex = result_temp;
    }

}


template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
Q Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::sum_over_internal_scalar(const VertexInput& input_external, const double vpp) const {

    VertexInput input_l = input_external, input_r = input_external;
    input_l.v2 = vpp; input_r.v1 = vpp;
    input_l.iK = i0_left; input_r.iK = i0_right;
    input_l.spin = 0; input_r.spin = 0;

    // create multiarrays to store loaded values
    buffer_type_vertex_l values_vertex_l;
    buffer_type_vertex_r values_vertex_r;
    auto Pi_matrix = Pi.template value_vectorized<channel>(input_external.w, vpp, input_external.i_in);
    // load vertex values:
    int spin_idx = 0;
    load_vertex_keldyshComponents_left_scalar<spin> (values_vertex_l, input_l);
    load_vertex_keldyshComponents_right_scalar<spin>(values_vertex_r, input_r);

    /// Comment in to compare with vectorized access to vertices:
    //buffer_type_vertex_l values_vertex_l_test;
    //buffer_type_vertex_r values_vertex_r_test;
    //load_vertex_keldyshComponents_left_vectorized (values_vertex_l_test, input_l);
    //load_vertex_keldyshComponents_right_vectorized(values_vertex_r_test, input_r);
    //if (std::abs((values_vertex_r - values_vertex_r_test).sum()) > 1e-10) {
    //    H5::H5File file(data_dir + "values_vertex_r", H5F_ACC_TRUNC);
    //    write_to_hdf(file, "scalar", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_r), false);
    //    write_to_hdf(file, "vector", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_r_test), false);
    //    file.close();
    //    assert(false);
    //}
    //if (std::abs((values_vertex_l - values_vertex_l_test).sum()) > 1e-10) {
    //    H5::H5File file(data_dir + "values_vertex_l", H5F_ACC_TRUNC);
    //    write_to_hdf(file, "scalar", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_l), false);
    //    write_to_hdf(file, "vector", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_l_test), false);
    //    file.close();
    //    assert(false);
    //}


    Q result;
    // load other spin component
    if constexpr ((channel == 't' and spin == 0)
#if DEBUG_SYMMETRIES
        or (channel == 'a' and spin == 1)
#endif
            ){

        buffer_type_vertex_l values_vertex_l_other;
        buffer_type_vertex_r values_vertex_r_other;

        load_vertex_keldyshComponents_left_scalar <1-spin>(values_vertex_l_other, input_l);
        load_vertex_keldyshComponents_right_scalar<1-spin>(values_vertex_r_other, input_r);
        if constexpr(KELDYSH)
        {
            result = ((values_vertex_l + values_vertex_l_other) * Pi_matrix * values_vertex_r +
                      values_vertex_l * Pi_matrix * (values_vertex_r + values_vertex_r_other)).eval()[0];
        }
        else {
            result = ((values_vertex_l + values_vertex_l_other) * Pi_matrix * values_vertex_r +
                      values_vertex_l * Pi_matrix * (values_vertex_r + values_vertex_r_other));
        }
    }
#if DEBUG_SYMMETRIES
        /*else if (channel == 'p' and input_external.spin == 0) {
        buffer_type_vertex_l values_vertex_l_other;
        buffer_type_vertex_r values_vertex_r_other;
        input_l.spin = 1 - input_external.spin;
        input_r.spin = 1 - input_external.spin;
        spin_idx = 1;

        load_vertex_keldyshComponents_left_scalar (values_vertex_l_other, input_l);
        load_vertex_keldyshComponents_right_scalar(values_vertex_r_other, input_r);

            result = (values_vertex_l * Pi_matrix * values_vertex_r +  values_vertex_l_other * Pi_matrix * values_vertex_r_other) * 0.5;
        }*/
    else  if constexpr(channel == 'p' and spin == 1) {
        buffer_type_vertex_l values_vertex_l_other;
        buffer_type_vertex_r values_vertex_r_other;

        load_vertex_keldyshComponents_left_scalar <1-spin>(values_vertex_l_other, input_l);
        load_vertex_keldyshComponents_right_scalar<1-spin>(values_vertex_r_other, input_r);

            result = (values_vertex_l * Pi_matrix * values_vertex_r_other);
        }
#endif
    else {
        result = values_vertex_l * Pi_matrix * values_vertex_r;
    }



    /*
    if constexpr(false) { // set to true for debugging the loading of vertex values

        multidimensional::multiarray<Q,2> values_vertex_l_alt(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2})), values_vertex_r_alt(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2}));

        // load values (basic inefficient version):
        for (int i = 0; i < glb_number_of_Keldysh_components_bubble; i++) {

            VertexInput input_l_alt = input_external, input_r_alt = input_external;
            input_l_alt.v2 = vpp; input_r_alt.v1 = vpp;

            int i2 = glb_non_zero_Keldysh_bubble[i];
            std::vector<int> indices_vertices = indices_sum(input_external.iK, i2, channel);
            input_l_alt.iK = indices_vertices[0];
            input_r_alt.iK = indices_vertices[1];

            values_vertex_l_alt.at(i,0) = value_vertex_l(input_l_alt);
            values_vertex_r_alt.at(i,0) = value_vertex_r(input_r_alt);

            // load other spin component
            if ((channel == 't' and input_external.spin == 0)
#if DEBUG_SYMMETRIES
                or (channel == 'a' and input_external.spin == 1) or (channel == 'p')
#endif
                    ){
                if(channel == 't') assert(0 == input_external.spin);
                if(channel == 'a') assert(1 == input_external.spin);
                input_l_alt.spin = 1 - input_external.spin;
                input_r_alt.spin = 1 - input_external.spin;

                values_vertex_l_alt.at(i,1) = value_vertex_l(input_l_alt);
                values_vertex_r_alt.at(i,1) = value_vertex_r(input_r_alt);

            }
        }



        multidimensional::multiarray<Q,2> values_vertex_l_diff = values_vertex_l - values_vertex_l_alt;
        multidimensional::multiarray<Q,2> values_vertex_r_diff = values_vertex_r - values_vertex_r_alt;
        auto diff_l = values_vertex_l_diff.maxabs();
        auto diff_r = values_vertex_r_diff.maxabs();
        assert(std::abs(diff_l) < 1e-10);
        assert(std::abs(diff_r) < 1e-10);


    }
    */
    return result;
}


template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
return_type Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::sum_over_internal_vectorized(const VertexInput& input_external, const double vpp) const {

    VertexInput input_l = input_external, input_r = input_external;
    input_l.v2 = vpp; input_r.v1 = vpp;
    input_l.iK = i0_left; input_r.iK = i0_right;
    input_l.spin = 0;
    input_r.spin = 0;

    // create multiarrays to store loaded values
    buffer_type_vertex_l values_vertex_l;
    buffer_type_vertex_r values_vertex_r;

    auto Pi_matrix = Pi.template value_vectorized<channel>(input_external.w, vpp, input_external.i_in);
    // load vertex values:
    int spin_idx = 0;
    load_vertex_keldyshComponents_left_vectorized <spin>(values_vertex_l, input_l);
    load_vertex_keldyshComponents_right_vectorized<spin>(values_vertex_r, input_r);

    /// Comment in to compare with vectorized access to vertices:
    //buffer_type_vertex_l values_vertex_l_test;
    //buffer_type_vertex_r values_vertex_r_test;
    //load_vertex_keldyshComponents_left_scalar (values_vertex_l_test, input_l);
    //load_vertex_keldyshComponents_right_scalar(values_vertex_r_test, input_r);
    //if (std::abs((values_vertex_r - values_vertex_r_test).sum()) > 1e-10) {
    //    H5::H5File file(data_dir + "values_vertex_r", H5F_ACC_TRUNC);
    //    write_to_hdf(file, "scalar", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_r_test), false);
    //    write_to_hdf(file, "vector", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_r), false);
    //    file.close();
    //    assert(false);
    //}
    //if (std::abs((values_vertex_l - values_vertex_l_test).sum()) > 1e-10) {
    //    H5::H5File file(data_dir + "values_vertex_l", H5F_ACC_TRUNC);
    //    write_to_hdf(file, "scalar", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_l), false);
    //    write_to_hdf(file, "vector", multidimensional::multiarray<Q,1>(std::array<size_t,1>({4}), values_vertex_l_test), false);
    //    file.close();
    //    assert(false);
    //}

    return_type result;
    // load other spin component
    if constexpr ((channel == 't' and spin == 0)
#if DEBUG_SYMMETRIES
        or (channel == 'a' and spin == 1)
#endif
            ){

        buffer_type_vertex_l values_vertex_l_other;
        buffer_type_vertex_r values_vertex_r_other;

        load_vertex_keldyshComponents_left_vectorized <1-spin>(values_vertex_l_other, input_l);
        load_vertex_keldyshComponents_right_vectorized<1-spin>(values_vertex_r_other, input_r);
        if constexpr(VECTORIZED_INTEGRATION) {
            assert(input_external.iK == 0);
            result = ((values_vertex_l + values_vertex_l_other) * Pi_matrix * values_vertex_r +
                      values_vertex_l * Pi_matrix * (values_vertex_r + values_vertex_r_other)).eval();
        }
        else {
            result = ((values_vertex_l + values_vertex_l_other) * Pi_matrix * values_vertex_r +
                      values_vertex_l * Pi_matrix * (values_vertex_r + values_vertex_r_other)).eval()[0];
        }
    }
#if DEBUG_SYMMETRIES
    //else if (channel == 'p' and input_external.spin == 0) {
    //    buffer_type_vertex_l values_vertex_l_other;
    //    buffer_type_vertex_r values_vertex_r_other;
    //    input_l.spin = 1 - input_external.spin;
    //    input_r.spin = 1 - input_external.spin;
    //    spin_idx = 1;
//
    //    load_vertex_keldyshComponents_left_scalar (values_vertex_l_other, input_l);
    //    load_vertex_keldyshComponents_right_scalar(values_vertex_r_other, input_r);
//
    //    result = (values_vertex_l * Pi_matrix * values_vertex_r +  values_vertex_l_other * Pi_matrix * values_vertex_r_other).eval()[0] * 0.5;
    //}
    else if constexpr(channel == 'p' and spin == 1) {
        buffer_type_vertex_l values_vertex_l_other;
        buffer_type_vertex_r values_vertex_r_other;

        load_vertex_keldyshComponents_left_vectorized <1-spin>(values_vertex_l_other, input_l);
        load_vertex_keldyshComponents_right_vectorized<1-spin>(values_vertex_r_other, input_r);
        if constexpr (VECTORIZED_INTEGRATION) {
            result = (values_vertex_l * Pi_matrix * values_vertex_r_other).eval();
        }
        else {
            result = (values_vertex_l * Pi_matrix * values_vertex_r_other).eval()[0];
        }
    }
#endif
    else {
        if constexpr(VECTORIZED_INTEGRATION) {
            result = (values_vertex_l * Pi_matrix * values_vertex_r).eval();
        }
        else {
            result = (values_vertex_l * Pi_matrix * values_vertex_r).eval()[0];
        }
    }



    /*
    if constexpr(false) { // set to true for debugging the loading of vertex values

        multidimensional::multiarray<Q,2> values_vertex_l_alt(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2})), values_vertex_r_alt(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2}));

        // load values (basic inefficient version):
        for (int i = 0; i < glb_number_of_Keldysh_components_bubble; i++) {

            VertexInput input_l_alt = input_external, input_r_alt = input_external;
            input_l_alt.v2 = vpp; input_r_alt.v1 = vpp;

            int i2 = glb_non_zero_Keldysh_bubble[i];
            std::vector<int> indices_vertices = indices_sum(input_external.iK, i2, channel);
            input_l_alt.iK = indices_vertices[0];
            input_r_alt.iK = indices_vertices[1];

            values_vertex_l_alt.at(i,0) = value_vertex_l(input_l_alt);
            values_vertex_r_alt.at(i,0) = value_vertex_r(input_r_alt);

            // load other spin component
            if ((channel == 't' and input_external.spin == 0)
#if DEBUG_SYMMETRIES
                or (channel == 'a' and input_external.spin == 1) or (channel == 'p')
#endif
                    ){
                if(channel == 't') assert(0 == input_external.spin);
                if(channel == 'a') assert(1 == input_external.spin);
                input_l_alt.spin = 1 - input_external.spin;
                input_r_alt.spin = 1 - input_external.spin;

                values_vertex_l_alt.at(i,1) = value_vertex_l(input_l_alt);
                values_vertex_r_alt.at(i,1) = value_vertex_r(input_r_alt);

            }
        }



        multidimensional::multiarray<Q,2> values_vertex_l_diff = values_vertex_l - values_vertex_l_alt;
        multidimensional::multiarray<Q,2> values_vertex_r_diff = values_vertex_r - values_vertex_r_alt;
        auto diff_l = values_vertex_l_diff.maxabs();
        auto diff_r = values_vertex_r_diff.maxabs();
        assert(std::abs(diff_l) < 1e-10);
        assert(std::abs(diff_r) < 1e-10);


    }
    */
    return result;
}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::get_integrand_vals(const rvec& freqs, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& integrand_vals, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& Pivals, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& vertex_vals1, Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>& vertex_vals2) const {
    int npoints = freqs.size();
    VertexInput input_external (i0, 0, w, v, vp, i_in, channel, diag_class, iw);

    Pivals = Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>(KELDYSH?16:1, npoints);
    integrand_vals = Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic>(VECTORIZED_INTEGRATION?16:1, npoints);
    for (int i=0; i<npoints; ++i) {

        double vpp = freqs[i];

        VertexInput input_l = input_external, input_r = input_external;
        input_l.v2 = vpp; input_r.v1 = vpp;
        input_l.iK = i0_left; input_r.iK = i0_right;
        input_l.spin = 0;
        input_r.spin = 0;

        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> Pival;
#if VECTORIZED_INTEGRATION
        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> integrand_value;
        Eigen::Matrix<Q,4,4> vertex_val1(4,4);
        Eigen::Matrix<Q,4,4> vertex_val2(4,4);
        load_vertex_keldyshComponents_left_vectorized<0> (vertex_val1, input_l);
        load_vertex_keldyshComponents_right_vectorized<0>(vertex_val2, input_r);
#else
        Q integrand_value;
#endif

        Pival = Pi.template value_vectorized<channel>(w, vpp, i_in);
        integrand_value = (*this)(vpp);


#if VECTORIZED_INTEGRATION
        integrand_value.resize(16,1);

        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> vertex_val1_temp = vertex_val1;
        Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> vertex_val2_temp = vertex_val2;
        vertex_val1_temp.resize(16,1);
        vertex_val2_temp.resize(16,1);
        integrand_vals.col(i) = integrand_value;
        vertex_vals1.col(i) = vertex_val1_temp;
        vertex_vals2.col(i) = vertex_val2_temp;
#else
        integrand_vals(i) = integrand_value;
#endif


        Pival.resize(KELDYSH?16:1,1);
        Pivals.col(i) = Pival;

    }


}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::save_integrand() const {
    /// Define standard frequency points on which to evaluate the integrand
    int npoints = 10000;//nBOS;
    if (diag_class == k2) {npoints = 1000;}
    else if (diag_class == k3) {npoints = 100;}

    rvec freqs (npoints);

    for (int i=0; i<npoints; ++i) {
        double wl, wu;
        switch (diag_class) {
            case k1:
                wl = vertex1.avertex().K1.frequencies.get_wupper_b();
                wu = vertex1.avertex().K1.frequencies.get_wupper_b();
                wl *= 2.;
                wu *= 2.;
                break;
            case k2:
                wl = vertex1.avertex().K2.frequencies.get_wlower_f();
                wu = vertex1.avertex().K2.frequencies.get_wupper_f();
                break;
            case k3:
                wl = vertex1.avertex().K3.frequencies.get_wlower_f();
                wu = vertex1.avertex().K3.frequencies.get_wupper_f();
                break;
            default:;
        }
        double vpp = wl + i * (wu - wl) / (npoints - 1);
        //if (diag_class == k1 and not HUBBARD_MODEL) {vertex1.avertex().K1.frequencies.get_freqs_w(vpp, i);}
        freqs[i] = vpp;
    }

    std::string filename_prefix = "";
    if (HUBBARD_MODEL) filename_prefix = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/integrands/";
    save_integrand(freqs, filename_prefix);

}

template<K_class diag_class, char channel, int spin, typename Q, vertexType symmetry_left, vertexType symmetry_right, class Bubble_Object,typename return_type>
void Integrand<diag_class,channel, spin, Q, symmetry_left, symmetry_right, Bubble_Object,return_type>::save_integrand(const rvec& freqs, const std::string& filename_prefix) const {
    /// evaluate the integrand on frequency points in freqs
    int npoints = freqs.size();

    rvec integrand_re (npoints);
    rvec integrand_im (npoints);
    rvec Pival_re (npoints);
    rvec Pival_im (npoints);
    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> integrand_vals(VECTORIZED_INTEGRATION ? 16 : 1, npoints);
    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> Pivals(KELDYSH?16:1,npoints);
    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> vertex_vals1(VECTORIZED_INTEGRATION ? 16 : 1, npoints);
    Eigen::Matrix<Q,Eigen::Dynamic,Eigen::Dynamic> vertex_vals2(VECTORIZED_INTEGRATION ? 16 : 1,npoints);

    get_integrand_vals(freqs, integrand_vals, Pivals, vertex_vals1, vertex_vals2);

    std::string filename = "";
    if (not HUBBARD_MODEL) filename += data_dir;
    filename += filename_prefix+"integrand_K" + (diag_class == k1 ? "1" : (diag_class == k2 ? "2" : (diag_class == k3 ? "3" : "2b")));
    filename += channel;
    filename += "_i0=" + std::to_string(i0_symmred)
                + "_i2=" + std::to_string(i2)
                + "_w=" + std::to_string(w);
    if (diag_class == k2) {filename += "_v=" + std::to_string(v);}
    else if (diag_class == k3) {filename += "_vp=" + std::to_string(vp);}
    filename += "_i_in=" + std::to_string(i_in);
    filename += + ".h5";

    H5::H5File file(filename, H5F_ACC_TRUNC);
    write_to_hdf(file, "v", freqs, false);
    write_to_hdf(file, "integrand", integrand_vals, false);
    write_to_hdf(file, "Pival", Pivals, false);
    write_to_hdf(file, "vertex1", vertex_vals1, false);
    write_to_hdf(file, "vertex2", vertex_vals2, false);
}


#endif //KELDYSH_MFRG_INTEGRAND_HPP
