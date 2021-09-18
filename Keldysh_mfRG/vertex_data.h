#ifndef FPP_MFRG_VERTEX_DATA_H
#define FPP_MFRG_VERTEX_DATA_H

/**
 * This header contains a class which is responsible for saving and retrieving vertex data
 */


#include "data_structures.h"          // real/complex vector classes
#include "parameters/master_parameters.h"               // system parameters (lengths of vectors etc.)
#include "symmetries/Keldysh_symmetries.h"       // transformations on Keldysh indices
#include "symmetries/internal_symmetries.h"      // symmetry transformations for internal indices (momentum etc.), currently trivial
#include "interpolations/vertex_interpolations.h"           // frequency interpolations for vertices
#include "symmetries/symmetry_transformations.h" // symmetry transformations of frequencies
#include "symmetries/symmetry_table.h"           // table containing information when to apply which symmetry transformations
#include "grids/momentum_grid.h"            // functionality for the internal structure of the Hubbard model

template <typename Q> class fullvert; // forward declaration of fullvert
template <typename Q> class State; // forward declaration of fullvert
template <typename Q, template <typename> class symmetry_type> class GeneralVertex;
template <typename Q>class symmetric;
template <typename Q>using Vertex = GeneralVertex<Q, symmetric>;
class Buffer;


template <typename Q>
class vertexDataContainer{
    friend void check_Kramers_Kronig(std::string filename);
    friend void test_PT4(double Lambda, bool write_flag);
    template<typename T> friend void susceptibilities_postprocessing(Vertex<T>& chi, Vertex<T>& chi_diff, const State<T>& state, double Lambda);
    template<typename T> friend void check_FDTs(const State<T>& state, bool verbose);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template <K_class k, typename T> friend class Interpolate;

private:
    vec<Q> empty_K2() { // for pure K1-calculation no memory should be allocated unnecessarily for K2
        if (MAX_DIAG_CLASS >= 2) return vec<Q> (nK_K2 * nw2 * nv2 * n_in);  // data points of K2;
        else                     return vec<Q> (0);                         // empty vector, never used in calculations
    }
    vec<Q> empty_K3() { // for  K2-calculation no memory should be allocated unnecessarily for K3
        if (MAX_DIAG_CLASS >= 3) return vec<Q> (nK_K3 * nw3 * nv3 * nv3 * n_in);    // data points of K3  // data points of K2;
        else                     return vec<Q> (0);                                 // empty vector, never used in calculations
    }
    VertexFrequencyGrid frequencies;    // frequency grid

protected:
    /** K1-functionality */
    vec<Q> K1 = vec<Q> (nK_K1 * nw1 * n_in);  // data points of K1
    /** K2 functionality */
    vec<Q> K2 = empty_K2();
    /** K3 functionality */
    vec<Q> K3 = empty_K3();

    void K2_convert2internalFreqs(double& w, double& v) const;
    void K2_convert2naturalFreqs(double& w, double& v) const;

public:



    explicit vertexDataContainer(double Lambda) : frequencies(Lambda) { };

    auto get_VertexFreqGrid() const -> VertexFrequencyGrid;
    void set_VertexFreqGrid(const VertexFrequencyGrid &frequencyGrid);

    /// Member functions for accessing the reducible vertex in channel r at arbitrary frequencies ///
    /// by interpolating stored data, in all possible channel-dependent frequency representations ///

    /// Member functions for accessing/setting values of the vector K1 ///

    /** Return the value of the vector K1 at index i. */
    auto K1_acc(int i) const -> Q;

    /** Set the value of the vector K1 at index i to "value". */
    void K1_direct_set(int i, Q value);

    /** Set the value of the vector K1 at Keldysh index iK, frequency index iw,
     * internal structure index i_in to "value". */
    void K1_setvert(int iK, int iw, int i_in, Q value);

    /** Add "value" to the value of the vector K1 at Keldysh index iK, frequency index iw,
     * internal structure index i_in. */
    void K1_addvert(int iK, int iw, int i_in, Q value);

    /** Return the value of the vector K1 at Keldysh index iK, frequency index iw,
     * internal structure index i_in. */
    auto K1_val(int iK, int iw, int i_in) const -> Q;

    void K1_add(vec<Q> summand);
    double K1_get_wlower() const;
    double K1_get_wupper() const;
    auto K1_get_freqGrid() const -> FrequencyGrid;
    void K1_get_freq_w(double& w, int i) const;
    auto K2_get_correction_MFfiniteT(int iw) const -> double;

    /// Member functions for accessing/setting values of the vector K2 ///

    /** Return the value of the vector K2 at index i. */
    auto K2_acc(int i) const -> Q;

    /** Set the value of the vector K2 at index i to "value". */
    void K2_direct_set(int i, Q value);

    /** Set the value of the vector K2 at Keldysh index iK, frequency indices iw, iv,
     * internal structure index i_in to "value". */
    void K2_setvert(int iK, int iw, int iv, int i_in, Q value);

    /** Add "value" to the value of the vector K2 at Keldysh index iK, frequency indices iw, iv,
     * internal structure index i_in. */
    void K2_addvert(int iK, int iw, int iv, int i_in, Q value);

    /** Return the value of the vector K2 at Keldysh index iK, frequency indices iw, iv,
     * internal structure index i_in. */
    auto K2_val(int iK, int iw, int iv, int i_in) const -> Q;

    void K2_add(vec<Q> summand);
    double K2_get_wlower_b() const;
    double K2_get_wupper_b() const;
    double K2_get_wlower_f() const;
    double K2_get_wupper_f() const;
    FrequencyGrid K2_get_freqGrid_b() const;
    FrequencyGrid K2_get_freqGrid_f() const;
    void K2_get_freqs_w(double& w, double& v, int iw, int iv) const;

    /// Member functions for accessing/setting values of the vector K3 ///

    /** Return the value of the vector K3 at index i. */
    auto K3_acc(int i) const -> Q;

    /** Set the value of the vector K3 at index i to "value". */
    void K3_direct_set(int i, Q value);

    /** Set the value of the vector K3 at Keldysh index iK, frequency indices iw, iv, ivp,
     * internal structure index i_in to "value". */
    void K3_setvert(int iK, int iw, int iv, int ivp, int i_in, Q);

    /** Add "value" to the value of the vector K3 at Keldysh index iK, frequency indices iw, iv, ivp,
     * internal structure index i_in. */
    void K3_addvert(int iK, int iw, int iv, int ivp, int i_in, Q);

    /** Return the value of the vector K3 at Keldysh index iK, frequency indices iw, iv, ivp,
     * internal structure index i_in. */
    auto K3_val(int iK, int iw, int iv, int ivp, int i_in) const -> Q;

    void K3_add(vec<Q> summand);
    double K3_get_wlower_b() const;
    double K3_get_wupper_b() const;
    double K3_get_wlower_f() const;
    double K3_get_wupper_f() const;
    FrequencyGrid K3_get_freqGrid_b() const;
    FrequencyGrid K3_get_freqGrid_f() const;
    void K3_get_freqs_w(double& w, double& v, double& vp, int iw, int iv, int ivp) const;

    auto K3_get_correction_MFfiniteT(int iw) const -> double;



    double get_deriv_maxK1() const;
    double get_deriv_maxK2() const;
    double get_deriv_maxK3() const;

};

/****************************************** MEMBER FUNCTIONS OF THE R-VERTEX ******************************************/

template<typename Q>
auto vertexDataContainer<Q>::get_VertexFreqGrid() const -> VertexFrequencyGrid {
    return frequencies;
}
template<typename Q>
void vertexDataContainer<Q>::set_VertexFreqGrid(const VertexFrequencyGrid& frequencyGrid) {
    frequencies = frequencyGrid;
}

template <typename Q> auto vertexDataContainer<Q>::K1_acc(int i) const -> Q {
    if (i >= 0 && i < K1.size())
        return K1[i];
    else
        print("Error: Tried to access value outside of K1 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<Q>::K1_direct_set(int i, Q value) {
    if (i >= 0 && i < K1.size())
        K1[i] = value;
    else
        print("Error: Tried to access value outside of K1 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<Q>::K1_setvert(int iK, int iw, int i_in, Q value) {
    K1[iK*nw1*n_in + iw*n_in + i_in] = value;
}
template <typename Q> void vertexDataContainer<Q>::K1_addvert(int iK, int iw, int i_in, Q value) {
    K1[iK*nw1*n_in + iw*n_in + i_in] += value;
}
template <typename Q> auto vertexDataContainer<Q>::K1_val(int iK, int iw, int i_in) const -> Q {
    return K1[iK*nw1*n_in + iw*n_in + i_in];
}
template<typename Q>
void vertexDataContainer<Q>::K1_add(vec<Q> summand) {
    K1 += summand;
}
template<typename Q>
double vertexDataContainer<Q>::K1_get_wlower() const {
    return frequencies.b_K1.w_lower;
}
template<typename Q>
double vertexDataContainer<Q>::K1_get_wupper() const {
    return frequencies.b_K1.w_upper;
}
template<typename Q>
auto vertexDataContainer<Q>::K1_get_freqGrid() const -> FrequencyGrid {
    return frequencies.b_K1;
}
template<typename Q>
void vertexDataContainer<Q>::K1_get_freq_w(double& w, const int i) const {
    w = frequencies.b_K1.ws[i];
}




template <typename Q> auto vertexDataContainer<Q>::K2_acc(int i) const -> Q {
    if (i >= 0 && i < K2.size())
        return K2[i];
    else
        print("Error: Tried to access value outside of K2 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<Q>::K2_direct_set(int i, Q value) {
    if (i >= 0 && i < K2.size())
        K2[i] = value;
    else
        print("Error: Tried to access value outside of K2 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<Q>::K2_setvert(int iK, int iw, int iv, int i_in, Q value) {
    K2[iK*nw2*nv2*n_in + iw*nv2*n_in + iv*n_in + i_in] = value;
}
template <typename Q> void vertexDataContainer<Q>::K2_addvert(int iK, int iw, int iv, int i_in, Q value) {
    K2[iK*nw2*nv2*n_in + iw*nv2*n_in + iv*n_in + i_in] += value;
}
template <typename Q> auto vertexDataContainer<Q>::K2_val(int iK, int iw, int iv, int i_in) const -> Q {
    return K2[iK * nw2 * nv2 * n_in + iw * nv2 * n_in + iv * n_in + i_in];
}
template<typename Q>
void vertexDataContainer<Q>::K2_add(vec<Q> summand) {
    K2 += summand;
}
template<typename Q>
double vertexDataContainer<Q>::K2_get_wlower_b() const {
    return frequencies.b_K2.w_lower;
}
template<typename Q>
double vertexDataContainer<Q>::K2_get_wupper_b() const {
    return frequencies.b_K2.w_upper;
}
template<typename Q>
double vertexDataContainer<Q>::K2_get_wlower_f() const {
    return frequencies.f_K2.w_lower;
}
template<typename Q>
double vertexDataContainer<Q>::K2_get_wupper_f() const {
    return frequencies.f_K2.w_upper;
}
template<typename Q>
auto vertexDataContainer<Q>::K2_get_freqGrid_b() const -> FrequencyGrid {
    return frequencies.b_K2;
}
template<typename Q>
auto vertexDataContainer<Q>::K2_get_freqGrid_f() const -> FrequencyGrid {
    return frequencies.f_K2;
}
template<typename Q>
void vertexDataContainer<Q>::K2_get_freqs_w(double &w, double &v, const int iw, const int iv) const {
    w = frequencies.b_K2.ws[iw];
    v = frequencies.f_K2.ws[iv];
    K2_convert2naturalFreqs(w, v);
}
template<typename Q>
auto vertexDataContainer<Q>::K2_get_correction_MFfiniteT(int iw) const -> double {
    return floor2bfreq(frequencies.b_K2.ws[iw] / 2) - ceil2bfreq(frequencies.b_K2.ws[iw] / 2);
}
template<typename Q>
void vertexDataContainer<Q>::K2_convert2internalFreqs(double &w, double &v) const {
    //const double w_tmp = w/2. + v;
    //const double v_tmp = w/2. - v;
    //w = w_tmp;
    //v = v_tmp;
}
template<typename Q>
void vertexDataContainer<Q>::K2_convert2naturalFreqs(double &w, double &v) const {
    //const double w_tmp = w + v;
    //const double v_tmp =(w - v)/2.;
    //w = w_tmp;
    //v = v_tmp;
}


template <typename Q> auto vertexDataContainer<Q>::K3_acc(int i) const -> Q {
    if (i >= 0 && i < K3.size())
        return K3[i];
    else
        print("Error: Tried to access value outside of K3 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<Q>::K3_direct_set(int i, Q value) {
    if (i >= 0 && i < K3.size())
        K3[i] = value;
    else
        print("Error: Tried to access value outside of K3 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<Q>::K3_setvert(int iK, int iw, int iv, int ivp, int i_in, Q value) {
    K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in] = value;
}
template <typename Q> void vertexDataContainer<Q>::K3_addvert(int iK, int iw, int iv, int ivp, int i_in, Q value) {
    K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in] += value;
}
template <typename Q> auto vertexDataContainer<Q>::K3_val(int iK, int iw, int iv, int ivp, int i_in) const -> Q {
    return K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in];
}
template<typename Q>
void vertexDataContainer<Q>::K3_add(vec<Q> summand) {
    K3 += summand;
}
template<typename Q>
double vertexDataContainer<Q>::K3_get_wlower_b() const {
    return frequencies.b_K3.w_lower;
}
template<typename Q>
double vertexDataContainer<Q>::K3_get_wupper_b() const {
    return frequencies.b_K3.w_upper;
}
template<typename Q>
double vertexDataContainer<Q>::K3_get_wlower_f() const {
    return frequencies.f_K3.w_lower;
}
template<typename Q>
double vertexDataContainer<Q>::K3_get_wupper_f() const {
    return frequencies.f_K3.w_upper;
}
template<typename Q>
auto vertexDataContainer<Q>::K3_get_freqGrid_b() const -> FrequencyGrid {
    return frequencies.b_K3;
}
template<typename Q>
auto vertexDataContainer<Q>::K3_get_freqGrid_f() const -> FrequencyGrid {
    return frequencies.f_K3;
}
template<typename Q>
void vertexDataContainer<Q>::K3_get_freqs_w(double &w, double &v, double& vp, const int iw, const int iv, const int ivp) const {
    w = frequencies.b_K3.ws[iw];
    v = frequencies.f_K3.ws[iv];
    vp= frequencies.f_K3.ws[ivp];
}
template<typename Q>
auto vertexDataContainer<Q>::K3_get_correction_MFfiniteT(int iw) const -> double {
    return floor2bfreq(frequencies.b_K3.ws[iw] / 2) - ceil2bfreq(frequencies.b_K3.ws[iw] / 2);
}

template <typename Q> auto vertexDataContainer<Q>::get_deriv_maxK1() const -> double {
    double max_K1 = ::power2(::get_finite_differences(K1)).max_norm();
    return max_K1;

}
template <typename Q> auto vertexDataContainer<Q>::get_deriv_maxK2() const -> double {
    double max_K2 = (::power2(::get_finite_differences<Q,2>(K2, {nBOS2, nFER2}, {0, 1}))
                     + ::power2(::get_finite_differences<Q,2>(K2, {nFER2}, {1, 0}))
    ).max_norm();
    return max_K2;

}

template <typename Q> auto vertexDataContainer<Q>::get_deriv_maxK3() const -> double {
    double max_K3 = (::power2(::get_finite_differences<Q,3>(K3, {nBOS3, nFER3, nFER3}, {0, 1, 2}))
                     + ::power2(::get_finite_differences<Q,3>(K3, {nFER3, nBOS3, nFER3}, {1, 2, 0}))
                     + ::power2(::get_finite_differences<Q,3>(K3, {nFER3, nFER3, nBOS3}, {2, 0, 1}))
    ).max_norm();
    return max_K3;

}



#endif //FPP_MFRG_VERTEX_DATA_H
