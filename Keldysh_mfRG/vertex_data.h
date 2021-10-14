#ifndef FPP_MFRG_VERTEX_DATA_H
#define FPP_MFRG_VERTEX_DATA_H

/**
 * This header contains a class which is responsible for saving and retrieving vertex data
 */

#include "utilities/template_utils.h"
#include "utilities/math_utils.h"
#include "symmetries/Keldysh_symmetries.h"
#include "data_structures.h"          // real/complex vector classes
#include "parameters/master_parameters.h"               // system parameters (lengths of vectors etc.)
#include "parameters/frequency_parameters.h"
#include "grids/frequency_grid.h"            // functionality for the internal structure of the Hubbard model



template <typename Q> class rvert; // forward declaration of rvert
template <typename Q> class fullvert; // forward declaration of fullvert
template <typename Q> class State; // forward declaration of State
template <typename Q, template <typename> class symmetry_type> class GeneralVertex;
template <typename Q>class symmetric;
template <typename Q>using Vertex = GeneralVertex<Q, symmetric>;
class Buffer;

template <typename Q, size_t rank>
class vertexContainerBase {
private:
    template <typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
    size_t flattenIndex(const Types &... i) const {
#if FREQ_PADDING == 0
        return getFlatIndex({static_cast<size_t>(i)...}, dims);
#else
        std::array<size_t,rank> idx = {static_cast<size_t>(i+1)...};
        idx[0] -= 1; idx[rank-1] -= 1;

        size_t flatidx = getFlatIndex(idx, dims);
        assert(flatidx < data.size());
        return flatidx;
#endif
    }

protected:
    std::array<size_t,rank> dims;
    vec<Q> data;

public:
    explicit vertexContainerBase(const std::vector<size_t> dims_in) {
        assert(dims_in.size() == rank);
        for (size_t i = 0; i < rank; i++) dims[i] = dims_in[i];
    };
    template <typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
    explicit vertexContainerBase(const Types &... dims) : vertexContainerBase(std::vector<size_t>({static_cast<size_t>(dims)...})) {};
    vertexContainerBase(const size_t dims_in[rank], const vec<Q> &data_in) : data(data_in) {};

    void reserve() { data = vec<Q>(getFlatSize<rank>(dims)); }

    Q acc(const size_t flatIndex) const {return data[flatIndex];}
    void direct_set(const size_t flatIndex, Q value) {assert(flatIndex < data.size()); data[flatIndex] = value;}


    template <typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
    Q val(const Types &... i) const {return acc(flattenIndex(i...));}
    template <typename... Types
            ,typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true
                    >
    void setvert(const Q value, const Types &... i) {data[flattenIndex(i...)] = value;
    }

    vec<Q> get_vec() const {return data;}
    /*
    vec<Q> add_padding(const vec<Q> &data_in) {
        std::array<size_t,rank> dims_no_padding = dims;
        for (size_t i = 1; i < rank-1; i++) dims_no_padding[i] -= 2*FREQ_PADDING;
        assert(getFlatSize<rank>(dims_no_padding) == data_in.size());

        vec<Q> data_new (getFlatSize<rank>(dims));
        for (size_t i = 0; i < data_in.size(); i++) {
            std::array<size_t,rank> multIndex;
            getMultIndex(multIndex, i, dims_no_padding);
            for (size_t j = 0; j < rank; j++) multIndex[j] += FREQ_PADDING;
            data_new[getFlatIndex(multIndex, dims)] = data_in[i];
        }
        return data_new;
    }
     */
    void set_vec(const vec<Q> &data_in) {assert(data.size() == data_in.size()); data = data_in;}
    void add_vec(const vec<Q> &summand) {
#if FREQ_PADDING == 0
        data += summand;
#else
        std::array<size_t,rank> dims_no_padding = dims;
        for (size_t i = 1; i < rank-1; i++) dims_no_padding[i] -= 2*FREQ_PADDING;
        assert(getFlatSize<rank>(dims_no_padding) == summand.size());

        //vec<Q> data_new (getFlatSize<rank>(dims));
        for (size_t i = 0; i < summand.size(); i++) {
            std::array<size_t,rank> multIndex;
            getMultIndex(multIndex, i, dims_no_padding);
            for (size_t j = 1; j < rank-1; j++) multIndex[j] += FREQ_PADDING;
            size_t idx = getFlatIndex(multIndex, dims);
            assert(idx < data.size());
            data[idx] += summand[i];
        }
#endif
    }


};

template <K_class k, typename Q>
class vertexDataContainer{};

template<typename Q>
class vertexDataContainer<k1, Q> : public vertexContainerBase<Q,3>{
    template <K_class k, typename T> friend class UpdateGrid;
    template<typename T> friend class CostFullvert_Wscale_b_K1;
    friend void check_Kramers_Kronig(std::string filename);
    friend void test_PT4(double Lambda, bool write_flag);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template<typename T> friend void susceptibilities_postprocessing(Vertex<T>& chi, Vertex<T>& chi_diff, const State<T>& state, double Lambda);
    template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator+= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator- (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator-= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator* (rvert<T> lhs, const double& alpha);
    template<typename T> friend rvert<T> rvert<T>::operator*= (double alpha);


protected:
    //vec<Q> K1 = vec<Q> (nK_K1 * nw1 * n_in);  // data points of K1

    //size_t dims[3] = {nK_K1, nBOS, n_in};

    VertexFrequencyGrid<k1> frequencies_K1;    // frequency grid
public:
    //std::array<size_t,3> dimsK1 = {nK_K1, nBOS, n_in};

    explicit vertexDataContainer(double Lambda) : frequencies_K1(Lambda), vertexContainerBase<Q,3>(dimsK1) { };

    auto K1_get_VertexFreqGrid() const -> const VertexFrequencyGrid<k1>&;
    void K1_set_VertexFreqGrid(const VertexFrequencyGrid<k1> &frequencyGrid);


    /** K1-functionality */
    /// Member functions for accessing/setting values of the vector K1 ///

    const double& K1_get_wlower() const;
    const double& K1_get_wupper() const;
    auto K1_get_freqGrid() const -> const FrequencyGrid&;
    void K1_get_freq_w(double& w, int i) const;
    const double& K1_get_tlower_aux() const;
    const double& K1_get_tupper_aux() const;
    void K1_get_freq_aux(double& w, int i) const;
    auto K1_gridtransf(double w) const -> double;
    auto K1_gridtransf_inv(double w) const -> double;


    vec<Q> get_deriv_K1_x() const;
    double get_deriv_maxK1() const;
    double get_curvature_maxK1() const;

    double analyze_tails_K1() const;

    auto shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k1>;
    //void  findBestFreqGrid(double Lambda);
};

template<typename Q>
class vertexDataContainer<k2, Q>: public vertexContainerBase<Q,4> {
    template <K_class k, typename T> friend class UpdateGrid;
    template<typename T> friend class CostFullvert_Wscale_b_K2;
    template<typename T> friend class CostFullvert_Wscale_f_K2;
    friend void test_PT4(double Lambda, bool write_flag);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template<typename T> friend void check_FDTs(const State<T>& state, bool verbose);
    template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator+= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator- (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator-= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator* (rvert<T> lhs, const double& alpha);
    template<typename T> friend rvert<T> rvert<T>::operator*= (double alpha);

private:
    /*
    vec<Q> empty_K2() { // for pure K1-calculation no memory should be allocated unnecessarily for K2
        if (MAX_DIAG_CLASS >= 2) return vec<Q> (nK_K2 * nw2 * nv2 * n_in);  // data points of K2;
        else                     return vec<Q> (0);                         // empty vector, never used in calculations
    }
    */

protected:
    //vec<Q> K2 = empty_K2();
    //size_t dims[4] = {nK_K2, nBOS2, nFER2, n_in};
    VertexFrequencyGrid<k2> frequencies_K2;    // frequency grid
public:
    //std::array<size_t,4> dimsK2 = {nK_K2, nBOS2, nFER2, n_in};
    explicit vertexDataContainer(double Lambda) : frequencies_K2(Lambda), vertexContainerBase<Q,4>(dimsK2) { };


    /** K2 functionality */

    /// Member functions for accessing the reducible vertex in channel r at arbitrary frequencies ///
    /// by interpolating stored data, in all possible channel-dependent frequency representations ///

    /// Member functions for accessing/setting values of the vector K2 ///

    auto K2_get_VertexFreqGrid() const -> const VertexFrequencyGrid<k2>&;
    void K2_set_VertexFreqGrid(const VertexFrequencyGrid<k2> &frequencyGrid);

    /** Return the value of the vector K2 at index i. */
    //auto K2_acc(int i) const -> Q;

    /** Set the value of the vector K2 at index i to "value". */
    //void K2_direct_set(int i, Q value);

    /** Set the value of the vector K2 at Keldysh index iK, frequency indices iw, iv,
     * internal structure index i_in to "value". */
    //void K2_setvert(int iK, int iw, int iv, int i_in, Q value);

    /** Add "value" to the value of the vector K2 at Keldysh index iK, frequency indices iw, iv,
     * internal structure index i_in. */
    //void K2_addvert(int iK, int iw, int iv, int i_in, Q value);

    /** Return the value of the vector K2 at Keldysh index iK, frequency indices iw, iv,
     * internal structure index i_in. */
    //auto val(int iK, int iw, int iv, int i_in) const -> Q;
    //auto get_K2() const -> vec<Q>;
    //void set_K2(vec<Q> data);

    //void K2_add(vec<Q> summand);
    const double& K2_get_wlower_b() const;
    const double& K2_get_wupper_b() const;
    const double& K2_get_wlower_f() const;
    const double& K2_get_wupper_f() const;
    const FrequencyGrid& K2_get_freqGrid_b() const;
    const FrequencyGrid& K2_get_freqGrid_f() const;
    void K2_get_freqs_w(double& w, double& v, int iw, int iv) const;
    const double& K2_get_tlower_b_aux() const;
    const double& K2_get_tupper_b_aux() const;
    const double& K2_get_tlower_f_aux() const;
    const double& K2_get_tupper_f_aux() const;
    void K2_get_freqs_aux(double& w, double& v, int iw, int iv) const;
    auto K2_gridtransf_b(double w) const -> double;
    auto K2_gridtransf_f(double w) const -> double;
    auto K2_gridtransf_inv_b(double w) const -> double;
    auto K2_gridtransf_inv_f(double w) const -> double;


    void K2_convert2internalFreqs(double& w, double& v) const;
    void K2_convert2naturalFreqs(double& w, double& v) const;

    auto K2_get_correction_MFfiniteT(int iw) const -> double;


    vec<Q> get_deriv_K2_x () const;
    vec<Q> get_deriv_K2_y () const;
    vec<Q> get_deriv_K2_xy() const;
    vec<Q> get_deriv_K2_xx() const;
    vec<Q> get_deriv_K2_yy() const;
    double get_deriv_maxK2() const;
    auto get_curvature_maxK2() const -> double;

    double analyze_tails_K2_x() const;

    double analyze_tails_K2_y() const;
    auto shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k2>;

};

template <typename Q>
class vertexDataContainer<k3, Q>: public vertexContainerBase<Q,5> {
    template <K_class k, typename T> friend class UpdateGrid;
    template<typename T> friend class CostFullvert_Wscale_b_K3;
    template<typename T> friend class CostFullvert_Wscale_f_K3;
    friend void test_PT4(double Lambda, bool write_flag);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template<typename T> friend void check_FDTs(const State<T>& state, bool verbose);
    template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator+= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator- (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator-= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator* (rvert<T> lhs, const double& alpha);
    template<typename T> friend rvert<T> rvert<T>::operator*= (double alpha);

private:
    /*vec<Q> empty_K3() { // for  K2-calculation no memory should be allocated unnecessarily for K3
        if (MAX_DIAG_CLASS >= 3) return vec<Q> (nK_K3 * nw3 * nv3 * nv3 * n_in);    // data points of K3  // data points of K2;
        else                     return vec<Q> (0);                                 // empty vector, never used in calculations
    }*/

protected:
    //vec<Q> K3 = empty_K3();
    VertexFrequencyGrid<k3> frequencies_K3;    // frequency grid


public:
    //std::array<size_t,5> dimsK3 = {nK_K3, nBOS3, nFER3, nFER3, n_in};

    explicit vertexDataContainer(double Lambda) : frequencies_K3(Lambda), vertexContainerBase<Q,5>(dimsK3) { };


    /** K3 functionality */



    /// Member functions for accessing/setting values of the vector K3 ///

    auto K3_get_VertexFreqGrid() const -> const VertexFrequencyGrid<k3>&;
    void K3_set_VertexFreqGrid(const VertexFrequencyGrid<k3> &frequencyGrid);

    /** Return the value of the vector K3 at index i. */
    //auto K3_acc(int i) const -> Q;
    //auto get_K3() const -> vec<Q>;
    //void set_K3(vec<Q> data);

    /** Set the value of the vector K3 at index i to "value". */
    //void K3_direct_set(int i, Q value);

    /** Set the value of the vector K3 at Keldysh index iK, frequency indices iw, iv, ivp,
     * internal structure index i_in to "value". */
    //void K3_setvert(int iK, int iw, int iv, int ivp, int i_in, Q);

    /** Add "value" to the value of the vector K3 at Keldysh index iK, frequency indices iw, iv, ivp,
     * internal structure index i_in. */
    //void K3_addvert(int iK, int iw, int iv, int ivp, int i_in, Q);

    /** Return the value of the vector K3 at Keldysh index iK, frequency indices iw, iv, ivp,
     * internal structure index i_in. */
    //auto val(int iK, int iw, int iv, int ivp, int i_in) const -> Q;

    //void K3_add(vec<Q> summand);
    const double& K3_get_wlower_b() const;
    const double& K3_get_wupper_b() const;
    const double& K3_get_wlower_f() const;
    const double& K3_get_wupper_f() const;
    const FrequencyGrid& K3_get_freqGrid_b() const;
    const FrequencyGrid& K3_get_freqGrid_f() const;
    void K3_get_freqs_w(double& w, double& v, double& vp, int iw, int iv, int ivp) const;
    const double& K3_get_tlower_b_aux() const;
    const double& K3_get_tupper_b_aux() const;
    const double& K3_get_tlower_f_aux() const;
    const double& K3_get_tupper_f_aux() const;
    void K3_get_freqs_aux(double& w, double& v, double& vp, int iw, int iv, int ivp) const;
    auto K3_gridtransf_b(double w) const -> double;
    auto K3_gridtransf_f(double w) const -> double;
    auto K3_gridtransf_inv_b(double w) const -> double;
    auto K3_gridtransf_inv_f(double w) const -> double;

    auto K3_get_correction_MFfiniteT(int iw) const -> double;


    vec<Q> get_deriv_K3_x  () const;
    vec<Q> get_deriv_K3_y  () const;
    vec<Q> get_deriv_K3_z  () const;
    vec<Q> get_deriv_K3_xy () const;
    vec<Q> get_deriv_K3_xz () const;
    vec<Q> get_deriv_K3_yz () const;
    vec<Q> get_deriv_K3_xx () const;
    vec<Q> get_deriv_K3_yy () const;
    vec<Q> get_deriv_K3_zz () const;
    vec<Q> get_deriv_K3_xyz() const;
    double get_deriv_maxK3() const;
    auto get_curvature_maxK3() const -> double;
    double analyze_tails_K3_x() const;
    double analyze_tails_K3_y() const;
    double analyze_tails_K3_z() const;

    auto shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k3>;

};

/************************************ MEMBER FUNCTIONS OF THE VERTEX Data Container************************************/
template<typename Q>
auto vertexDataContainer<k1,Q>::K1_get_VertexFreqGrid() const -> const VertexFrequencyGrid<k1>& {
    return frequencies_K1;
}
template<typename Q>
void vertexDataContainer<k1,Q>::K1_set_VertexFreqGrid(const VertexFrequencyGrid<k1>& frequencyGrid) {
    frequencies_K1 = frequencyGrid;
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_get_VertexFreqGrid() const -> const VertexFrequencyGrid<k2>& {
    return frequencies_K2;
}
template<typename Q>
void vertexDataContainer<k2,Q>::K2_set_VertexFreqGrid(const VertexFrequencyGrid<k2>& frequencyGrid) {
    frequencies_K2 = frequencyGrid;
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_VertexFreqGrid() const -> const VertexFrequencyGrid<k3> & {
    return frequencies_K3;
}
template<typename Q>
void vertexDataContainer<k3,Q>::K3_set_VertexFreqGrid(const VertexFrequencyGrid<k3>& frequencyGrid) {
    frequencies_K3 = frequencyGrid;
}
/*
template <typename Q> auto vertexDataContainer<k1,Q>::K1_acc(int i) const -> Q {
    if (i >= 0 && i < K1.size())
        return K1[i];
    else
        print("Error: Tried to access value outside of K1 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<k1,Q>::K1_direct_set(int i, Q value) {
    if (i >= 0 && i < K1.size())
        K1[i] = value;
    else
        print("Error: Tried to access value outside of K1 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<k1,Q>::K1_setvert(int iK, int iw, int i_in, Q value) {
    K1[iK*nw1*n_in + iw*n_in + i_in] = value;
}
template <typename Q> void vertexDataContainer<k1,Q>::K1_addvert(int iK, int iw, int i_in, Q value) {
    K1[iK*nw1*n_in + iw*n_in + i_in] += value;
}
template <typename Q> auto vertexDataContainer<k1,Q>::val(int iK, int iw, int i_in) const -> Q {
    return K1[iK*nw1*n_in + iw*n_in + i_in];
}
template<typename Q>
void vertexDataContainer<k1,Q>::K1_add(vec<Q> summand) {
    K1 += summand;
}
template<typename Q> auto vertexDataContainer<k1,Q>::get_K1() const -> vec<Q> {
    return K1;
}
*/
/*
template<typename Q> void vertexDataContainer<k1,Q>::set_K1(const vec<Q> data) {
    return K1 = data;
}*/
template<typename Q>
const double& vertexDataContainer<k1,Q>::K1_get_wlower() const {
    return frequencies_K1.b.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k1,Q>::K1_get_wupper() const {
    return frequencies_K1.b.w_upper;
}
template<typename Q>
auto vertexDataContainer<k1,Q>::K1_get_freqGrid() const -> const FrequencyGrid& {
    return frequencies_K1.b;
}
template<typename Q>
void vertexDataContainer<k1,Q>::K1_get_freq_w(double& w, const int i) const {
    w = frequencies_K1.b.get_ws(i);
}
template<typename Q>
const double& vertexDataContainer<k1,Q>::K1_get_tlower_aux() const {
    return frequencies_K1.b.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k1,Q>::K1_get_tupper_aux() const {
    return frequencies_K1.b.t_upper;
}
template<typename Q>
void vertexDataContainer<k1,Q>::K1_get_freq_aux(double& w, const int i) const {
    w = frequencies_K1.b.get_ts(i);
}
template<typename Q>
auto vertexDataContainer<k1,Q>::K1_gridtransf(double w) const -> double {
    return frequencies_K1.b.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k1,Q>::K1_gridtransf_inv(double w) const -> double {
    return frequencies_K1.b.grid_transf_inv(w);
}



template <typename Q> auto vertexDataContainer<k1,Q>::get_deriv_K1_x() const -> vec<Q> {
    vec<Q> result = ::partial_deriv<Q,3>(vertexContainerBase<Q,3>::data, frequencies_K1.b.ts, vertexContainerBase<Q,3>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k1,Q>::get_deriv_maxK1() const -> double {
    double max_K1 = ::power2(get_deriv_K1_x()).max_norm();
    return max_K1;
}
template <typename Q> auto vertexDataContainer<k1,Q>::get_curvature_maxK1() const -> double {
    double max_K1 = ::power2( ::partial_deriv<Q,3>(get_deriv_K1_x(), frequencies_K1.b.ts, vertexContainerBase<Q,3>::dims, 1)).max_norm();
    return max_K1;
}


template <typename Q> double vertexDataContainer<k1,Q>::analyze_tails_K1() const {
    double maxabs_K1_total = vertexContainerBase<Q,3>::data.max_norm();
    vec<double> maxabsK1_along_w = maxabs(vertexContainerBase<Q,3>::data, vertexContainerBase<Q,3>::dims, 1);

    return maxabsK1_along_w[FREQ_PADDING] / maxabs_K1_total;
}

template <typename Q> auto vertexDataContainer<k1,Q>::shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k1> {
    vec<double> maxabsK1_along_w = maxabs(vertexContainerBase<Q,3>::data, vertexContainerBase<Q,3>::dims, 1);
    double maxmax = maxabsK1_along_w.max_norm();

    VertexFrequencyGrid<k1> frequencies_new = frequencies_K1;

    //const double rel_tail_threshold = 1e-4;
    size_t index = 0;
    while (true) {
        if (maxabsK1_along_w[index] > rel_tail_threshold * maxmax) break;
        index++;
    }
    frequencies_new.b.set_w_upper(std::abs(frequencies_K1.b.get_ws(index)));
    frequencies_new.b.initialize_grid();

    return frequencies_new;
}


/*
template <typename Q> auto vertexDataContainer<k2,Q>::K2_acc(int i) const -> Q {
    if (i >= 0 && i < K2.size())
        return K2[i];
    else
        print("Error: Tried to access value outside of K2 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<k2,Q>::K2_direct_set(int i, Q value) {
    if (i >= 0 && i < K2.size())
        K2[i] = value;
    else
        print("Error: Tried to access value outside of K2 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<k2,Q>::K2_setvert(int iK, int iw, int iv, int i_in, Q value) {
    K2[iK*nw2*nv2*n_in + iw*nv2*n_in + iv*n_in + i_in] = value;
}
template <typename Q> void vertexDataContainer<k2,Q>::K2_addvert(int iK, int iw, int iv, int i_in, Q value) {
    K2[iK*nw2*nv2*n_in + iw*nv2*n_in + iv*n_in + i_in] += value;
}
template <typename Q> auto vertexDataContainer<k2,Q>::val(int iK, int iw, int iv, int i_in) const -> Q {
    return K2[iK * nw2 * nv2 * n_in + iw * nv2 * n_in + iv * n_in + i_in];
}
template<typename Q>
void vertexDataContainer<k2,Q>::K2_add(vec<Q> summand) {
    K2 += summand;
}
template<typename Q> auto vertexDataContainer<k2,Q>::get_K2() const -> vec<Q> {
    return K2;
}
*/
/*
template<typename Q> void vertexDataContainer<k2,Q>::set_K2(const vec<Q> data) {
    return K2 = data;
}*/
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_wlower_b() const {
    return frequencies_K2.b.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_wupper_b() const {
    return frequencies_K2.b.w_upper;
}
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_wlower_f() const {
    return frequencies_K2.f.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_wupper_f() const {
    return frequencies_K2.f.w_upper;
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_get_freqGrid_b() const -> const FrequencyGrid& {
    return frequencies_K2.b;
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_get_freqGrid_f() const -> const FrequencyGrid& {
    return frequencies_K2.f;
}
template<typename Q>
void vertexDataContainer<k2,Q>::K2_get_freqs_w(double &w, double &v, const int iw, const int iv) const {
    w = frequencies_K2.b.get_ws(iw);
    v = frequencies_K2.f.get_ws(iv);
    K2_convert2naturalFreqs(w, v);
}

template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_tlower_b_aux() const {
    return frequencies_K2.b.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_tupper_b_aux() const {
    return frequencies_K2.b.t_upper;
}
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_tlower_f_aux() const {
    return frequencies_K2.f.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k2,Q>::K2_get_tupper_f_aux() const {
    return frequencies_K2.f.t_upper;
}
template<typename Q>
void vertexDataContainer<k2,Q>::K2_get_freqs_aux(double &w, double &v, const int iw, const int iv) const {
    w = frequencies_K2.b.get_ts(iw);
    v = frequencies_K2.f.get_ts(iv);
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_gridtransf_b(double w) const -> double {
    return frequencies_K2.b.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_gridtransf_f(double w) const -> double {
    return frequencies_K2.f.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_gridtransf_inv_b(double w) const -> double {
    return frequencies_K2.b.grid_transf_inv(w);
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_gridtransf_inv_f(double w) const -> double {
    return frequencies_K2.f.grid_transf_inv(w);
}
template<typename Q>
auto vertexDataContainer<k2,Q>::K2_get_correction_MFfiniteT(int iw) const -> double {
    if (not KELDYSH and not ZERO_T)
        return floor2bfreq(frequencies_K2.b.get_ws(iw) / 2) - ceil2bfreq(frequencies_K2.b.get_ws(iw) / 2);
    else return 0.;
}
template<typename Q>
void vertexDataContainer<k2,Q>::K2_convert2internalFreqs(double &w, double &v) const {
    //const double w_tmp = w/2. + v;
    //const double v_tmp = w/2. - v;
    //w = w_tmp;
    //v = v_tmp;
}
template<typename Q>
void vertexDataContainer<k2,Q>::K2_convert2naturalFreqs(double &w, double &v) const {
    //const double w_tmp = w + v;
    //const double v_tmp =(w - v)/2.;
    //w = w_tmp;
    //v = v_tmp;
}




template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_x() const -> vec<Q> {
    vec<Q> result = ::partial_deriv<Q,4>(vertexContainerBase<Q,4>::data, frequencies_K2.b.ts, vertexContainerBase<Q,4>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_y() const -> vec<Q> {
    vec<Q> result = ::partial_deriv<Q,4>(vertexContainerBase<Q,4>::data, frequencies_K2.f.ts, vertexContainerBase<Q,4>::dims, 2);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_xy() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,4>(vertexContainerBase<Q,4>::data, frequencies_K2.f.ts, vertexContainerBase<Q,4>::dims, 2);
    vec<Q> result       = ::partial_deriv<Q,4>(inter_result, frequencies_K2.b.ts, vertexContainerBase<Q,4>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_xx() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,4>(vertexContainerBase<Q,4>::data, frequencies_K2.b.ts, vertexContainerBase<Q,4>::dims, 1);
    vec<Q> result       = ::partial_deriv<Q,4>(inter_result, frequencies_K2.b.ts, vertexContainerBase<Q,4>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_yy() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,4>(vertexContainerBase<Q,4>::data, frequencies_K2.f.ts, vertexContainerBase<Q,4>::dims, 2);
    vec<Q> result       = ::partial_deriv<Q,4>(inter_result, frequencies_K2.f.ts, vertexContainerBase<Q,4>::dims, 2);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_maxK2() const -> double {
    double max_K2 = (     ::power2(get_deriv_K2_x())
                        + ::power2(get_deriv_K2_y())
                    ).max_norm();
    return max_K2;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_curvature_maxK2() const -> double {
    double max_K1 = (    ::power2( get_deriv_K2_xx())
                        +::power2( get_deriv_K2_yy())
                        +::power2( get_deriv_K2_xy())
                    ).max_norm();
    return max_K1;
}

template <typename Q> double vertexDataContainer<k2,Q>::analyze_tails_K2_x() const {
    double maxabs_K2_total = vertexContainerBase<Q,4>::data.max_norm();
    vec<double> maxabsK2_along_w = maxabs(vertexContainerBase<Q,4>::data, vertexContainerBase<Q,4>::dims, 1);

    return maxabsK2_along_w[FREQ_PADDING] / maxabs_K2_total;
}
template <typename Q> double vertexDataContainer<k2,Q>::analyze_tails_K2_y() const {
    double maxabs_K2_total = vertexContainerBase<Q,4>::data.max_norm();
    vec<double> maxabsK2_along_v = maxabs(vertexContainerBase<Q,4>::data, vertexContainerBase<Q,4>::dims, 2);

    return maxabsK2_along_v[FREQ_PADDING] / maxabs_K2_total;
}

template <typename Q> auto vertexDataContainer<k2,Q>::shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k2> {

    VertexFrequencyGrid<k2> frequencies_new = frequencies_K2;

    vec<double> maxabsK2_along_w = maxabs(vertexContainerBase<Q,4>::data, vertexContainerBase<Q,4>::dims, 1);
    double maxmax = maxabsK2_along_w.max_norm();

    //const double rel_tail_threshold = 1e-4;
    size_t index = 0;
    while (true) {
        if (maxabsK2_along_w[index] > rel_tail_threshold * maxmax) break;
        index++;
    }
    frequencies_new.b.set_w_upper(std::abs(frequencies_K2.b.get_ws(index)));
    frequencies_new.b.initialize_grid();



    vec<double> maxabsK2_along_v = maxabs(vertexContainerBase<Q,4>::data, vertexContainerBase<Q,4>::dims, 2);

    //const double rel_tail_threshold = 1e-4;
    index = 0;
    while (true) {
        if (maxabsK2_along_v[index] > rel_tail_threshold * maxmax) break;
        index++;
    }
    frequencies_new.f.set_w_upper(std::abs(frequencies_K2.f.get_ws(index)));
    frequencies_new.f.initialize_grid();

    return frequencies_new;
}

/*
template <typename Q> auto vertexDataContainer<k3,Q>::K3_acc(int i) const -> Q {
    if (i >= 0 && i < K3.size())
        return K3[i];
    else
        print("Error: Tried to access value outside of K3 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<k3,Q>::K3_direct_set(int i, Q value) {
    if (i >= 0 && i < K3.size())
        K3[i] = value;
    else
        print("Error: Tried to access value outside of K3 vertex in a-channel", true);
}
template <typename Q> void vertexDataContainer<k3,Q>::K3_setvert(int iK, int iw, int iv, int ivp, int i_in, Q value) {
    K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in] = value;
}
template <typename Q> void vertexDataContainer<k3,Q>::K3_addvert(int iK, int iw, int iv, int ivp, int i_in, Q value) {
    K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in] += value;
}
template <typename Q> auto vertexDataContainer<k3,Q>::val(int iK, int iw, int iv, int ivp, int i_in) const -> Q {
    return K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in];
}
template<typename Q>
void vertexDataContainer<k3,Q>::K3_add(vec<Q> summand) {
    K3 += summand;
}
template<typename Q> auto vertexDataContainer<k3,Q>::get_K3() const -> vec<Q> {
    return K3;
}
*/
/*
template<typename Q> void vertexDataContainer<k3,Q>::set_K3(const vec<Q> data) {
    return K3 = data;
}*/
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wlower_b() const {
    return frequencies_K3.b.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wupper_b() const {
    return frequencies_K3.b.w_upper;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wlower_f() const {
    return frequencies_K3.f.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wupper_f() const {
    return frequencies_K3.f.w_upper;
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_freqGrid_b() const -> const FrequencyGrid& {
    return frequencies_K3.b;
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_freqGrid_f() const -> const FrequencyGrid& {
    return frequencies_K3.f;
}
template<typename Q>
void vertexDataContainer<k3,Q>::K3_get_freqs_w(double &w, double &v, double& vp, const int iw, const int iv, const int ivp) const {
    w = frequencies_K3.b.get_ws(iw);
    v = frequencies_K3.f.get_ws(iv);
    vp= frequencies_K3.f.get_ws(ivp);
}

template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tlower_b_aux() const {
    return frequencies_K3.b.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tupper_b_aux() const {
    return frequencies_K3.b.t_upper;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tlower_f_aux() const {
    return frequencies_K3.f.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tupper_f_aux() const {
    return frequencies_K3.f.t_upper;
}
template<typename Q>
void vertexDataContainer<k3,Q>::K3_get_freqs_aux(double &w, double &v, double& vp, const int iw, const int iv, const int ivp) const {
    w = frequencies_K3.b.get_ts(iw);
    v = frequencies_K3.f.get_ts(iv);
    vp= frequencies_K3.f.get_ts(ivp);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_b(double w) const -> double {
    return frequencies_K3.b.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_f(double w) const -> double {
    return frequencies_K3.f.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_inv_b(double w) const -> double {
    return frequencies_K3.b.grid_transf_inv(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_inv_f(double w) const -> double {
    return frequencies_K3.f.grid_transf_inv(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_correction_MFfiniteT(int iw) const -> double {
    if (not KELDYSH and not ZERO_T)
    return floor2bfreq(frequencies_K3.b.get_ws(iw) / 2) - ceil2bfreq(frequencies_K3.b.get_ws(iw) / 2);
    else return 0.;
}


template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_x() const -> vec<Q> {
    vec<Q> result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_y() const -> vec<Q> {
    vec<Q> result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_z() const -> vec<Q> {
    vec<Q> result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xy() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xz() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_yz() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xx() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_yy() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_zz() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xyz() const -> vec<Q> {
    vec<Q> inter_result = ::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3);
    vec<Q> inter_result2= ::partial_deriv<Q,5>(inter_result, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2);
    vec<Q> result       = ::partial_deriv<Q,5>(inter_result2, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1);
    return result;
}

template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_maxK3() const -> double {
    double max_K3 = (::power2(::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 3))
                   + ::power2(::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.f.ts, vertexContainerBase<Q,5>::dims, 2))
                   + ::power2(::partial_deriv<Q,5>(vertexContainerBase<Q,5>::data, frequencies_K3.b.ts, vertexContainerBase<Q,5>::dims, 1))
    ).max_norm();
    return max_K3;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_curvature_maxK3() const -> double {
    double max_K3 = (   ::power2(get_deriv_K3_xx())
                      + ::power2(get_deriv_K3_yy())
                      + ::power2(get_deriv_K3_zz())
                      + ::power2(get_deriv_K3_xy())
                      + ::power2(get_deriv_K3_xz())
                      + ::power2(get_deriv_K3_yz())
                    ).max_norm();
    return max_K3;
}

template <typename Q> double vertexDataContainer<k3,Q>::analyze_tails_K3_x() const {
    double maxabs_K3_total = vertexContainerBase<Q,5>::data.max_norm();
    vec<double> maxabsK3_along_w = maxabs(vertexContainerBase<Q,5>::data, vertexContainerBase<Q,5>::dims, 1);

    return maxabsK3_along_w[1] / maxabs_K3_total;
}
template <typename Q> double vertexDataContainer<k3,Q>::analyze_tails_K3_y() const {
    double maxabs_K3_total = vertexContainerBase<Q,5>::data.max_norm();
    vec<double> maxabsK3_along_v = maxabs(vertexContainerBase<Q,5>::data, vertexContainerBase<Q,5>::dims, 2);

    return maxabsK3_along_v[1] / maxabs_K3_total;
}
template <typename Q> double vertexDataContainer<k3,Q>::analyze_tails_K3_z() const {
    double maxabs_K3_total = vertexContainerBase<Q,5>::data.max_norm();
    vec<double> maxabsK3_along_vp = maxabs(vertexContainerBase<Q,5>::data, vertexContainerBase<Q,5>::dims, 3);

    return maxabsK3_along_vp[1] / maxabs_K3_total;
}

template <typename Q> auto vertexDataContainer<k3,Q>::shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k3> {

    VertexFrequencyGrid<k3> frequencies_new = frequencies_K3;

    vec<double> maxabsK3_along_w = maxabs(vertexContainerBase<Q,5>::data, vertexContainerBase<Q,5>::dims, 1);
    double maxmax = maxabsK3_along_w.max_norm();

    //const double rel_tail_threshold = 1e-4;
    size_t index = 0;
    while (true) {
        if (maxabsK3_along_w[index] > rel_tail_threshold * maxmax) break;
        index++;
    }
    frequencies_new.b.set_w_upper(std::abs(frequencies_K3.b.get_ws(index)));
    frequencies_new.b.initialize_grid();



    vec<double> maxabsK3_along_v = maxabs(vertexContainerBase<Q,5>::data, vertexContainerBase<Q,5>::dims, 2);
    vec<double> maxabsK3_along_vp= maxabs(vertexContainerBase<Q,5>::data, vertexContainerBase<Q,5>::dims, 3);
    vec<double> maxabsK3_along_f = elementwise([](const double& l, const double & r) -> double {return std::max(l, r);}, maxabsK3_along_v, maxabsK3_along_vp);

            //const double rel_tail_threshold = 1e-4;
    index = 0;
    while (true) {
        if (maxabsK3_along_f[index] > rel_tail_threshold * maxmax) break;
        index++;
    }
    frequencies_new.f.set_w_upper(std::abs(frequencies_K3.f.get_ws(index)));
    frequencies_new.f.initialize_grid();

    return frequencies_new;
}

#endif //FPP_MFRG_VERTEX_DATA_H
