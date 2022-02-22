#ifndef FPP_MFRG_VERTEX_DATA_H
#define FPP_MFRG_VERTEX_DATA_H

/**
 * This header contains a class which is responsible for saving and retrieving vertex data
 */

#include "../../utilities/template_utils.hpp"
#include "../../utilities/math_utils.hpp"
#include "../../symmetries/Keldysh_symmetries.hpp"
#include "../../data_structures.hpp"          // real/complex vector classes
#include "../../parameters/master_parameters.hpp"               // system parameters (lengths of vectors etc.)
#include "../../parameters/frequency_parameters.hpp"
#include "../../grids/frequency_grid.hpp"            // functionality for the internal structure of the Hubbard model
#include "H5Cpp.h"

/// TODO: treat VertexFrequencyGrid as member of vertexContainerBase via template

template <typename Q> class rvert; // forward declaration of rvert
template <typename Q> class fullvert; // forward declaration of fullvert
template <typename Q> class State; // forward declaration of State
template <typename Q, vertexType symm_type> class GeneralVertex;
//template <typename Q>class symmetric_full;
template <typename Q>using Vertex = GeneralVertex<Q, symmetric_full>;
class Buffer;

/**
 * Offers basic functionality that is identical for all K_classes
 * @tparam Q        data type of vertex data
 * @tparam rank     rank of data tensor
 */
template <typename Q, std::size_t rank>
class vertexContainerBase {
    friend class State<Q>;
    template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator+= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator- (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator-= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator/ (rvert<T> lhs, const rvert<T>& rhs);
    template<typename T> friend rvert<T> rvert<T>::operator/= (const rvert<T>& rhs);
    template<typename T> friend rvert<T> operator* (rvert<T> lhs, const double& alpha);
    template<typename T> friend rvert<T> rvert<T>::operator*= (double alpha);
    template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const double& alpha);
    template<typename T> friend rvert<T> rvert<T>::operator+= (double alpha);
    friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);


protected:
    using buffer_type = multidimensional:: multiarray<Q,rank>;
    using index_type = typename buffer_type::index_type;
    using dimensions_type = typename buffer_type::dimensions_type;

    buffer_type data;

public:
    /// constructor:
    vertexContainerBase() = default;
    explicit vertexContainerBase(const dimensions_type dims) : data(dims) {};
    template <typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
    explicit vertexContainerBase(const Types &... dims) : vertexContainerBase(index_type({static_cast<size_t>(dims)...})) {};
    explicit vertexContainerBase(const buffer_type &data_in) : buffer_type(data_in) {};

    /// Access vertex data via flattened index.
    Q acc(const size_t flatIndex) const {assert(flatIndex<data.size()); return data.flat_at(flatIndex);}
    void direct_set(const size_t flatIndex, Q value) {assert(flatIndex < data.size()); data.flat_at(flatIndex) = value;}

    /// Returns value for a multiIndex
    template <typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
    Q val(const Types &... i) const {return data(i...);}
    Q val(const index_type& idx) const {return data.at(idx);}
    /// Returns reference to a value for a multiIndex
    template <typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
    const Q& at(const Types &... i) const {return data.at(i...);}
    template <std::size_t freqrank, std::size_t vecsize, typename... Types,
            typename std::enable_if_t<(sizeof...(Types) == freqrank+pos_first_freq+1) and (are_all_integral<size_t, Types...>::value), bool> = true>
                    auto val_vectorized(const Types &... i) const -> Eigen::Matrix<Q,vecsize,1>{
        return data.template at_vectorized<pos_first_freq, freqrank, vecsize>(i...);
    }
    /// Sets a value at a multiIndex
    template <typename... Types
            ,typename std::enable_if_t<(sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true
    >
    void setvert(const Q value, const Types &... i) {data.at(i...) = value;}
    void setvert(const Q value, const index_type & idx) {data.at(idx) = value;}

    auto get_dims() const {return data.length();}

    /// Returns the buffer "data" containing the vertex data
    const buffer_type& get_vec() const {return data;}
    /// Sets the the buffer "data"
    template <typename container,
            std::enable_if_t<std::is_same_v<typename container::value_type, Q> && !std::is_same_v<container, buffer_type>, bool> = true
    >
    void set_vec(const container &data_in) {assert(data.size() == data_in.size()); data = buffer_type(data.length(), data_in);}
    void set_vec(const buffer_type &data_in) {assert(data.is_same_length(data_in)); data = data_in;}
    void set_vec(const buffer_type &&data_in) {assert(data.is_same_length(data_in)); data = data_in;}
    /// Adds a vector to the data
    template <typename container,
            std::enable_if_t< std::is_same_v<typename container::value_type, Q> && !std::is_same_v<container, buffer_type>, bool> = true
    >
    void add_vec(const container &summand) {
        assert(data.size() == summand.size()); /// Check that summand has the right length
        data += buffer_type(data.length(), summand);
    }
    void add_vec(const buffer_type &summand) {
        assert(data.is_same_length(summand)); /// Check that summand has the right length
        data += summand;
    }

    constexpr auto eigen_segment(const index_type &start, const index_type &end)
    {
        return data.eigen_segment(start, end);
    }

    constexpr auto eigen_segment(const index_type &start, const index_type &end) const
    {
        return data.eigen_segment(start, end);
    }

};


template<K_class k, typename Q>
class vertexDataContainer: public vertexContainerBase<Q,5> {
    friend void test_PT4(double Lambda, bool write_flag);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template<typename T> friend void check_FDTs(const State<T>& state, bool verbose);
    friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);

protected:
    using base_class = vertexContainerBase<Q,5>;
public:
    using index_type = typename base_class::index_type;
    using dimensions_type = typename base_class::dimensions_type;
    using buffer_type = typename base_class::buffer_type;


    VertexFrequencyGrid<k2> frequencies;    // frequency grid
    vertexDataContainer() = default;
    explicit vertexDataContainer(double Lambda, dimensions_type dims) : frequencies(Lambda), base_class(dims) { };

    /// Functions for getting and setting the frequency grid and its members

    /**
     * gets the frequency corresponding to the frequency index i
     */
    //void get_freqs_w(double& w, double& v, int iw, int iv) const;
    //void K2_get_freqs_aux(double& w, double& v, int iw, int iv) const;

    auto get_VertexFreqGrid() const -> const VertexFrequencyGrid<k2>&;
    void set_VertexFreqGrid(const VertexFrequencyGrid<k2> frequencyGrid);

    //const double& K2_get_wlower_b() const;
    //const double& K2_get_wupper_b() const;
    //const double& frequencies.get_wlower_f() const;
    //const double& K2_get_wupper_f() const;
    //const FrequencyGrid& K2_get_freqGrid_b() const;
    //const FrequencyGrid& K2_get_freqGrid_f() const;
    //const double& K2_get_tlower_b_aux() const;
    //const double& K2_get_tupper_b_aux() const;
    //const double& K2_get_tlower_f_aux() const;
    //const double& K2_get_tupper_f_aux() const;
    //auto K2_gridtransf_b(double w) const -> double;
    //auto K2_gridtransf_f(double w) const -> double;
    //auto K2_gridtransf_inv_b(double w) const -> double;
    //auto K2_gridtransf_inv_f(double w) const -> double;


    auto K2_get_correction_MFfiniteT(int iw) const -> double;


    buffer_type get_deriv_K2_x (int order=5) const;
    buffer_type get_deriv_K2_y (int order=5) const;
    buffer_type get_deriv_K2_xy(int order=5) const;
    buffer_type get_deriv_K2_xx(int order=5) const;
    buffer_type get_deriv_K2_yy(int order=5) const;
    double get_deriv_maxK2() const;
    auto get_curvature_maxK2() const -> double;

    double analyze_tails_K2_x() const;

    double analyze_tails_K2_y() const;
    auto shrink_freq_box(const double rel_tail_threshold, bool verbose=true) const -> VertexFrequencyGrid<k2>;

};


/**
 * vertex data container for K1
 */
template<typename Q>
class vertexDataContainer<k1, Q> : public vertexContainerBase<Q,4>{
    friend void check_Kramers_Kronig(std::string filename);
    friend void test_PT4(double Lambda, bool write_flag);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template<typename T> friend void susceptibilities_postprocessing(Vertex<T>& chi, Vertex<T>& chi_diff, const State<T>& state, double Lambda);
    friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);


protected:
    using base_class = vertexContainerBase<Q,4>;
public:
    using index_type = typename base_class::index_type;
    using dimensions_type = typename base_class::dimensions_type;
    using buffer_type = typename base_class::buffer_type;

    VertexFrequencyGrid<k1> frequencies;    // frequency grid

    vertexDataContainer() = default;
    explicit vertexDataContainer(double Lambda, dimensions_type dims) : frequencies(Lambda), base_class(dims) { };

    /// Functions for getting and setting frequency grid and its members;
    /// TODO: Can probably be shifted to vertexContainerBase    (problems: non-existing fermionic grid for K1 ==> use std::enable_if; need multiple parameter packs for get_freqs()
    /// Or shift these to the VertexFrequencyGrid class?
    /// For now I keep this to retain flexibility, e.g. to store data on very different grids in K1/K2/K3

    /**
     * gets the frequency corresponding to the frequency index i
     */
    //void frequencies.get_freqs_w(double& w, int i) const;     /// returns regular frequency
    //void K1_get_freq_aux(double& w, int i) const;   /// returns frequency on the auxiliary grid
//
    auto get_VertexFreqGrid() const -> const VertexFrequencyGrid<k1>&;
    void set_VertexFreqGrid(const VertexFrequencyGrid<k1> frequencyGrid);
//
    //auto frequencies.get_freqGrid_b() const -> const FrequencyGrid&;
//
    //const double& frequencies.get_wlower() const;
    //const double& frequencies.get_wupper_b() const;
    //const double& K1_get_tlower_aux() const;
    //const double& K1_get_tupper_aux() const;
    //auto K1_gridtransf(double w) const -> double;
    //auto K1_gridtransf_inv(double w) const -> double;

    /// Compute derivative in w-direction using a finite-differences method
    buffer_type get_deriv_K1_x(int order=5) const;
    /// Compute the maximum norm of the gradient
    double get_deriv_maxK1() const;
    /// Compute the maximum norm of the curvature
    double get_curvature_maxK1() const;

    double analyze_tails_K1() const;

    /// shrink the frequency box if the data on the outermost gridpoints is smaller than data.maxnorm()*rel_tail_threshold
    auto shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k1>;
};


template <typename Q>
class vertexDataContainer<k3, Q>: public vertexContainerBase<Q,6> {
    friend void test_PT4(double Lambda, bool write_flag);
    template <typename T> friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);
    template <typename T> friend void result_set_frequency_grids(State<T>& result, Buffer& buffer);
    template<typename T> friend void check_FDTs(const State<T>& state, bool verbose);
    friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);


protected:
    using base_class = vertexContainerBase<Q,6>;


public:
    using index_type = typename base_class::index_type;
    using dimensions_type = typename base_class::dimensions_type;
    using buffer_type = typename base_class::buffer_type;

    VertexFrequencyGrid<k3> frequencies;    // frequency grid

    vertexDataContainer() = default;
    explicit vertexDataContainer(double Lambda, dimensions_type dims) : frequencies(Lambda), base_class(dims) { };


    //void K3.frequencies.get_freqs_w(double& w, double& v, double& vp, int iw, int iv, int ivp, char channel) const;
    //void K3_get_freqs_aux(double& w, double& v, double& vp, int iw, int iv, int ivp) const;

    auto get_VertexFreqGrid() const -> const VertexFrequencyGrid<k3>&;
    void set_VertexFreqGrid(const VertexFrequencyGrid<k3> frequencyGrid);
    //const FrequencyGrid& frequencies.get_freqGrid_b() const;
    //const FrequencyGrid& frequencies.get_freqGrid_f() const;
//
    //const double& K3_get_wlower_b() const;
    //const double& K3_get_wupper_b() const;
    //const double& K3_get_wlower_f() const;
    //const double& K3_get_wupper_f() const;
    //const double& K3_get_tlower_b_aux() const;
    //const double& K3_get_tupper_b_aux() const;
    //const double& K3_get_tlower_f_aux() const;
    //const double& K3_get_tupper_f_aux() const;
    //auto K3_gridtransf_b(double w) const -> double;
    //auto K3_gridtransf_f(double w) const -> double;
    //auto K3_gridtransf_inv_b(double w) const -> double;
    //auto K3_gridtransf_inv_f(double w) const -> double;

    auto K3_get_correction_MFfiniteT(int iw) const -> double;


    buffer_type get_deriv_K3_x  (int order=5) const;
    buffer_type get_deriv_K3_y  (int order=5) const;
    buffer_type get_deriv_K3_z  (int order=5) const;
    buffer_type get_deriv_K3_xy (int order=5) const;
    buffer_type get_deriv_K3_xz (int order=5) const;
    buffer_type get_deriv_K3_yz (int order=5) const;
    buffer_type get_deriv_K3_xx (int order=5) const;
    buffer_type get_deriv_K3_yy (int order=5) const;
    buffer_type get_deriv_K3_zz (int order=5) const;
    buffer_type get_deriv_K3_xyz(int order=5) const;
    double get_deriv_maxK3() const;
    auto get_curvature_maxK3() const -> double;
    double analyze_tails_K3_x() const;
    double analyze_tails_K3_y() const;
    double analyze_tails_K3_z() const;

    auto shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k3>;

};

/************************************ MEMBER FUNCTIONS OF THE VERTEX Data Container************************************/

/// K1:

template<typename Q>
auto vertexDataContainer<k1,Q>::get_VertexFreqGrid() const -> const VertexFrequencyGrid<k1>& {
    return frequencies;
}
template<typename Q>
void vertexDataContainer<k1,Q>::set_VertexFreqGrid(const VertexFrequencyGrid<k1> frequencyGrid) {
    frequencies = frequencyGrid;
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::get_VertexFreqGrid() const -> const VertexFrequencyGrid<k2>& {
    return frequencies;
}
template<K_class k, typename Q>
void vertexDataContainer<k,Q>::set_VertexFreqGrid(const VertexFrequencyGrid<k2> frequencyGrid) {
    frequencies = frequencyGrid;
}
template<typename Q>
auto vertexDataContainer<k3,Q>::get_VertexFreqGrid() const -> const VertexFrequencyGrid<k3> & {
    return frequencies;
}
template<typename Q>
void vertexDataContainer<k3,Q>::set_VertexFreqGrid(const VertexFrequencyGrid<k3> frequencyGrid) {
    frequencies = frequencyGrid;
}

/*
template<typename Q>
const double& vertexDataContainer<k1,Q>::frequencies.get_wlower() const {
    return frequencies.b.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k1,Q>::frequencies.get_wupper_b() const {
    return frequencies.b.w_upper;
}
template<typename Q>
auto vertexDataContainer<k1,Q>::frequencies.get_freqGrid_b() const -> const FrequencyGrid& {
    return frequencies.b;
}
template<typename Q>
void vertexDataContainer<k1,Q>::frequencies.get_freqs_w(double& w, const int i) const {
    w = frequencies.b.get_ws(i);
}
template<typename Q>
const double& vertexDataContainer<k1,Q>::K1_get_tlower_aux() const {
    return frequencies.b.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k1,Q>::K1_get_tupper_aux() const {
    return frequencies.b.t_upper;
}
template<typename Q>
void vertexDataContainer<k1,Q>::K1_get_freq_aux(double& w, const int i) const {
    w = frequencies.b.get_ts(i);
}
template<typename Q>
auto vertexDataContainer<k1,Q>::K1_gridtransf(double w) const -> double {
    return frequencies.b.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k1,Q>::K1_gridtransf_inv(double w) const -> double {
    return frequencies.b.grid_transf_inv(w);
}
*/

template <typename Q> auto vertexDataContainer<k1,Q>::get_deriv_K1_x(const int order) const -> buffer_type {
    buffer_type result = ::partial_deriv<Q,4>(base_class::data, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k1,Q>::get_deriv_maxK1() const -> double {
    int order = 3;
    double dt = frequencies.get_freqGrid_b().dt;
    double Kmax = base_class::get_vec().max_norm();
    double max_K1 = ::power2(get_deriv_K1_x(order)*dt*(1/Kmax)).max_norm();  // normalize by magnitude of vertex contribution and grid spacing
    return max_K1;
}
template <typename Q> auto vertexDataContainer<k1,Q>::get_curvature_maxK1() const -> double {
    int order = 3;
    double dt = frequencies.get_freqGrid_b().dt;
    double Kmax = base_class::get_vec().max_norm();
    double max_K1 = ::power2(::partial_deriv<Q,4>(::partial_deriv<Q,4>(base_class::data, frequencies.b.ts, base_class::data.length(), 2, order), frequencies.b.ts, base_class::data.length(), 2, order) * dt * dt * (1 / Kmax)).max_norm();
    return max_K1;
}


template <typename Q> double vertexDataContainer<k1,Q>::analyze_tails_K1() const {
    double maxabs_K1_total = base_class::data.max_norm();
    vec<double> maxabsK1_along_w = maxabs(base_class::data, base_class::data.length(), 2);

    return maxabsK1_along_w[0] / maxabs_K1_total;
}


template <typename Q> auto vertexDataContainer<k1,Q>::shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k1>
{
    double maxmax = base_class::data.max_norm();
    vec<double> maxabsK1_along_w = maxabs(base_class::data, base_class::data.length(), 2) * (1/maxmax);
    VertexFrequencyGrid<k1> frequencies_new = frequencies;
    if (std::abs(maxmax) > 1e-30) { frequencies_new.b = freqGrid::shrink_freq_box(frequencies.b, rel_tail_threshold, maxabsK1_along_w); }

    return frequencies_new;
}


/// K2:
/*
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_wlower_b() const {
    return frequencies.b.w_lower;
}
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_wupper_b() const {
    return frequencies.b.w_upper;
}
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::frequencies.get_wlower_f() const {
    return frequencies.f.w_lower;
}
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_wupper_f() const {
    return frequencies.f.w_upper;
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_get_freqGrid_b() const -> const FrequencyGrid& {
    return frequencies.b;
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_get_freqGrid_f() const -> const FrequencyGrid& {
    return frequencies.f;
}
template<K_class k, typename Q>
void vertexDataContainer<k,Q>::frequencies.get_freqs_w(double &w, double &v, const int iw, const int iv) const {
    frequencies.get_freqs_w(w, v, iw, iv);
}

template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_tlower_b_aux() const {
    return frequencies.b.t_lower;
}
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_tupper_b_aux() const {
    return frequencies.b.t_upper;
}
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_tlower_f_aux() const {
    return frequencies.f.t_lower;
}
template<K_class k, typename Q>
const double& vertexDataContainer<k,Q>::K2_get_tupper_f_aux() const {
    return frequencies.f.t_upper;
}
template<K_class k, typename Q>
void vertexDataContainer<k,Q>::K2_get_freqs_aux(double &w, double &v, const int iw, const int iv) const {
    frequencies.get_freqs_aux(w, v, iw, iv);
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_gridtransf_b(double w) const -> double {
    return frequencies.b.grid_transf(w);
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_gridtransf_f(double w) const -> double {
    return frequencies.f.grid_transf(w);
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_gridtransf_inv_b(double w) const -> double {
    return frequencies.b.grid_transf_inv(w);
}
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_gridtransf_inv_f(double w) const -> double {
    return frequencies.f.grid_transf_inv(w);
}
 */
template<K_class k, typename Q>
auto vertexDataContainer<k,Q>::K2_get_correction_MFfiniteT(int iw) const -> double {
    if (not KELDYSH and not ZERO_T)
        return signFlipCorrection_MF(frequencies.b.get_ws(iw));
    else return 0.;
}

template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_deriv_K2_x(const int order) const -> buffer_type {
    buffer_type result = ::partial_deriv<Q,5>(base_class::data, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_deriv_K2_y(const int order) const -> buffer_type {
    buffer_type result = ::partial_deriv<Q,5>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    return result;
}
template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_deriv_K2_xy(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,5>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    buffer_type result       = ::partial_deriv<Q,5>(inter_result, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_deriv_K2_xx(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,5>(base_class::data, frequencies.b.ts, base_class::data.length(), 2, order);
    buffer_type result       = ::partial_deriv<Q,5>(inter_result, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_deriv_K2_yy(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,5>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    buffer_type result       = ::partial_deriv<Q,5>(inter_result, frequencies.f.ts, base_class::data.length(), 3, order);
    return result;
}
template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_deriv_maxK2() const -> double {
    const int order = 3;

    double Kmax = base_class::get_vec().max_norm();
    double dtb = frequencies.get_freqGrid_b().dt;
    double dtf = frequencies.get_freqGrid_f().dt;
    double max_K2 = (     ::power2(get_deriv_K2_x(order)*dtb*(1/Kmax)) // multiply by dt to make the result independent of the grid spacing
                        + ::power2(get_deriv_K2_y(order)*dtf*(1/Kmax))
                    ).max_norm();
    return max_K2;
}
template <K_class k, typename Q> auto vertexDataContainer<k,Q>::get_curvature_maxK2() const -> double {
    const int order = 3;

    double Kmax = base_class::get_vec().max_norm();
    double dtb = frequencies.get_freqGrid_b().dt;
    double dtf = frequencies.get_freqGrid_f().dt;
    double max_K1 = (    ::power2( get_deriv_K2_xx(order)*dtb*dtb*(1/Kmax))
                        +::power2( get_deriv_K2_yy(order)*dtf*dtf*(1/Kmax))
                        +::power2( get_deriv_K2_xy(order)*dtb*dtf*(1/Kmax))*2.
                    ).max_norm();
    return max_K1;
}

template <K_class k, typename Q> double vertexDataContainer<k,Q>::analyze_tails_K2_x() const {
    double maxmax = base_class::data.max_norm();
    vec<double> maxabsK2_along_w = maxabs(base_class::data, base_class::data.length(), 2) * (1/maxmax);

    return maxabsK2_along_w[0];
}
template <K_class k, typename Q> double vertexDataContainer<k,Q>::analyze_tails_K2_y() const {
    double maxmax = base_class::data.max_norm();
    vec<double> maxabsK2_along_v = maxabs(base_class::data, base_class::data.length(), 3) * (1/maxmax);

    return maxabsK2_along_v[0];
}

template <K_class k, typename Q> auto vertexDataContainer<k,Q>::shrink_freq_box(const double rel_tail_threshold, const bool verbose) const -> VertexFrequencyGrid<k2> {

    VertexFrequencyGrid<k2> frequencies_new = frequencies;

    double maxmax = base_class::data.max_norm();
    vec<double> maxabsK2_along_w = maxabs(base_class::data, base_class::data.length(), 2) * (1/maxmax);
    vec<double> maxabsK2_along_v = maxabs(base_class::data, base_class::data.length(), 3) * (1/maxmax);


    frequencies_new.b = freqGrid::shrink_freq_box(frequencies.b, rel_tail_threshold, maxabsK2_along_w, verbose);
    if (verbose and mpi_world_rank() == 0) std::cout << "in direction w to " << frequencies_new.b.w_upper << std::endl;
    frequencies_new.f = freqGrid::shrink_freq_box(frequencies.f, rel_tail_threshold, maxabsK2_along_v, verbose);
    if (verbose and mpi_world_rank() == 0) std::cout << "in direction v to " << frequencies_new.f.w_upper << std::endl;

    return frequencies_new;
}


/// K3:
/*
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wlower_b() const {
    return frequencies.b.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wupper_b() const {
    return frequencies.b.w_upper;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wlower_f() const {
    return frequencies.f.w_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_wupper_f() const {
    return frequencies.f.w_upper;
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_freqGrid_b() const -> const FrequencyGrid& {
    return frequencies.b;
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_freqGrid_f() const -> const FrequencyGrid& {
    return frequencies.f;
}
template<typename Q>
void vertexDataContainer<k3,Q>::K3.frequencies.get_freqs_w(double &w, double &v, double& vp, const int iw, const int iv, const int ivp, const char channel) const {
    frequencies.get_freqs_w(w, v, vp, iw, iv, ivp, channel);
}

template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tlower_b_aux() const {
    return frequencies.b.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tupper_b_aux() const {
    return frequencies.b.t_upper;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tlower_f_aux() const {
    return frequencies.f.t_lower;
}
template<typename Q>
const double& vertexDataContainer<k3,Q>::K3_get_tupper_f_aux() const {
    return frequencies.f.t_upper;
}
template<typename Q>
void vertexDataContainer<k3,Q>::K3_get_freqs_aux(double &w, double &v, double& vp, const int iw, const int iv, const int ivp) const {
    frequencies.get_freqs_aux(w, v, vp, iw, iv, ivp);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_b(double w) const -> double {
    return frequencies.b.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_f(double w) const -> double {
    return frequencies.f.grid_transf(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_inv_b(double w) const -> double {
    return frequencies.b.grid_transf_inv(w);
}
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_gridtransf_inv_f(double w) const -> double {
    return frequencies.f.grid_transf_inv(w);
}
 */
template<typename Q>
auto vertexDataContainer<k3,Q>::K3_get_correction_MFfiniteT(int iw) const -> double {
    if (not KELDYSH and not ZERO_T)
    return floor2bfreq(frequencies.b.get_ws(iw) / 2) - ceil2bfreq(frequencies.b.get_ws(iw) / 2);
    else return 0.;
}


template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_x(const int order) const -> buffer_type {
    buffer_type result = ::partial_deriv<Q,6>(base_class::data, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_y(const int order) const -> buffer_type {
    buffer_type result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_z(const int order) const -> buffer_type {
    buffer_type result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 4, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xy(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xz(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 4, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_yz(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result, frequencies.f.ts, base_class::data.length(), 4, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xx(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.b.ts, base_class::data.length(), 2, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_yy(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 3, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result, frequencies.f.ts, base_class::data.length(), 3, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_zz(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 4, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result, frequencies.f.ts, base_class::data.length(), 4, order);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xyz(const int order) const -> buffer_type {
    buffer_type inter_result = ::partial_deriv<Q,6>(base_class::data, frequencies.f.ts, base_class::data.length(), 4, order);
    buffer_type inter_result2= ::partial_deriv<Q,6>(inter_result, frequencies.f.ts, base_class::data.length(), 3, order);
    buffer_type result       = ::partial_deriv<Q,6>(inter_result2, frequencies.b.ts, base_class::data.length(), 2, order);
    return result;
}

template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_maxK3() const -> double {
    const int order = 3;

    double Kmax = base_class::get_vec().max_norm();
    double dtb = frequencies.get_freqGrid_b().dt;
    double dtf = frequencies.get_freqGrid_f().dt;
    double max_K3 = (::power2(get_deriv_K3_z(order) * dtf * (1/Kmax))
                   + ::power2(get_deriv_K3_y(order) * dtf * (1/Kmax))
                   + ::power2(get_deriv_K3_x(order) * dtb * (1/Kmax))
    ).max_norm();
    return max_K3;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_curvature_maxK3() const -> double {
    const int order = 3;

    double Kmax = base_class::get_vec().max_norm();
    double dtb = frequencies.get_freqGrid_b().dt;
    double dtf = frequencies.get_freqGrid_f().dt;
    double max_K3 = (   ::power2(get_deriv_K3_xx(order)*dtb*dtb*(1/Kmax))
                      + ::power2(get_deriv_K3_yy(order)*dtf*dtf*(1/Kmax))
                      + ::power2(get_deriv_K3_zz(order)*dtf*dtf*(1/Kmax))
                      +(::power2(get_deriv_K3_xy(order)*dtb*dtf*(1/Kmax))
                      + ::power2(get_deriv_K3_xz(order)*dtb*dtf*(1/Kmax))
                      + ::power2(get_deriv_K3_yz(order)*dtf*dtf*(1/Kmax)))*2.
                    ).max_norm();
    return max_K3;
}

template <typename Q> double vertexDataContainer<k3,Q>::analyze_tails_K3_x() const {
    double maxabs_K3_total = base_class::data.max_norm();
    vec<double> maxabsK3_along_w = maxabs(base_class::data, base_class::data.length(), 2);

    return maxabsK3_along_w[1] / maxabs_K3_total;
}
template <typename Q> double vertexDataContainer<k3,Q>::analyze_tails_K3_y() const {
    double maxabs_K3_total = base_class::data.max_norm();
    vec<double> maxabsK3_along_v = maxabs(base_class::data, base_class::data.length(), 3);

    return maxabsK3_along_v[1] / maxabs_K3_total;
}
template <typename Q> double vertexDataContainer<k3,Q>::analyze_tails_K3_z() const {
    double maxabs_K3_total = base_class::data.max_norm();
    vec<double> maxabsK3_along_vp = maxabs(base_class::data, base_class::data.length(), 4);

    return maxabsK3_along_vp[1] / maxabs_K3_total;
}

template <typename Q> auto vertexDataContainer<k3,Q>::shrink_freq_box(const double rel_tail_threshold) const -> VertexFrequencyGrid<k3> {

    VertexFrequencyGrid<k3> frequencies_new = frequencies;

    double maxmax = base_class::data.max_norm();
    vec<double> maxabsK3_along_w = maxabs(base_class::data, base_class::data.length(), 2) * (1/maxmax);
    vec<double> maxabsK3_along_v = maxabs(base_class::data, base_class::data.length(), 3) * (1/maxmax);

    frequencies_new.b = freqGrid::shrink_freq_box(frequencies.b, rel_tail_threshold, maxabsK3_along_w);
    frequencies_new.f = freqGrid::shrink_freq_box(frequencies.f, rel_tail_threshold, maxabsK3_along_v);

    return frequencies_new;
}

#endif //FPP_MFRG_VERTEX_DATA_H
