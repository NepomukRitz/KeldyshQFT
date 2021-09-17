#ifndef FPP_MFRG_VERTEX_DATA_H
#define FPP_MFRG_VERTEX_DATA_H

/**
 * This header contains a class which is responsible for saving and retrieving vertex data
 */


#include "data_structures.h"          // real/complex vector classes
#include "parameters/master_parameters.h"               // system parameters (lengths of vectors etc.)
#include "grids/frequency_grid.h"            // functionality for the internal structure of the Hubbard model

template <typename Q> class fullvert; // forward declaration of fullvert
template <typename Q> class State; // forward declaration of fullvert

template <int k, typename Q>
class vertexDataContainer{
    friend void check_Kramers_Kronig(std::string filename);
    friend void test_PT4(double Lambda, bool write_flag);
    //friend void susceptibilities_postprocessing<>(Vertex<Q>& chi, Vertex<Q>& chi_diff,
    //                                            const State<Q>& state, const double Lambda);

};

template<typename Q>
class vertexDataContainer<k1, Q> {

public:
    VertexFrequencyGrid<k1> frequencies_K1;    // frequency grid

    explicit vertexDataContainer(double Lambda) : frequencies_K1(Lambda) { };

    vec<Q> K1 = vec<Q> (nK_K1 * nw1 * n_in);  // data points of K1

    /** K1-functionality */
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

    vec<Q> get_deriv_K1_x(bd_type left, bd_type right, Q value_left, Q value_right) const;
    double get_deriv_maxK1() const;
};

template<typename Q>
class vertexDataContainer<k2, Q> {

private:
    vec<Q> empty_K2() { // for pure K1-calculation no memory should be allocated unnecessarily for K2
        if (MAX_DIAG_CLASS >= 2) return vec<Q> (nK_K2 * nw2 * nv2 * n_in);  // data points of K2;
        else                     return vec<Q> (0);                         // empty vector, never used in calculations
    }

    /// Member functions for accessing the reducible vertex in channel r at arbitrary frequencies ///
    /// by interpolating stored data, in all possible channel-dependent frequency representations ///

public:
    size_t dimsK2[4] = {nK_K2, nBOS2, nFER2, n_in};
    VertexFrequencyGrid<k2> frequencies_K2;    // frequency grid

    explicit vertexDataContainer(double Lambda) : frequencies_K2(Lambda) { };

    /** K2 functionality */
    vec<Q> K2 = empty_K2();

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

    vec<Q> get_deriv_K2_x(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K2_y(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K2_xy(bd_type left, bd_type right, Q value_left, Q value_right) const;
    double get_deriv_maxK2() const;
};

template <typename Q>
class vertexDataContainer<k3, Q>{

private:
    vec<Q> empty_K3() { // for  K2-calculation no memory should be allocated unnecessarily for K3
        if (MAX_DIAG_CLASS >= 3) return vec<Q> (nK_K3 * nw3 * nv3 * nv3 * n_in);    // data points of K3  // data points of K2;
        else                     return vec<Q> (0);                                 // empty vector, never used in calculations
    }


public:
    size_t dimsK3[5] = {nK_K3, nBOS3, nFER3, nFER3, n_in};
    VertexFrequencyGrid<k3> frequencies_K3;    // frequency grid

    explicit vertexDataContainer(double Lambda) : frequencies_K3(Lambda) { };


    /** K3 functionality */
    vec<Q> K3 = empty_K3();



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


    vec<Q> get_deriv_K3_x(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K3_y(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K3_z(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K3_xy(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K3_xz(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K3_yz(bd_type left, bd_type right, Q value_left, Q value_right) const;
    vec<Q> get_deriv_K3_xyz(bd_type left, bd_type right, Q value_left, Q value_right) const;
    double get_deriv_maxK3() const;
};

/****************************************** MEMBER FUNCTIONS OF THE R-VERTEX ******************************************/
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
template <typename Q> auto vertexDataContainer<k1,Q>::K1_val(int iK, int iw, int i_in) const -> Q {
    return K1[iK*nw1*n_in + iw*n_in + i_in];
}


template <typename Q> auto vertexDataContainer<k1,Q>::get_deriv_maxK1() const -> double {
    double max_K1 = ::power2(::get_finite_differences(K1)).max_norm();
    return max_K1;

}



template <typename Q> auto vertexDataContainer<k1,Q>::get_deriv_K1_x(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_x[3] = {n_in, nK_K1, nBOS};
    const size_t perm_x[3] = {1, 2, 0};
    vec<Q> result = ::get_finite_differences<Q,3>(K1, frequencies_K1.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;

}



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
template <typename Q> auto vertexDataContainer<k2,Q>::K2_val(int iK, int iw, int iv, int i_in) const -> Q {
    return K2[iK * nw2 * nv2 * n_in + iw * nv2 * n_in + iv * n_in + i_in];
}




template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_x(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_x[4] = {nFER2, n_in, nK_K2, nBOS2};
    const size_t perm_x[4] = {2, 3, 0, 1};
    vec<Q> result = ::get_finite_differences<Q,4>(K2, frequencies_K2.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_y(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_y[4] = {n_in, nK_K2, nBOS2, nFER2};
    const size_t perm_y[4] = {1, 2, 3, 0};
    vec<Q> result = ::get_finite_differences<Q,4>(K2, frequencies_K2.f.ts, dims_y, perm_y, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_K2_xy(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_y[4] = {n_in, nK_K2, nBOS2, nFER2};
    const size_t perm_y[4] = {1, 2, 3, 0};
    vec<Q> inter_result = ::get_finite_differences<Q,4>(K2, frequencies_K2.f.ts, dims_y, perm_y, left, right, value_left, value_right);
    const size_t dims_x[4] = {nFER2, n_in, nK_K2, nBOS2};
    const size_t perm_x[4] = {2, 3, 0, 1};
    vec<Q> result = ::get_finite_differences<Q,4>(inter_result, frequencies_K2.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k2,Q>::get_deriv_maxK2() const -> double {
    size_t dims1[4] = {n_in, nK_K2, nBOS2, nFER2};
    size_t dims2[4] = {nFER2, n_in, nK_K2, nBOS2};
    size_t perm1[4] = {1, 2, 3, 0};
    size_t perm2[4] = {2, 3, 0, 1};
    double max_K2 = (::power2(::get_finite_differences<Q,4>(K2, frequencies_K2.f.ts, dims1, perm1))
                   + ::power2(::get_finite_differences<Q,4>(K2, frequencies_K2.b.ts, dims2, perm2))
    ).max_norm();
    return max_K2;

}

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
template <typename Q> auto vertexDataContainer<k3,Q>::K3_val(int iK, int iw, int iv, int ivp, int i_in) const -> Q {
    return K3[iK*nw3*nv3*nv3*n_in + iw*nv3*nv3*n_in + iv*nv3*n_in + ivp*n_in + i_in];
}



template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_x(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_x[5] = {nFER3, nFER3, n_in, nK_K3, nBOS3};
    const size_t perm_x[5] = {3, 4, 0, 1, 2};
    vec<Q> result = ::get_finite_differences<Q,5>(K3, frequencies_K3.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_y(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_y[5] = {nFER3, n_in, nK_K3, nBOS3, nFER3};
    const size_t perm_y[5] = {2, 3, 4, 0, 1};
    vec<Q> result = ::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims_y, perm_y, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_z(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_z[5] = {n_in, nK_K3, nBOS3, nFER3, nFER3};
    const size_t perm_z[5] = {1, 2, 3, 4, 0};
    vec<Q> result = ::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims_z, perm_z, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xy(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_y[5] = {nFER3, n_in, nK_K3, nBOS3, nFER3};
    const size_t perm_y[5] = {2, 3, 4, 0, 1};
    vec<Q> inter_result = ::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims_y, perm_y, left, right, value_left, value_right);
    const size_t dims_x[5] = {nFER3, nFER3, n_in, nK_K3, nBOS3};
    const size_t perm_x[5] = {3, 4, 0, 1, 2};
    vec<Q> result = ::get_finite_differences<Q,5>(inter_result, frequencies_K3.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xz(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_z[5] = {n_in, nK_K3, nBOS3, nFER3, nFER3};
    const size_t perm_z[5] = {1, 2, 3, 4, 0};
    vec<Q> inter_result = ::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims_z, perm_z, left, right, value_left, value_right);
    const size_t dims_x[5] = {nFER3, nFER3, n_in, nK_K3, nBOS3};
    const size_t perm_x[5] = {3, 4, 0, 1, 2};
    vec<Q> result = ::get_finite_differences<Q,5>(inter_result, frequencies_K3.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_yz(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_y[5] = {nFER3, n_in, nK_K3, nBOS3, nFER3};
    const size_t perm_y[5] = {2, 3, 4, 0, 1};
    vec<Q> inter_result = ::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims_y, perm_y, left, right, value_left, value_right);
    const size_t dims_z[5] = {n_in, nK_K3, nBOS3, nFER3, nFER3};
    const size_t perm_z[5] = {1, 2, 3, 4, 0};
    vec<Q> result = ::get_finite_differences<Q,5>(inter_result, frequencies_K3.f.ts, dims_z, perm_z, left, right, value_left, value_right);
    return result;
}
template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_K3_xyz(const bd_type left, const bd_type right, const Q value_left, const Q value_right) const -> vec<Q> {
    const size_t dims_z[5] = {n_in, nK_K3, nBOS3, nFER3, nFER3};
    const size_t perm_z[5] = {1, 2, 3, 4, 0};
    vec<Q> inter_result = ::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims_z, perm_z, left, right, value_left, value_right);
    const size_t dims_y[5] = {nFER3, n_in, nK_K3, nBOS3, nFER3};
    const size_t perm_y[5] = {2, 3, 4, 0, 1};
    vec<Q> inter_result2= ::get_finite_differences<Q,5>(inter_result, frequencies_K3.f.ts, dims_y, perm_y, left, right, value_left, value_right);
    const size_t dims_x[5] = {nFER3, nFER3, n_in, nK_K3, nBOS3};
    const size_t perm_x[5] = {3, 4, 0, 1, 2};
    vec<Q> result = ::get_finite_differences<Q,5>(inter_result2, frequencies_K3.b.ts, dims_x, perm_x, left, right, value_left, value_right);
    return result;
}

template <typename Q> auto vertexDataContainer<k3,Q>::get_deriv_maxK3() const -> double {
    size_t dims1[5] = {n_in, nK_K3, nBOS3, nFER3, nFER3};
    size_t dims2[5] = {nFER3, n_in, nK_K3, nBOS3, nFER3};
    size_t dims3[5] = {nFER3, nFER3, n_in, nK_K3, nBOS3};
    size_t perm1[5] = {1, 2, 3, 4, 0};
    size_t perm2[5] = {2, 3, 4, 0, 1};
    size_t perm3[5] = {3, 4, 0, 1, 2};

    double max_K3 = (::power2(::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims1, perm1))
                   + ::power2(::get_finite_differences<Q,5>(K3, frequencies_K3.f.ts, dims2, perm2))
                   + ::power2(::get_finite_differences<Q,5>(K3, frequencies_K3.b.ts, dims3, perm3))
    ).max_norm();
    return max_K3;

}

#endif //FPP_MFRG_VERTEX_DATA_H
