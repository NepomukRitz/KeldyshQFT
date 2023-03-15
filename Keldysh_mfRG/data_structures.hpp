/**
 * Define essential data types:
 * comp        : complex number (complex<double>)
 * glb_i       : imaginary unit
 * vec         : vector class with additional functionality such as element-wise operations, real/imag part etc.
 * VertexInput : auxiliary struct that contains all input variables of vertices
 */

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <complex>          // for usage of complex numbers
#include <cmath>            // for math. operations (real, imag, abs etc.)
#include <vector>           // vec class is derived from vector class
#include <initializer_list> // to initialize vec class with initializer list
#include <functional>
#include <cassert>
#include <numeric>
#include "parameters/master_parameters.hpp"
#include "Eigen/Dense"
//#include "utilities/math_utils.h"

template<typename Q, bool differentiated> class State;

typedef std::complex<double> comp; // Complex number

constexpr comp glb_i (0., 1.);    // Imaginary unit

template <typename T> double myabs(const T& x) {
    if constexpr(std::is_same_v<T, comp> or std::is_same_v<T, double>) {
        return std::abs(x);
    }
    else if constexpr(std::is_same_v<T, State<comp,false>> or std::is_same_v<T, State<double,false>>) {
        return x.abs();
    }
    else {
        return x.template lpNorm<Eigen::Infinity>();
    }
}

template <typename T> T myzero() {
    if constexpr(std::is_same_v<T, comp> or std::is_same_v<T, double> or std::is_same_v<T, int>) {
        return (T)0.;
    }
    else {
        return T::Zero();
    }
}

template <typename T>
constexpr T get_Id_Matrix_4x4(){
    T res;
    res << 1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1, 0,
           0, 0, 0, 1;
    return res;
}

template <typename T> T myIdentity() {
    if constexpr(std::is_same_v<T, comp> or std::is_same_v<T, double> or std::is_same_v<T, int>) {
        return (T)1;
    }
    else {
        const T id = get_Id_Matrix_4x4<T>();
        return id;
    }
}

template <typename T> constexpr int myRowsAtCompileTime() {
    if constexpr(std::is_same_v<T,double> || std::is_same_v<T,comp>) {
        return 1;
    }
    else {
        return T::RowsAtCompileTime;
    }
}

template <typename T> constexpr int myColsAtCompileTime() {
    if constexpr(std::is_same_v<T,double> || std::is_same_v<T,comp>) {
        return 1;
    }
    else {
        return T::ColsAtCompileTime;
    }
}

#if PARTICLE_HOLE_SYMM and not KELDYSH_FORMALISM and not defined(HUBBARD)
using state_datatype = double;
#else
using state_datatype = comp;
#endif
template<typename Q>
inline auto myreal(const Q x) -> double {return x.real();};
template<>
inline auto myreal<double>(const double x) -> double {return x;};
template<typename Q>
inline auto myimag(const Q x) -> double {return x.imag();};
template<>
inline auto myimag<double>(const double x) -> double {return x;};
template<typename Q>
inline auto myconj(const Q x) -> Q {return conj(x);};
template<>
inline auto myconj<double>(const double x) -> double {return x;};

auto isfinite(comp z) -> bool;

template<typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr >
auto isfinite(T z) -> bool {
    return std::isfinite(z);
};


/// DECLARATIONS ///

/// General vector class, defining element-wise addition, subtraction and multiplication, as well as real/imag part etc.
template <typename T>
class vec : public std::vector<T> {
public:
    vec() : std::vector<T> () {}; 						 // trivial constructor
    explicit vec(int n) : std::vector<T> (n) {};				     // constructor with number of elements
    vec(int n, T value) : std::vector<T> (n, value) {};   // constructor with number of elements and value
    template <class InputIterator>
    vec (InputIterator first, InputIterator last)
     : std::vector<T> (first, last) {};                   // constructor from iterators to copy parts of existing vector
    vec(std::initializer_list<T> m) : std::vector<T> (m) {};   // constructor from initializer list
    vec(std::vector<T> vect) : std::vector<T> (vect) {};

    T operator() (int i) const {return (*this)[i]; }	     // operator for element access
    vec<T> operator() (int i1, int i2) const;              // get a subvector {x[i1], ..., x[i2]}

    vec<T> inv() const;   // element-wise inverse
    vec<double> real() const;   // element-wise real part
    vec<double> imag() const;   // element-wise imaginary part
    template<typename T2=double> vec<T2> abs() const;    // element-wise absolute value
    vec<T> conj() const;        // element-wise complex conjugate
    double max_norm() const;    // maximum norm
    vec<T> diff() const;        // vector of differences between adjacent elements
    T sum() const;              // sum of all elements

    template<typename Q>
    vec<T> & operator+= (const vec<Q>& m);     // element-wise addition of two vectors
    template<typename Q>
    vec<T> & operator+= (const Q& c);          // addition of a constant
    template<typename Q>
    vec<T> & operator-= (const vec<Q>& m);     // element-wise subtraction of two vectors
    template<typename Q>
    vec<T> & operator-= (const Q& c);          // subtraction of a constant
    template<typename Q>
    vec<T> & operator*= (const vec<Q>& m);     // element-wise multiplication of two vectors
    template<typename Q>
    vec<T> & operator*= (const Q& c);          // multiplication with a constant
    template<typename Q>
    vec<T> & operator/= (const vec<Q>& m);     // element-wise division of two vectors

    // elementwise arithmetics-assignment op's
    template <typename Q, typename R>
    vec<T> &elementwise_map_assign(
            const Q &op,
            const vec<R> &rhs)
    {
        if (rhs.size() != this->size())
        {
            throw std::length_error("Cannot perform pairwise operations on vectors of different size.");
        }
        for (size_t i = 0; i < this->size(); i++)
        {
            (*this)[i] = op((*this)[i], rhs[i]);
        }
        return *this;
    }
    // scalar arithmetics-assignment op's
    template <typename Q, typename R>
    vec<T> &scalar_map_assign(
            const Q &op,
            const R &rhs)
    {
        for (size_t i = 0; i < this->size(); i++)
        {
            (*this)[i] = op((*this)[i], rhs);
        }
        return *this;
    }

};


template <typename Q, typename L, typename R>
auto scalar_map(const Q &op, const vec<L> &lhs, const R &rhs)
{
    vec<decltype(op(std::declval<L>(), std::declval<R>()))> res(lhs.size());
    for (size_t i = 0; i < lhs.size(); i++)
    {
        res[i] = op(lhs[i], rhs);
    }
    return res;
}


template <typename binaryOp, typename L, typename R>
auto elementwise(const binaryOp &op, const vec<L> &lhs, const vec<R> &rhs) {
    if (lhs.size() != rhs.size())
    {
        throw std::length_error("Cannot perform pairwise operations on vectors of different size.");
    }
    vec<decltype(op(std::declval<L>(), std::declval<R>()))> res(lhs.size());
    for (size_t i = 0; i < res.size(); i++)
    {
        res[i] = op(lhs[i], rhs[i]);
    }
    assert(res.size() == lhs.size());
    return res;
}

template <typename T, typename O>
auto transform_vec(const std::function<O(T)> &op, const vec<T> &vect) {
    vec<O> res(vect.size());
    for (size_t i = 0; i < vect.size(); i++)
    {
        res[i] = op(vect[i]);
    }
    assert(res.size() == vect.size());
    return res;
}

template <typename T, typename R>
void convert_vec_to_type(const vec<T> &vect, R& res) {
    for (size_t i = 0; i < vect.size(); i++)
    {
        res[i] = vect[i];
    }
}




// define aliases for real and complex vector
typedef vec<double> rvec;
typedef vec<comp> cvec;


/// DEFINITIONS -- MEMBER FUNCTIONS ///

// element-wise addition of two vectors
template <typename T> template<typename Q>
vec<T> & vec<T>::operator+= (const vec<Q>& m) {
    return elementwise_map_assign([](const T &l, const Q &r) { return l + r; }, m);
}

// addition of a constant
template <typename T> template<typename Q>
vec<T> & vec<T>::operator+= (const Q& c) {
    return scalar_map_assign([](const T &l, const Q &r) { return l + r; }, c);
}

// element-wise subtraction of two vectors
template <typename T> template<typename Q>
vec<T> & vec<T>::operator-= (const vec<Q>& m) {
    return elementwise_map_assign([](const T &l, const Q &r) { return l - r; }, m);
}

// subtraction of a constant
template <typename T> template<typename Q>
vec<T> & vec<T>::operator-= (const Q& c) {
    return scalar_map_assign([](const T &l, const Q &r) { return l - r; }, c);

}

// element-wise division of two vectors
template <typename T> template<typename Q>
vec<T> & vec<T>::operator/= (const vec<Q>& m) {
    return elementwise_map_assign([](const T &l, const Q &r) { return l / r; }, m);
}
// element-wise multiplication of two vectors
template <typename T> template<typename Q>
vec<T> & vec<T>::operator*= (const vec<Q>& m) {
    return elementwise_map_assign([](const T &l, const Q &r) { return l * r; }, m);

}

// multiplication with a constant
template <typename T> template<typename Q>
vec<T> & vec<T>::operator*= (const Q& c) {
    return scalar_map_assign([](const T &l, const Q &r) { return l * r; }, c);

}

// get a subvector {x[i1], ..., x[i2]}
// if indices are negative, count from the end
template <typename T>
vec<T> vec<T>::operator() (int i1, int i2) const {
    auto it1 = (i1 >= 0) ? this->begin() : this->end();
    auto it2 = (i2 >= 0) ? this->begin() : this->end();
    vec<T> subvector (it1 + i1, it2 + i2 + 1);
    return subvector;
}

template<typename L, typename R>
auto operator+ (const vec<L>& lhs, const vec<R>& rhs) { // element-wise addition of two vectors
    return elementwise([](const L &l, const R &r) { return l + r; }, lhs, rhs);
};
template<typename L, typename R>
auto operator+ (const vec<L>& lhs, const R& rhs) {      // addition of a constant
    return scalar_map([](const L &l, const R &r) { return l + r; }, lhs, rhs);
};
template<typename L, typename R>
auto operator+ (const L& lhs, const vec<R>& rhs) {      // addition of a constant from the left (commutative)
    return scalar_map([](const R &r, const L &l) { return l + r; }, rhs, lhs);
};
template<typename L, typename R>
auto operator- (const vec<L>& lhs, const vec<R>& rhs) { // element-wise subtraction of two vectors
    return elementwise([](const L &l, const R &r) { return l - r; }, lhs, rhs);
};
template<typename L, typename R>
auto operator- (const vec<L>& lhs, const R& rhs) {      // subtraction of a constant
    return scalar_map([](const L &l, const R &r) { return l - r; }, lhs, rhs);
};
template<typename L, typename R>
auto operator- (const L& lhs, const vec<R>& rhs) {      // subtraction of a constant from the left
    return scalar_map([](const R &r, const L &l) { return l - r; }, rhs, lhs);
};
template<typename L, typename R>
auto operator* (const vec<L> lhs, const vec<R>& rhs) { // element-wise multiplication of two vectors
    return elementwise([](const L &l, const R &r) { return l * r; }, lhs, rhs);
};
template<typename L, typename R>
auto operator* (const vec<L> lhs, const R& rhs) {      // scalar multiplication with a constant from the right
    return scalar_map([](const L &l, const R &r) { return l * r; }, lhs, rhs);
};
template<typename L, typename R>
auto operator* (const L& lhs, const vec<R> rhs) {      // scalar multiplication with a constant from the left
    return scalar_map([](const R &r, const L &l) { return l * r; }, rhs, lhs);
};
template<typename L, typename R>
auto operator/ (const vec<L> lhs, const vec<R>& rhs) { // element-wise division of two vectors
    return elementwise([](const L &l, const R &r) { return l / r; }, lhs, rhs);
};

// element-wise inverse
template <typename T>
vec<T> vec<T>::inv() const {
    return transform_vec<T,T>([](const T& x) -> T {return 1./x;}, *this);
}

// element-wise real part
template <typename T>
vec<double> vec<T>::real() const {                    // if T != comp or double, vector of zeros is returned
    return transform_vec<T,double>(myreal<T>, *this);
}

// element-wise imaginary part
template <typename T>
vec<double> vec<T>::imag() const {                    // if T != comp, vector of zeros is returned
    return transform_vec<T,double>(myimag<T>, *this);
}


// element-wise absolute value
template <typename T>
template<typename T2>
vec<T2> vec<T>::abs() const {
    vec<T2> temp (this->size());
#pragma omp parallel for
    for (unsigned int i=0; i<this->size(); ++i) {
        temp[i] = std::abs((*this)[i]);
    }
    return temp;
}

// element-wise complex conjugate
template <typename T>
vec<T> vec<T>::conj() const {                         // if T != comp, return input
    return transform_vec<T,T>(myconj<T>, *this);
}


// maximum norm
template <typename T>
double vec<T>::max_norm() const {
    double out = 0.;
    for (size_t i=0; i<this->size(); ++i) {
        out = std::max(out, (double)std::abs((*this)[i]));
    }
    return out;
}

// vector of differences between adjacent elements
template <typename T>
vec<T> vec<T>::diff() const {
    vec<T> xp (this->begin() + 1, this->end()); // second -> last element
    vec<T> xm (this->begin(), this->end() - 1); // first -> second to last element
    return xp - xm;                             // compute differences
}

// sum of all elements
template <typename T>
T vec<T>::sum() const {
    return std::accumulate(this->begin(), this->end(), myzero<T>());
}

/// NON-MEMBER FUNCTIONS ///

// The functions below are necessary for operations concerning a complex vector and a double constant.

// addition of a double constant to comp vector
template <typename T>
vec<T> operator+= (vec<T>& lhs, const double& rhs) {
#pragma omp parallel for
    for (unsigned int i=0; i<lhs.size(); ++i) {
        lhs[i] += rhs;
    }
    return lhs;
}
template <typename T>
vec<T> operator+ (vec<T> lhs, const double& rhs) {
    lhs += rhs; return lhs;
}
template <typename T>
vec<T> operator+ (const double& rhs, vec<T> lhs) {
    lhs += rhs; return lhs;
}

// subtraction of a double constant to comp vector
template <typename T>
vec<T> operator-= (vec<T>& lhs, const double& rhs) {
//#pragma omp parallel for
    for (size_t i=0; i<lhs.size(); ++i) {
        lhs[i] -= rhs;
    }
    return lhs;
}
template <typename T>
vec<T> operator- (vec<T> lhs, const double& rhs) {
    lhs -= rhs; return lhs;
}

// multiplication of a double constant to comp vector
template <typename T>
vec<T> operator*= (vec<T>& lhs, const double& rhs) {
#pragma omp parallel for
    for (unsigned int i=0; i<lhs.size(); ++i) {
        lhs[i] *= rhs;
    }
    return lhs;
}
template <typename T>
vec<T> operator* (vec<T> lhs, const double& rhs) {
    lhs *= rhs; return lhs;
}
template <typename T>
vec<T> operator* (const double& rhs, vec<T> lhs) {
    lhs *= rhs; return lhs;
}


// multiplication of a comp constant to double vector
template <typename T>
vec<comp> operator* (vec<T> lhs, const comp& rhs) {
    size_t size = lhs.size();
    vec<comp> result(size);
#pragma omp parallel for
    for (int i=0; i<lhs.size(); ++i) {
        result[i] *= lhs[i] * rhs;
    }
    return result;
}







/// Keldysh index parameters ///



using my_index_t = std::size_t;
enum K_class {selfenergy = -1, k1=0, k2=1, k2b=2, k3=3, k3_sbe=4};
std::ostream& operator << (std::ostream& out, K_class k);
enum vertexType {symmetric_full, symmetric_r_irred, non_symmetric_diffleft, non_symmetric_diffright};
std::ostream& operator << (std::ostream& out, vertexType symmtype);

constexpr size_t rank_SE = 3;
constexpr size_t rank_K1 = 4;
constexpr size_t rank_K2 = 5;
constexpr size_t rank_K3 = 6;
namespace my_defs {
    namespace SE {
        constexpr my_index_t rank = 3;
        enum names {keldysh, nu, internal};
        using index_type = std::array<my_index_t, rank>;
        using dimensions_type = std::array<my_index_t, rank>;
        using frequencies_type = std::array<freqType, 1>;
    }
    namespace K1 {
        constexpr my_index_t rank = 4;
        enum names {spin,  omega,  keldysh,  internal};
        using index_type = std::array<my_index_t, rank>;
        using dimensions_type = std::array<my_index_t, rank>;
        using frequencies_type = std::array<freqType, 1>;
    }
    namespace K2 {
        constexpr my_index_t rank = 5;
        enum names {spin,  omega,  nu,  keldysh,  internal};
        using index_type = std::array<my_index_t, rank>;
        using dimensions_type = std::array<my_index_t, rank>;
        using frequencies_type = std::array<freqType, 2>;
    }
    namespace K2b{
        constexpr my_index_t rank = 5;
        enum names {spin,  omega,  nup,  keldysh,  internal};
        using index_type = std::array<my_index_t, rank>;
        using dimensions_type = std::array<my_index_t, rank>;
        using frequencies_type = std::array<freqType, 2>;
    }
    namespace K3 {
        constexpr  my_index_t rank = 6;
        enum names {spin,  omega,  nu,  nup,  keldysh,  internal};
        using index_type = std::array<my_index_t, rank>;
        using dimensions_type = std::array<my_index_t, rank>;
        using frequencies_type = std::array<freqType, 3>;
    }
}



/** auxiliary struct that contains all input variables of vertices
 * @param iK       :   integer from 0 to 15 (Keldysh indices expressed as one integer)
 * @param w.v1,v2  :   frequency arguments
 * @param i_in     :   additional internal index (currently unused)
 * @param spin     :   0 or 1 for the two comfigurations (0 for V, 1 for V^)
 * @param channel  :   'a', 't' or 'p', or 'f'
 * */
struct VertexInput{
    int iK;
    freqType w, v1, v2;
    my_index_t i_in;
    my_index_t spin;
    char channel_bubble;
    K_class kClass_aim;
    my_index_t iw_r;

    VertexInput(int iK_in, my_index_t spin_in, freqType w_in, freqType v1_in, freqType v2_in, my_index_t i_in_in, char channel_in, K_class k_in=k1, my_index_t iw_in=0)
            :
//#if KELDYSH_FORMALISM
            iK(iK_in),
//#else
//            iK(0),
//#endif
            w(w_in), v1(v1_in), v2(v2_in), i_in(i_in_in), spin(spin_in), channel_bubble(channel_in), kClass_aim(k_in), iw_r(iw_in)
    {assert(iK < 16);}

    template<K_class k>
    auto get_freqs() const {
        if constexpr(k == k1) {
            my_defs::K1::frequencies_type freqs = {w};
            return freqs;
        }
        else
        if constexpr(k == k2) {
            my_defs::K2::frequencies_type freqs = {w, v1};
            return freqs;
        }
        else
        if constexpr(k == k2b) {
            my_defs::K2b::frequencies_type freqs = {w, v2};
            return freqs;
        }
        else
        if constexpr(k == k3) {
            my_defs::K3::frequencies_type freqs = {w, v1, v2};
            return freqs;
        }
        else assert(false);
    }

    template<K_class k>
    auto get_indices() const {
        if constexpr(k == k1) {
            my_defs::K1::index_type idx;
            idx[my_defs::K1::spin] = spin;
            idx[my_defs::K1::keldysh] = iK;
            idx[my_defs::K1::internal] = i_in;

            return idx;
        }
        else
        if constexpr(k == k2) {
            my_defs::K2::index_type idx;
            idx[my_defs::K2::spin] = spin;
            idx[my_defs::K2::keldysh] = iK;
            idx[my_defs::K2::internal] = i_in;

            return idx;
        }
        else
        if constexpr(k == k2b) {
            my_defs::K2b::index_type idx;
            idx[my_defs::K2b::spin] = spin;
            idx[my_defs::K2b::keldysh] = iK;
            idx[my_defs::K2b::internal] = i_in;
            return idx;
        }
        else
        if constexpr(k == k3) {
            my_defs::K3::index_type idx;
            idx[my_defs::K3::spin] = spin;
            idx[my_defs::K3::keldysh] = iK;
            idx[my_defs::K3::internal] = i_in;
            return idx;
        }
        else assert(false);
    }
};

/**
 * specifies the diagrammatic contribution ( + modifications by prefactor, complex conjugation)
 */
struct IndicesSymmetryTransformations: VertexInput{
    //int iK;
    double prefactor = 1.; // fermionic sign factor; comes in effect for T1, T2 (and sometimes Tc)
    bool conjugate = false;
    bool asymmetry_transform = false;
    //int spin;
    //int iw_r;
    //double w, v1, v2; int i_in;
    //K_class kClass_aim;     // we only distinguish kClass==k3 from kClass!=k3 --> important for interpolation in K3
    char channel_rvert;
    //char channel_bubble;
    char channel_parametrization = channel_bubble; // W.r.t. which channel is the vertex parametrized? Used for the Hubbard model.

    IndicesSymmetryTransformations(int iK_in, my_index_t spin_in, double w_in, double v1_in, double v2_in, my_index_t i_in_in, char channel_rvert_in, K_class k_in, my_index_t iw_in, char channel_bubble_in)
            : VertexInput(iK_in, spin_in, w_in, v1_in, v2_in, i_in_in, channel_bubble_in, k_in, iw_in), channel_rvert(channel_rvert_in)
    {assert(iK < 16);}

    IndicesSymmetryTransformations(VertexInput input, char channel_in)
            : VertexInput(input), channel_rvert(channel_in)
    {assert(iK < 16);}
    IndicesSymmetryTransformations(VertexInput&& input, char channel_in)
            : VertexInput(input), channel_rvert(channel_in)
    {assert(iK < 16);}
};



struct fRG_config {
    int nODE_;
    double epsODE_abs_;
    double epsODE_rel_;
    double U = 1.;           // Impurity on-site interaction strength
    double T = 1.;           // temperature
    double Gamma = 1.;       // hybridization strength
    double epsilon = -0.5 * U;     // Impurity on-site energy
    bool save_intermediateResults = false;
    int nloops=1;
    int number_of_nodes=1;
};



#endif // DATA_STRUCTURES_H