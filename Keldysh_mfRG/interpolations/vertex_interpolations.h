#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

//#include <gsl/gsl_interp2d.h>
//#include <gsl/gsl_spline2d.h>
#include "InterpolatorSpline1D.h"
#include "InterpolatorSpline2D.h"
#include "InterpolatorSpline3D.h"
#include "../grids/frequency_grid.h"
#include "../parameters/master_parameters.h"
#include "../data_structures.h"
#include "../symmetries/symmetry_transformations.h"
#include "interpolation_functions.h"
#include "../vertex_data.h"
#include "../selfenergy.h"

#include "InterpolatorSpline1D.h"

// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;
// forward declaration of vertexInterpolator (see below)
template <typename Q> class vertexInterpolator;


/* linearly interpolate vertices */
template<int k, typename Q>
class Interpolate {
    explicit Interpolate(double Lambda) {
        assert(false);
    }
};

#if INTERPOLATION!=1 and INTERPOLATION!=3
namespace linOrSloppy {

    /* linearly interpolate vertices */
    template<int k, typename Q>
    class Interpolate {
        explicit Interpolate(double Lambda) {
            assert(false);
        }
    };
#endif

/** Template specialization for K1 */
    template<typename Q>
    class Interpolate<k1, Q> : public vertexDataContainer<k1, Q> {
    public:
        explicit Interpolate<k1, Q>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};

        void initializeK1() {};

        auto interpolK1(const IndicesSymmetryTransformations &indices) const -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
            //{
            Q result = indices.prefactor * interpolate1D<Q>(indices.w, vertexDataContainer<k1, Q>::frequencies_K1.b,
                                                            [&](int i) -> Q {
                                                                return vertexDataContainer<k1, Q>::K1_val(indices.iK, i,
                                                                                                          indices.i_in);
                                                            });
            // Lambda function (aka anonymous function) in last argument
            return result;
            //} else {
            //    return 0.;  // asymptotic value
            //}
        };
    };

    template<typename Q>
    class Interpolate<k2, Q> : public vertexDataContainer<k2, Q> {
    public:
        explicit Interpolate<k2, Q>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};

        void initializeK2() {};
        // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
        auto interpolK2(const IndicesSymmetryTransformations &indices) const -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            //if (    std::abs(indices.w ) < vertex.frequencies_K2.b.w_upper + inter_tol
            //        && std::abs(indices.v1) < vertex.frequencies_K2.f.w_upper + inter_tol )
            //{
            Q result = indices.prefactor * interpolate2D<Q>(indices.w, indices.v1,
                                                            vertexDataContainer<k2, Q>::frequencies_K2.b,
                                                            vertexDataContainer<k2, Q>::frequencies_K2.f,
                                                            [&](int i, int j) -> Q {
                                                                return vertexDataContainer<k2, Q>::K2_val(indices.iK, i,
                                                                                                          j,
                                                                                                          indices.i_in);
                                                            });
            return result;
            //}
            //else {
            //    return 0.;      // asymptotic value
            //}
        }
    };



/** Template specialization for K3 */
    template<typename Q>
    class Interpolate<k3, Q> : public vertexDataContainer<k3, Q> {
    public:
        void initializeK3() {};
        explicit Interpolate<k3, Q>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

        auto interpolK3(const IndicesSymmetryTransformations &indices) const -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol
            //    && std::abs(indices.v1) < vertex.frequencies_K3.f.w_upper + inter_tol
            //    && std::abs(indices.v2) < vertex.frequencies_K3.f.w_upper + inter_tol)
            //{
            Q result = indices.prefactor * interpolate3D<Q>(indices.w, indices.v1, indices.v2,
                                                            vertexDataContainer<k3, Q>::frequencies_K3.b,
                                                            vertexDataContainer<k3, Q>::frequencies_K3.f,
                                                            vertexDataContainer<k3, Q>::frequencies_K3.f,
                                                            [&](int i, int j, int k) -> Q {
                                                                return vertexDataContainer<k3, Q>::K3_val(indices.iK, i,
                                                                                                          j, k,
                                                                                                          indices.i_in);
                                                            });
            return result;
            //} else {
            //    return 0.;  // asymptotic value
            //}

        }
    };

#if INTERPOLATION!=1 and INTERPOLATION!=3
}
#endif



#if INTERPOLATION!=4
namespace spline {

    /* linearly interpolate vertices */
    template<int k, typename Q>
    class Interpolate {
        explicit Interpolate(double Lambda) {
            assert(false);
        }
    };
#endif

    template<typename Q>
    class Interpolate<k1,Q>: public SplineK1<vertexDataContainer<k1,Q>, Q> {
    public:
        explicit Interpolate<k1, Q>(double Lambda) : SplineK1<vertexDataContainer<k1,Q>, Q>(Lambda) {};
        auto interpolK1(const IndicesSymmetryTransformations &indices) const -> Q {

            // Check if the frequency runs out of the box; if yes: return asymptotic value
            //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
            //{
            Q result = indices.prefactor * SplineK1<vertexDataContainer<k1,Q>, Q>::interpolK1 (indices.iK, indices.w, indices.i_in);
            // Lambda function (aka anonymous function) in last argument
            return result;
            //} else {
            //    return 0.;  // asymptotic value
            //}
        };
    };



template<typename Q>
class Interpolate<k2,Q>: public SplineK2<vertexDataContainer<k2,Q>, Q> {
public:
    explicit Interpolate<k2, Q>(double Lambda) : SplineK2<vertexDataContainer<k2,Q>, Q>(Lambda) {};
    auto interpolK2(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K2.b.w_upper + inter_tol)
        //{
        Q result = indices.prefactor * SplineK2<vertexDataContainer<k2,Q>, Q>::interpolK2 (indices.iK, indices.w, indices.v1, indices.i_in);
        // Lambda function (aka anonymous function) in last argument
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };
};


template<typename Q>
class Interpolate<k3,Q>: public SplineK3<vertexDataContainer<k3,Q>, Q> {
public:
    explicit Interpolate<k3, Q>(double Lambda) : SplineK3<vertexDataContainer<k3,Q>, Q>(Lambda) {};
    auto interpolK3(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol)
        //{
        Q result = indices.prefactor * SplineK3<vertexDataContainer<k3,Q>, Q>::interpolK3 (indices.iK, indices.w, indices.v1, indices.v2, indices.i_in);
        // Lambda function (aka anonymous function) in last argument
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };
};


#if INTERPOLATION!=4
}
#endif


namespace {

    //template<typename Q> class vertexInterpolator;

    template<K_class k, typename Q>
    class Select {
    public:
        auto operator() (const IndicesSymmetryTransformations &indices, const vertexInterpolator<Q>& vertex) {
            return vertex.interpolK2(indices);
        }
    };

    template<typename Q>
    class Select<k1,Q> {
    public:
        auto operator() (const IndicesSymmetryTransformations &indices, const vertexInterpolator<Q>& vertex) {
            Q value = vertex.interpolK1(indices);
            assert(isfinite(value));
            return value;
        }
    };

    template<typename Q>
    class Select<k3,Q> {
    public:
        auto operator() (const IndicesSymmetryTransformations &indices, const vertexInterpolator<Q>& vertex) {
            return vertex.interpolK3(indices);
        }
    };


}


template <int k, typename Q> class vertexDataContainer; // forward declaration of vertexDataContainer

template<typename Q>
class vertexInterpolator: public Interpolate<k1,Q>,  public Interpolate<k2,Q>,  public Interpolate<k3,Q>  {
    friend fullvert<Q>;
    bool initialized = false;

public:

    explicit vertexInterpolator(double Lambda) : Interpolate<k1,Q>(Lambda), Interpolate<k2,Q>(Lambda), Interpolate<k3,Q>(Lambda) {};

    template <K_class k>
    auto interpolate(IndicesSymmetryTransformations& indices) const -> Q {

        assert(initialized);
        Q result = Select<k,Q>() (indices, *this);
        return result;
    }

    void initInterpolator() {
        if (not initialized) {      /// TODO: Is the member 'initialized' copied by the default copy constructor?
                                    ///       Or is it set to false?
            double t_start = get_time();
            //double tK3 = get_time() - t_start;
            Interpolate<k1,Q>::initializeK1();
            Interpolate<k2,Q>::initializeK2();
            Interpolate<k3,Q>::initializeK3();
            //print("K1 and K2 Interpolator initialized - ");
            //get_time(t_start);

            initialized = true;
        }
    }


};

#endif //KELDYSH_MFRG_INTERPOLATIONS_H