#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H


/// Interpolator class
template<int k, typename Q>
class Interpolate {
    explicit Interpolate(double Lambda) {
        assert(false);
    }
};

#include "InterpolatorSpline1D.h"
#include "InterpolatorSpline2D.h"
#include "InterpolatorSpline3D.h"
#include "InterpolatorLinOrSloppy.h"
#include "../grids/frequency_grid.h"
#include "../parameters/master_parameters.h"
#include "../data_structures.h"
#include "../symmetries/symmetry_transformations.h"
#include "../vertex_data.h"
#include "../selfenergy.h"


// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;




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
protected:
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
                                    /// This becomes relevant if I want to reuse interpolation coefficients.
                                    /// Currently initialization takes < 1 second
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