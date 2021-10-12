#ifndef FPP_MFRG_VERTEX_BUFFER_H
#define FPP_MFRG_VERTEX_BUFFER_H

#include <cassert>
//#include "interpolations/vertex_interpolations.h"
#include "interpolations/InterpolatorSpline1D.h"
#include "interpolations/InterpolatorSpline2D.h"
#include "interpolations/InterpolatorSpline3D.h"
#include "interpolations/InterpolatorLinOrSloppy.h"
#include "symmetries/symmetry_table.h"
#include "symmetries/symmetry_transformations.h"



template<K_class k, typename Q, interpolMethod inter>
class vertexBuffer {
    explicit vertexBuffer(double Lambda) {
        assert(false);
    }
};

template<typename Q>
class vertexBuffer<k1,Q, cubic>: public SplineK1<vertexDataContainer<k1,Q>, Q> {

public:
    explicit vertexBuffer<k1, Q, cubic>(double Lambda) : SplineK1<vertexDataContainer<k1,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {
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
class vertexBuffer<k2,Q, cubic>: public SplineK2<vertexDataContainer<k2,Q>, Q> {
public:
    explicit vertexBuffer<k2, Q, cubic>(double Lambda) : SplineK2<vertexDataContainer<k2,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

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
class vertexBuffer<k3,Q, cubic>: public SplineK3<vertexDataContainer<k3,Q>, Q> {
public:
    explicit vertexBuffer<k3, Q, cubic>(double Lambda) : SplineK3<vertexDataContainer<k3,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

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


/// Interpolator class


/** Template specialization for K1 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k1, Q, linear> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, linear>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    bool initialized = false;
    void initInterpolator() {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (std::abs(indices.w) < vertexDataContainer<k1, Q>::frequencies_K1.b.w_upper + inter_tol)
        {
        Q result = indices.prefactor * interpolate_lin1D<Q>(indices.w, vertexDataContainer<k1, Q>::K1_get_VertexFreqGrid().b,
                                                            [&](int i) -> Q {return vertexDataContainer<k1, Q>::val(indices.iK, i, indices.i_in);});
        // Lambda function (aka anonymous function) in last argument
        return result;
        } else {
            return 0.;  // asymptotic value
        }
    };
};

/** Template specialization for K2 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k2, Q, linear> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, linear>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    bool initialized = false;

    void initInterpolator() {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (    std::abs(indices.w ) < vertexDataContainer<k2, Q>::frequencies_K2.b.w_upper + inter_tol
                && std::abs(indices.v1) < vertexDataContainer<k2, Q>::frequencies_K2.f.w_upper + inter_tol )
        {
        Q result = indices.prefactor * interpolate_lin2D<Q>(indices.w, indices.v1,
                                                            vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().b,
                                                            vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().f,
                                                            [&](int i, int j) -> Q {return vertexDataContainer<k2, Q>::val(indices.iK, i, j, indices.i_in);});
        return result;
        }
        else {
            return 0.;      // asymptotic value
        }
    }
};


/** Template specialization for K3 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k3, Q, linear> : public vertexDataContainer<k3, Q> {
public:
    bool initialized = false;

    void initInterpolator() {initialized = true;};
    explicit vertexBuffer<k3, Q, linear>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (std::abs(indices.w) < vertexDataContainer<k2, Q>::frequencies_K3.b.w_upper + inter_tol
            && std::abs(indices.v1) < vertexDataContainer<k2, Q>::frequencies_K3.f.w_upper + inter_tol
            && std::abs(indices.v2) < vertexDataContainer<k2, Q>::frequencies_K3.f.w_upper + inter_tol)
        {
        Q result = indices.prefactor * interpolate_lin3D<Q>(indices.w, indices.v1, indices.v2,
                                                            vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().b,
                                                            vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                            vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                            [&](int i, int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, i, j, k, indices.i_in);});
        return result;
        } else {
            return 0.;  // asymptotic value
        }

    }
};



/** Template specialization for K1 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k1, Q, linear_on_aux> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, linear_on_aux>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    bool initialized = false;

    void initInterpolator() {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
        //{
        Q result = indices.prefactor * interpolate_lin_on_aux1D<Q>(indices.w, vertexDataContainer<k1, Q>::K1_get_VertexFreqGrid().b,
                                                                   [&](int i) -> Q {
                                                                       return vertexDataContainer<k1, Q>::val(indices.iK, i,
                                                                                                              indices.i_in);
                                                                   });
        // Lambda function (aka anonymous function) in last argument
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };
};

/** Template specialization for K2 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k2, Q, linear_on_aux> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, linear_on_aux>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    bool initialized = false;

    void initInterpolator() {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (    std::abs(indices.w ) < vertex.frequencies_K2.b.w_upper + inter_tol
        //        && std::abs(indices.v1) < vertex.frequencies_K2.f.w_upper + inter_tol )
        //{
        Q result = indices.prefactor * interpolate_lin_on_aux2D<Q>(indices.w, indices.v1,
                                                                   vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().b,
                                                                   vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().f,
                                                                   [&](int i, int j) -> Q {
                                                                       return vertexDataContainer<k2, Q>::val(indices.iK, i,
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


/** Template specialization for K3 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k3, Q, linear_on_aux> : public vertexDataContainer<k3, Q> {
public:
    bool initialized = false;

    void initInterpolator() {initialized = true;};
    explicit vertexBuffer<k3, Q, linear_on_aux>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol
        //    && std::abs(indices.v1) < vertex.frequencies_K3.f.w_upper + inter_tol
        //    && std::abs(indices.v2) < vertex.frequencies_K3.f.w_upper + inter_tol)
        //{
        Q result = indices.prefactor * interpolate_lin_on_aux3D<Q>(indices.w, indices.v1, indices.v2,
                                                                   vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().b,
                                                                   vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                                   vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                                   [&](int i, int j, int k) -> Q {
                                                                       return vertexDataContainer<k3, Q>::val(indices.iK, i,
                                                                                                              j, k,
                                                                                                              indices.i_in);
                                                                   });
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}

    }
};




/** Template specialization for K1 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k1, Q, sloppycubic> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, sloppycubic>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    bool initialized = false;

    void initInterpolator() {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
        //{
        Q result = indices.prefactor * interpolate_sloppycubic1D<Q>(indices.w, vertexDataContainer<k1, Q>::K1_get_VertexFreqGrid().b,
                                                                    [&](int i) -> Q {return vertexDataContainer<k1, Q>::val(indices.iK, i, indices.i_in);

                                                                    });
        // Lambda function (aka anonymous function) in last argument
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };
};

/** Template specialization for K2 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k2, Q, sloppycubic> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, sloppycubic>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    bool initialized = false;

    void initInterpolator() {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (    std::abs(indices.w ) < vertex.frequencies_K2.b.w_upper + inter_tol
        //        && std::abs(indices.v1) < vertex.frequencies_K2.f.w_upper + inter_tol )
        //{
        Q result = indices.prefactor * interpolate_sloppycubic2D<Q>(indices.w, indices.v1,
                                                                    vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().b,
                                                                    vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().f,
                                                                    [&](int i, int j) -> Q {return vertexDataContainer<k2, Q>::val(indices.iK, i, j, indices.i_in);});
        return result;
        //}
        //else {
        //    return 0.;      // asymptotic value
        //}
    }
};


/** Template specialization for K3 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k3, Q, sloppycubic> : public vertexDataContainer<k3, Q> {
public:
    bool initialized = false;

    void initInterpolator() {initialized = true;};
    explicit vertexBuffer<k3, Q, sloppycubic>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol
        //    && std::abs(indices.v1) < vertex.frequencies_K3.f.w_upper + inter_tol
        //    && std::abs(indices.v2) < vertex.frequencies_K3.f.w_upper + inter_tol)
        //{
        Q result = indices.prefactor * interpolate_sloppycubic3D<Q>(indices.w, indices.v1, indices.v2,
                                                                    vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().b,
                                                                    vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                                    vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                                    [&](int i, int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, i,j, k,indices.i_in);});
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}

    }
};



#endif //FPP_MFRG_VERTEX_BUFFER_H
