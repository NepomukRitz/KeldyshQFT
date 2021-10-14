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
        Q result =  SplineK1<vertexDataContainer<k1,Q>, Q>::interpolK1 (indices.iK, indices.w, indices.i_in);
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };

    auto operator+= (const vertexBuffer<k1,Q,cubic>& rhs) -> vertexBuffer<k1,Q,cubic> {SplineK1<vertexDataContainer<k1,Q>,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k1,Q,cubic>& rhs) -> vertexBuffer<k1,Q,cubic> {SplineK1<vertexDataContainer<k1,Q>,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k1,Q,cubic>& operator+ (vertexBuffer<k1,Q,cubic>& lhs, const vertexBuffer<k1,Q,cubic>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k1,Q,cubic>& operator- (vertexBuffer<k1,Q,cubic>& lhs, const vertexBuffer<k1,Q,cubic>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};



template<typename Q>
class vertexBuffer<k2,Q, cubic>: public SplineK2<vertexDataContainer<k2,Q>, Q> {
public:
    explicit vertexBuffer<k2, Q, cubic>(double Lambda) : SplineK2<vertexDataContainer<k2,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K2.b.w_upper + inter_tol)
        //{
        Q result =  SplineK2<vertexDataContainer<k2,Q>, Q>::interpolK2 (indices.iK, indices.w, indices.v1, indices.i_in);
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };

    auto operator+= (const vertexBuffer<k2,Q,cubic>& rhs) -> vertexBuffer<k2,Q,cubic> {SplineK1<vertexDataContainer<k2,Q>,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k2,Q,cubic>& rhs) -> vertexBuffer<k2,Q,cubic> {SplineK1<vertexDataContainer<k2,Q>,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k2,Q,cubic>& operator+ (vertexBuffer<k2,Q,cubic>& lhs, const vertexBuffer<k2,Q,cubic>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k2,Q,cubic>& operator- (vertexBuffer<k2,Q,cubic>& lhs, const vertexBuffer<k2,Q,cubic>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


template<typename Q>
class vertexBuffer<k3,Q, cubic>: public SplineK3<vertexDataContainer<k3,Q>, Q> {
public:
    explicit vertexBuffer<k3, Q, cubic>(double Lambda) : SplineK3<vertexDataContainer<k3,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol)
        //{
        Q result =  SplineK3<vertexDataContainer<k3,Q>, Q>::interpolK3 (indices.iK, indices.w, indices.v1, indices.v2, indices.i_in);
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };

    auto operator+= (const vertexBuffer<k3,Q,cubic>& rhs) -> vertexBuffer<k3,Q,cubic> {SplineK1<vertexDataContainer<k3,Q>,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k3,Q,cubic>& rhs) -> vertexBuffer<k3,Q,cubic> {SplineK1<vertexDataContainer<k3,Q>,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k3,Q,cubic>& operator+ (vertexBuffer<k3,Q,cubic>& lhs, const vertexBuffer<k3,Q,cubic>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k3,Q,cubic>& operator- (vertexBuffer<k3,Q,cubic>& lhs, const vertexBuffer<k3,Q,cubic>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


/// Interpolator class


/** Template specialization for K1 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k1, Q, linear> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, linear>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    mutable bool initialized = false;
    void initInterpolator() const {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (std::abs(indices.w) < vertexDataContainer<k1, Q>::frequencies_K1.b.w_upper + inter_tol)
        {
        Q result =  interpolate_lin1D<Q>(indices.w, vertexDataContainer<k1, Q>::K1_get_VertexFreqGrid().b,
                                                            [&](int i) -> Q {return vertexDataContainer<k1, Q>::val(indices.iK, i, indices.i_in);});
        // Lambda function (aka anonymous function) in last argument
        return result;
        } else {
            return 0.;  // asymptotic value
        }
    };

    auto operator+= (const vertexBuffer<k1,Q,linear>& rhs) -> vertexBuffer<k1,Q,linear> {vertexDataContainer<k1,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k1,Q,linear>& rhs) -> vertexBuffer<k1,Q,linear> {vertexDataContainer<k1,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k1,Q,linear>& operator+ (vertexBuffer<k1,Q,linear>& lhs, const vertexBuffer<k1,Q,linear>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k1,Q,linear>& operator- (vertexBuffer<k1,Q,linear>& lhs, const vertexBuffer<k1,Q,linear>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** Template specialization for K2 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k2, Q, linear> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, linear>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (    std::abs(indices.w ) < vertexDataContainer<k2, Q>::frequencies_K2.b.w_upper + inter_tol
                && std::abs(indices.v1) < vertexDataContainer<k2, Q>::frequencies_K2.f.w_upper + inter_tol )
        {
        Q result =  interpolate_lin2D<Q>(indices.w, indices.v1,
                                                            vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().b,
                                                            vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().f,
                                                            [&](int i, int j) -> Q {return vertexDataContainer<k2, Q>::val(indices.iK, i, j, indices.i_in);});
        return result;
        }
        else {
            return 0.;      // asymptotic value
        }
    }

    auto operator+= (const vertexBuffer<k2,Q,linear>& rhs) -> vertexBuffer<k2,Q,linear> {vertexDataContainer<k2,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k2,Q,linear>& rhs) -> vertexBuffer<k2,Q,linear> {vertexDataContainer<k2,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k2,Q,linear>& operator+ (vertexBuffer<k2,Q,linear>& lhs, const vertexBuffer<k2,Q,linear>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k2,Q,linear>& operator- (vertexBuffer<k2,Q,linear>& lhs, const vertexBuffer<k2,Q,linear>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


/** Template specialization for K3 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k3, Q, linear> : public vertexDataContainer<k3, Q> {
public:
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    explicit vertexBuffer<k3, Q, linear>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (std::abs(indices.w) < vertexDataContainer<k3, Q>::K3_get_wupper_b() + inter_tol
            && std::abs(indices.v1) < vertexDataContainer<k3, Q>::K3_get_wupper_f() + inter_tol
            && std::abs(indices.v2) < vertexDataContainer<k3, Q>::K3_get_wupper_f() + inter_tol)
        {
        Q result =  interpolate_lin3D<Q>(indices.w, indices.v1, indices.v2,
                                                            vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().b,
                                                            vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                            vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                            [&](int i, int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, i, j, k, indices.i_in);});
        return result;
        } else {
            return 0.;  // asymptotic value
        }

    }

    auto operator+= (const vertexBuffer<k3,Q,linear>& rhs) -> vertexBuffer<k3,Q,linear> {vertexDataContainer<k3,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k3,Q,linear>& rhs) -> vertexBuffer<k3,Q,linear> {vertexDataContainer<k3,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k3,Q,linear>& operator+ (vertexBuffer<k3,Q,linear>& lhs, const vertexBuffer<k3,Q,linear>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k3,Q,linear>& operator- (vertexBuffer<k3,Q,linear>& lhs, const vertexBuffer<k3,Q,linear>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};



/** Template specialization for K1 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k1, Q, linear_on_aux> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, linear_on_aux>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
        //{
        Q result =  interpolate_lin_on_aux1D<Q>(indices.w, vertexDataContainer<k1, Q>::K1_get_VertexFreqGrid().b,
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

    auto operator+= (const vertexBuffer<k1,Q,linear_on_aux>& rhs) -> vertexBuffer<k1,Q,linear_on_aux> {vertexDataContainer<k1,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k1,Q,linear_on_aux>& rhs) -> vertexBuffer<k1,Q,linear_on_aux> {vertexDataContainer<k1,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k1,Q,linear_on_aux>& operator+ (vertexBuffer<k1,Q,linear_on_aux>& lhs, const vertexBuffer<k1,Q,linear_on_aux>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k1,Q,linear_on_aux>& operator- (vertexBuffer<k1,Q,linear_on_aux>& lhs, const vertexBuffer<k1,Q,linear_on_aux>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** Template specialization for K2 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k2, Q, linear_on_aux> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, linear_on_aux>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (    std::abs(indices.w ) < vertex.frequencies_K2.b.w_upper + inter_tol
        //        && std::abs(indices.v1) < vertex.frequencies_K2.f.w_upper + inter_tol )
        //{
        Q result =  interpolate_lin_on_aux2D<Q>(indices.w, indices.v1,
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

    auto operator+= (const vertexBuffer<k2,Q,linear_on_aux>& rhs) -> vertexBuffer<k2,Q,linear_on_aux> {vertexDataContainer<k2,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k2,Q,linear_on_aux>& rhs) -> vertexBuffer<k2,Q,linear_on_aux> {vertexDataContainer<k2,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k2,Q,linear_on_aux>& operator+ (vertexBuffer<k2,Q,linear_on_aux>& lhs, const vertexBuffer<k2,Q,linear_on_aux>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k2,Q,linear_on_aux>& operator- (vertexBuffer<k2,Q,linear_on_aux>& lhs, const vertexBuffer<k2,Q,linear_on_aux>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


/** Template specialization for K3 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k3, Q, linear_on_aux> : public vertexDataContainer<k3, Q> {
public:
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    explicit vertexBuffer<k3, Q, linear_on_aux>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol
        //    && std::abs(indices.v1) < vertex.frequencies_K3.f.w_upper + inter_tol
        //    && std::abs(indices.v2) < vertex.frequencies_K3.f.w_upper + inter_tol)
        //{
        Q result =  interpolate_lin_on_aux3D<Q>(indices.w, indices.v1, indices.v2,
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

    auto operator+= (const vertexBuffer<k3,Q,linear_on_aux>& rhs) -> vertexBuffer<k3,Q,linear_on_aux> {vertexDataContainer<k3,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k3,Q,linear_on_aux>& rhs) -> vertexBuffer<k3,Q,linear_on_aux> {vertexDataContainer<k3,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k3,Q,linear_on_aux>& operator+ (vertexBuffer<k3,Q,linear_on_aux>& lhs, const vertexBuffer<k3,Q,linear_on_aux>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k3,Q,linear_on_aux>& operator- (vertexBuffer<k3,Q,linear_on_aux>& lhs, const vertexBuffer<k3,Q,linear_on_aux>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};




/** Template specialization for K1 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k1, Q, sloppycubic> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, sloppycubic>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
        //{
        Q result =  interpolate_sloppycubic1D<Q>(indices.w, vertexDataContainer<k1, Q>::K1_get_VertexFreqGrid().b,
                                                                    [&](int i) -> Q {return vertexDataContainer<k1, Q>::val(indices.iK, i, indices.i_in);

                                                                    });
        // Lambda function (aka anonymous function) in last argument
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}
    };

    auto operator+= (const vertexBuffer<k1,Q,sloppycubic>& rhs) -> vertexBuffer<k1,Q,sloppycubic> {vertexDataContainer<k1,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k1,Q,sloppycubic>& rhs) -> vertexBuffer<k1,Q,sloppycubic> {vertexDataContainer<k1,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k1,Q,sloppycubic>& operator+ (vertexBuffer<k1,Q,sloppycubic>& lhs, const vertexBuffer<k1,Q,sloppycubic>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k1,Q,sloppycubic>& operator- (vertexBuffer<k1,Q,sloppycubic>& lhs, const vertexBuffer<k1,Q,sloppycubic>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** Template specialization for K2 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k2, Q, sloppycubic> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, sloppycubic>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (    std::abs(indices.w ) < vertex.frequencies_K2.b.w_upper + inter_tol
        //        && std::abs(indices.v1) < vertex.frequencies_K2.f.w_upper + inter_tol )
        //{
        Q result = interpolate_sloppycubic2D<Q>(indices.w, indices.v1,
                                                                    vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().b,
                                                                    vertexDataContainer<k2, Q>::K2_get_VertexFreqGrid().f,
                                                                    [&](int i, int j) -> Q {return vertexDataContainer<k2, Q>::val(indices.iK, i, j, indices.i_in);});
        return result;
        //}
        //else {
        //    return 0.;      // asymptotic value
        //}
    }

    auto operator+= (const vertexBuffer<k2,Q,sloppycubic>& rhs) -> vertexBuffer<k2,Q,sloppycubic> {vertexDataContainer<k2,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k2,Q,sloppycubic>& rhs) -> vertexBuffer<k2,Q,sloppycubic> {vertexDataContainer<k2,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k2,Q,sloppycubic>& operator+ (vertexBuffer<k2,Q,sloppycubic>& lhs, const vertexBuffer<k2,Q,sloppycubic>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k2,Q,sloppycubic>& operator- (vertexBuffer<k2,Q,sloppycubic>& lhs, const vertexBuffer<k2,Q,sloppycubic>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


/** Template specialization for K3 (linear or sloppy cubic interpolation) */
template<typename Q>
class vertexBuffer<k3, Q, sloppycubic> : public vertexDataContainer<k3, Q> {
public:
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    explicit vertexBuffer<k3, Q, sloppycubic>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol
        //    && std::abs(indices.v1) < vertex.frequencies_K3.f.w_upper + inter_tol
        //    && std::abs(indices.v2) < vertex.frequencies_K3.f.w_upper + inter_tol)
        //{
        Q result = interpolate_sloppycubic3D<Q>(indices.w, indices.v1, indices.v2,
                                                                    vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().b,
                                                                    vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                                    vertexDataContainer<k3, Q>::K3_get_VertexFreqGrid().f,
                                                                    [&](int i, int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, i,j, k,indices.i_in);});
        return result;
        //} else {
        //    return 0.;  // asymptotic value
        //}

    }

    auto operator+= (const vertexBuffer<k3,Q,sloppycubic>& rhs) -> vertexBuffer<k3,Q,sloppycubic> {vertexDataContainer<k3,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k3,Q,sloppycubic>& rhs) -> vertexBuffer<k3,Q,sloppycubic> {vertexDataContainer<k3,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k3,Q,sloppycubic>& operator+ (vertexBuffer<k3,Q,sloppycubic>& lhs, const vertexBuffer<k3,Q,sloppycubic>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k3,Q,sloppycubic>& operator- (vertexBuffer<k3,Q,sloppycubic>& lhs, const vertexBuffer<k3,Q,sloppycubic>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};



#endif //FPP_MFRG_VERTEX_BUFFER_H
