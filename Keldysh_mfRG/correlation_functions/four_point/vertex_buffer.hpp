#ifndef FPP_MFRG_VERTEX_BUFFER_H
#define FPP_MFRG_VERTEX_BUFFER_H

#include <cassert>
//#include "interpolations/vertex_interpolations.h"
#include "../../interpolations/InterpolatorSpline1D.hpp"
#include "../../interpolations/InterpolatorSpline2D.hpp"
#include "../../interpolations/InterpolatorSpline3D.hpp"
#include "../../interpolations/InterpolatorLinOrSloppy.hpp"
#include "../../symmetries/symmetry_table.hpp"
#include "../../symmetries/symmetry_transformations.hpp"


/**
 * Vertex buffers store the vertex data and frequency grids and interpolates data points on the grid.
 */
template<K_class k, typename Q, interpolMethod inter>
class vertexBuffer {
    explicit vertexBuffer(double Lambda) {
        assert(false);
    }
};


/** Template specialization for K1 (linear interpolation on the auxiliary grid) */
template<typename Q, interpolMethod inter>
class vertexBuffer<k1, Q, inter> : public vertexDataContainer<k1, Q> {
public:
    explicit vertexBuffer<k1, Q, inter>(double Lambda) : vertexDataContainer<k1, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {


        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (std::abs(indices.w) < vertexDataContainer<k1, Q>::frequencies_K1.b.w_upper + inter_tol)
        {

            Q result;
#ifdef DENSEGRID
            result =  interpolate_nearest1D<Q>(indices.w, vertexDataContainer<k1, Q>::get_VertexFreqGrid().b,
                                                  [&](int i) -> Q {
                                                      return vertexDataContainer<k1, Q>::val(indices.iK, indices.spin, i,
                                                                                             indices.i_in);
                                                  });
#else
            if constexpr(inter == linear_on_aux) {
                result =  interpolate_lin_on_aux1D<Q>(indices.w, vertexDataContainer<k1, Q>::get_VertexFreqGrid().b,
                                                        [&](int i) -> Q {
                                                            return vertexDataContainer<k1, Q>::val(indices.iK, indices.spin, i,
                                                                                                   indices.i_in);
                                                        });
            }
            else if constexpr(inter == linear) {
                result =  interpolate_lin1D<Q>(indices.w, vertexDataContainer<k1, Q>::get_VertexFreqGrid().b,
                                                      [&](int i) -> Q {
                                                          return vertexDataContainer<k1, Q>::val(indices.iK, indices.spin, i,
                                                                                                 indices.i_in);
                                                      });
            }
            else if constexpr(inter == sloppycubic) {
                result =  interpolate_sloppycubic1D<Q>(indices.w, vertexDataContainer<k1, Q>::get_VertexFreqGrid().b,
                                                   [&](int i) -> Q {
                                                       return vertexDataContainer<k1, Q>::val(indices.iK, indices.spin, i,
                                                                                              indices.i_in);
                                                   });
            }
            else assert(false);
            // Lambda function (aka anonymous function) in last argument
#endif

            assert(isfinite(result));
            return result;

        } else {
            return 0.;  // asymptotic value
        }

    };

    auto operator+= (const vertexBuffer<k1,Q,inter>& rhs) -> vertexBuffer<k1,Q,inter> {vertexDataContainer<k1,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k1,Q,inter>& rhs) -> vertexBuffer<k1,Q,inter> {vertexDataContainer<k1,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k1,Q,inter>& operator+ (vertexBuffer<k1,Q,inter>& lhs, const vertexBuffer<k1,Q,inter>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k1,Q,inter>& operator- (vertexBuffer<k1,Q,inter>& lhs, const vertexBuffer<k1,Q,inter>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

/** Template specialization for K2 (linear interpolation on the auxiliary grid) */
template<typename Q, interpolMethod inter>
class vertexBuffer<k2, Q, inter> : public vertexDataContainer<k2, Q> {
public:
    explicit vertexBuffer<k2, Q, inter>(double Lambda) : vertexDataContainer<k2, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        double w_temp = indices.w;
        double v_temp = indices.v1;
        K2_convert2internalFreqs(w_temp, v_temp); // convert natural frequency parametrization in channel r to internal parametrization


        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (    std::abs(indices.w ) < vertexDataContainer<k2, Q>::frequencies_K2.b.w_upper + inter_tol
                && std::abs(indices.v1) < vertexDataContainer<k2, Q>::frequencies_K2.f.w_upper + inter_tol )
        {

            Q result;
#ifdef DENSEGRID
            result =  interpolate_nearest2D<Q>(w_temp, v_temp,
                                               vertexDataContainer<k2, Q>::get_VertexFreqGrid().b,
                                               vertexDataContainer<k2, Q>::get_VertexFreqGrid().f,
                                               [&](int i, int j) -> Q {
                                                   return vertexDataContainer<k2, Q>::val(indices.iK, indices.spin, i,
                                                                                           j,
                                                                                           indices.i_in);
                                               });
#else
            if constexpr(inter == linear_on_aux) {
                result =  interpolate_lin_on_aux2D<Q>(w_temp, v_temp,
                                                      vertexDataContainer<k2, Q>::get_VertexFreqGrid().b,
                                                      vertexDataContainer<k2, Q>::get_VertexFreqGrid().f,
                                                        [&](int i, int j) -> Q {
                                                            return vertexDataContainer<k2, Q>::val(indices.iK, indices.spin, i,
                                                                                                   j,
                                                                                                   indices.i_in);
                                                        });
            }
            else if constexpr(inter == linear) {
                result =  interpolate_lin2D<Q>(w_temp, v_temp,
                                               vertexDataContainer<k2, Q>::get_VertexFreqGrid().b,
                                               vertexDataContainer<k2, Q>::get_VertexFreqGrid().f,
                                                        [&](int i, int j) -> Q {
                                                            return vertexDataContainer<k2, Q>::val(indices.iK, indices.spin, i,
                                                                                                   j,
                                                                                                   indices.i_in);
                                                        });
            }
            else if constexpr(inter == sloppycubic) {

                result =  interpolate_sloppycubic2D<Q>(w_temp, v_temp,
                                                       vertexDataContainer<k2, Q>::get_VertexFreqGrid().b,
                                                       vertexDataContainer<k2, Q>::get_VertexFreqGrid().f,
                                                     [&](int i, int j) -> Q {
                                                         return vertexDataContainer<k2, Q>::val(indices.iK, indices.spin, i,
                                                                                                j,
                                                                                                indices.i_in);
                                                     });
            }
            else assert(false);
#endif // DENSEGRID
            assert(isfinite(result));
            return result;

        }
        else {
            return 0.;      // asymptotic value
        }

    }

    auto operator+= (const vertexBuffer<k2,Q,inter>& rhs) -> vertexBuffer<k2,Q,inter> {vertexDataContainer<k2,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k2,Q,inter>& rhs) -> vertexBuffer<k2,Q,inter> {vertexDataContainer<k2,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k2,Q,inter>& operator+ (vertexBuffer<k2,Q,inter>& lhs, const vertexBuffer<k2,Q,inter>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k2,Q,inter>& operator- (vertexBuffer<k2,Q,inter>& lhs, const vertexBuffer<k2,Q,inter>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

#ifdef DEBUG_SYMMETRIES
/** Template specialization for K2b(linear interpolation on the auxiliary grid) */
template<typename Q, interpolMethod inter>
class vertexBuffer<k2b, Q, inter> : public vertexDataContainer<k2b, Q> {
public:
    explicit vertexBuffer<k2b, Q, inter>(double Lambda) : vertexDataContainer<k2b, Q>(Lambda) {};
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    // Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below)
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        double w_temp = indices.w;
        double v_temp = indices.v2;
        K2_convert2internalFreqs(w_temp, v_temp); // convert natural frequency parametrization in channel r to internal parametrization


        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (    std::abs(indices.w ) < vertexDataContainer<k2b, Q>::frequencies_K2.b.w_upper + inter_tol
                && std::abs(indices.v2) < vertexDataContainer<k2b, Q>::frequencies_K2.f.w_upper + inter_tol )
        {

            Q result;

#ifdef DENSEGRID
            result =  interpolate_nearest2D<Q>(w_temp, v_temp,
                                               vertexDataContainer<k2b, Q>::get_VertexFreqGrid().b,
                                               vertexDataContainer<k2b, Q>::get_VertexFreqGrid().f,
                                               [&](int i, int j) -> Q {
                                                   return vertexDataContainer<k2b, Q>::val(indices.iK, indices.spin, i,
                                                                                           j,
                                                                                           indices.i_in);
                                               });
#else
            if constexpr(inter == linear_on_aux){
                result =  interpolate_lin_on_aux2D<Q>(w_temp, v_temp,
                                                      vertexDataContainer<k2b, Q>::get_VertexFreqGrid().b,
                                                      vertexDataContainer<k2b, Q>::get_VertexFreqGrid().f,
                                                      [&](int i, int j) -> Q {
                                                          return vertexDataContainer<k2b, Q>::val(indices.iK, indices.spin, i,
                                                                                                  j,
                                                                                                  indices.i_in);
                                                      });
            }
            else if constexpr(inter == linear) {
                result =  interpolate_lin2D<Q>(w_temp, v_temp,
                                               vertexDataContainer<k2b, Q>::get_VertexFreqGrid().b,
                                               vertexDataContainer<k2b, Q>::get_VertexFreqGrid().f,
                                                      [&](int i, int j) -> Q {
                                                          return vertexDataContainer<k2b, Q>::val(indices.iK, indices.spin, i,
                                                                                                  j,
                                                                                                  indices.i_in);
                                                      });
            }
            else if constexpr(inter == sloppycubic) {
                result =  interpolate_sloppycubic2D<Q>(w_temp, v_temp,
                                                       vertexDataContainer<k2b, Q>::get_VertexFreqGrid().b,
                                                       vertexDataContainer<k2b, Q>::get_VertexFreqGrid().f,
                                               [&](int i, int j) -> Q {
                                                   return vertexDataContainer<k2b, Q>::val(indices.iK, indices.spin, i,
                                                                                           j,
                                                                                           indices.i_in);
                                               });
            }
            else assert(false);
#endif // DENSEGRID
            assert(isfinite(result));
            return result;

        }
        else {
            return 0.;      // asymptotic value
        }

    }

    auto operator+= (const vertexBuffer<k2b,Q,inter>& rhs) -> vertexBuffer<k2b,Q,inter> {vertexDataContainer<k2b,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k2b,Q,inter>& rhs) -> vertexBuffer<k2b,Q,inter> {vertexDataContainer<k2b,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k2b,Q,inter>& operator+ (vertexBuffer<k2b,Q,inter>& lhs, const vertexBuffer<k2b,Q,inter>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k2b,Q,inter>& operator- (vertexBuffer<k2b,Q,inter>& lhs, const vertexBuffer<k2b,Q,inter>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};
#endif

/** Template specialization for K3 (linear interpolation on the auxiliary grid) */
template<typename Q, interpolMethod inter>
class vertexBuffer<k3, Q, inter> : public vertexDataContainer<k3, Q> {
public:
    mutable bool initialized = false;

    void initInterpolator() const {initialized = true;};
    explicit vertexBuffer<k3, Q, inter>(double Lambda) : vertexDataContainer<k3, Q>(Lambda) {};

    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {
        double w_temp = indices.w;
        double v_temp = indices.v1;
        double vp_temp= indices.v2;
        if (BOSONIC_PARAM_FOR_K3) {
            if (indices.channel == 'a') { switch2bosonicFreqs<'a'>(w_temp, v_temp, vp_temp); }
            else if (indices.channel == 'p') { switch2bosonicFreqs<'p'>(w_temp, v_temp, vp_temp); }
            else if (indices.channel == 't') { switch2bosonicFreqs<'t'>(w_temp, v_temp, vp_temp); }
        }

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (std::abs(indices.w) < vertexDataContainer<k3, Q>::frequencies_K3.b.w_upper + inter_tol
            && std::abs(indices.v1) < vertexDataContainer<k3, Q>::frequencies_K3.f.w_upper + inter_tol
            && std::abs(indices.v2) < vertexDataContainer<k3, Q>::frequencies_K3.f.w_upper + inter_tol)
        {

            Q result;


#ifdef DENSEGRID
            result =  interpolate_nearest3D<Q>(w_temp, v_temp, vp_temp,
                                               vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                               vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                               vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                               [&](int i, int j, int k) -> Q {
                                                   return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,
                                                                                          j, k,
                                                                                          indices.i_in);
                                               });
#else
            if (not INTERPOL2D_FOR_K3 or indices.kClass_aim != k3) {

                if constexpr(inter == linear_on_aux) {
                    result = interpolate_lin_on_aux3D<Q>(w_temp, v_temp, vp_temp,
                                                         vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                         vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                                         vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                                         [&](int i, int j, int k) -> Q {
                                                             return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,
                                                                                                    j, k,
                                                                                                    indices.i_in);
                                                         });
                }
                else if constexpr(inter == linear) {
                    result = interpolate_lin3D<Q>(w_temp, v_temp, vp_temp,
                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                                         [&](int i, int j, int k) -> Q {
                                                             return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,
                                                                                                    j, k,
                                                                                                    indices.i_in);
                                                         });
                }
                else if constexpr(inter == sloppycubic) {
                    result = interpolate_sloppycubic3D<Q>(w_temp, v_temp, vp_temp,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().f,
                                                  [&](int i, int j, int k) -> Q {
                                                      return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,
                                                                                             j, k,
                                                                                             indices.i_in);
                                                  });
                }
                else assert(false);


            }

            else {

                if constexpr(inter == linear_on_aux) {
                    switch (indices.channel_bubble) {
                        case 'a':
                            result = interpolate_lin_on_aux2D<Q>(v_temp, vp_temp,
                                                                 vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, indices.iw,j, k,indices.i_in);});
                            break;
                        case 'p':
                            result = interpolate_lin_on_aux2D<Q>(w_temp, vp_temp,
                                                                 vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int i, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i, indices.iw, k,indices.i_in);});

                            break;
                        case 't':

                            result = interpolate_lin_on_aux2D<Q>(w_temp, v_temp,
                                                                 vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int i, int j) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,j, indices.iw,indices.i_in);});
                            break;
                        default:;
                    }
                }
                else if constexpr(inter == linear) {
                    switch (indices.channel_bubble) {
                        case 'a':
                            result = interpolate_lin2D<Q>(v_temp, vp_temp,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, indices.iw,j, k,indices.i_in);});
                            break;
                        case 'p':
                            result = interpolate_lin2D<Q>(w_temp, vp_temp,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int i, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i, indices.iw, k,indices.i_in);});
                            break;
                        case 't':

                            result = interpolate_lin2D<Q>(w_temp, v_temp,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                          vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int i, int j) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,j, indices.iw,indices.i_in);});
                            break;
                        default:;
                    }
                }
                else if constexpr(inter == sloppycubic) {
                    switch (indices.channel_bubble) {
                        case 'a':
                            result = interpolate_sloppycubic2D<Q>(v_temp, vp_temp,
                                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int j, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, indices.iw,j, k,indices.i_in);});
                            break;
                        case 'p':
                            result = interpolate_sloppycubic2D<Q>(w_temp, vp_temp,
                                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int i, int k) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i, indices.iw, k,indices.i_in);});

                            break;
                        case 't':
                            result = interpolate_sloppycubic2D<Q>(w_temp, v_temp,
                                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                  vertexDataContainer<k3, Q>::get_VertexFreqGrid().b,
                                                                 [&](int i, int j) -> Q {return vertexDataContainer<k3, Q>::val(indices.iK, indices.spin, i,j, indices.iw,indices.i_in);});
                            break;
                        default:;
                    }
                }
                else assert(false);
            }
#endif // DENSEGRID
            assert(isfinite(result));
            return result;


        } else {
            return 0.;  // asymptotic value
        }

    }

    auto operator+= (const vertexBuffer<k3,Q,inter>& rhs) -> vertexBuffer<k3,Q,inter> {vertexDataContainer<k3,Q>::data += rhs.data; return *this;}
    auto operator-= (const vertexBuffer<k3,Q,inter>& rhs) -> vertexBuffer<k3,Q,inter> {vertexDataContainer<k3,Q>::data -= rhs.data; return *this;}
    friend vertexBuffer<k3,Q,inter>& operator+ (vertexBuffer<k3,Q,inter>& lhs, const vertexBuffer<k3,Q,inter>& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vertexBuffer<k3,Q,inter>& operator- (vertexBuffer<k3,Q,inter>& lhs, const vertexBuffer<k3,Q,inter>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};



/** Template specialization for K1 (cubic interpolation) */
template<typename Q>
class vertexBuffer<k1,Q, cubic>: public SplineK1<vertexDataContainer<k1,Q>, Q> {

public:
    explicit vertexBuffer<k1, Q, cubic>(double Lambda) : SplineK1<vertexDataContainer<k1,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {
        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K1.b.w_upper + inter_tol)
        //{
        Q result =  SplineK1<vertexDataContainer<k1,Q>, Q>::interpolK1 (indices.iK, indices.spin, indices.w, indices.i_in);
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

/** Template specialization for K2 (cubic interpolation) */
template<typename Q>
class vertexBuffer<k2,Q, cubic>: public SplineK2<vertexDataContainer<k2,Q>, Q> {
public:
    explicit vertexBuffer<k2, Q, cubic>(double Lambda) : SplineK2<vertexDataContainer<k2,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        double w_temp = indices.w;
        double v_temp = indices.v1;
        K2_convert2internalFreqs(w_temp, v_temp); // convert natural frequency parametrization in channel r to internal parametrization

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K2.b.w_upper + inter_tol)
        //{

        Q result =  SplineK2<vertexDataContainer<k2,Q>, Q>::interpolK2 (indices.iK, indices.spin, w_temp, v_temp, indices.i_in);
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

/** Template specialization for K3 (cubic interpolation) */
template<typename Q>
class vertexBuffer<k3,Q, cubic>: public SplineK3<vertexDataContainer<k3,Q>, Q> {
public:
    explicit vertexBuffer<k3, Q, cubic>(double Lambda) : SplineK3<vertexDataContainer<k3,Q>, Q>(Lambda) {};
    auto interpolate(const IndicesSymmetryTransformations &indices) const -> Q {

        double w_temp = indices.w;
        double v_temp = indices.v1;
        double vp_temp= indices.v2;
#ifdef BOSONIC_PARAM_FOR_K3
        if (indices.channel == 'a') {switch2bosonicFreqs<'a'>(w_temp, v_temp, vp_temp);}
        else if (indices.channel == 'p') {switch2bosonicFreqs<'p'>(w_temp, v_temp, vp_temp);}
        else if (indices.channel == 't') {switch2bosonicFreqs<'t'>(w_temp, v_temp, vp_temp);}
#endif
        // Check if the frequency runs out of the box; if yes: return asymptotic value
        //if (std::abs(indices.w) < vertex.frequencies_K3.b.w_upper + inter_tol)
        //{
        Q result =  SplineK3<vertexDataContainer<k3,Q>, Q>::interpolK3 (indices.iK, indices.spin, w_temp, v_temp, vp_temp, indices.i_in);
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






#endif //FPP_MFRG_VERTEX_BUFFER_H
