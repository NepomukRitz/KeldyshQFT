#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include "frequency_grid.h"
#include "parameters.h"
#include "data_structures.h"
#include "symmetry_transformations.h"

//TODO improve to return the edge values
//TODO: references to indices instead of copy (for speed)??

// forward declaration of rvert from r_vertex.h
template <typename Q> class rvert;

/* linearly interpolate vertices */
template <K_class k, typename Q>
class Interpolate {
public:
    const double small = 1e-12;
    /** Template class call operator: used for K2 and K2b. For K1 and K3: template specializations (below) */
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
//    assert(vertex.frequencies.b_K2.w_lower <= w && w <= vertex.frequencies.b_K2.w_upper); // give error message if w out of range
//    assert(vertex.frequencies.f_K2.w_lower <= v && v <= vertex.frequencies.f_K2.w_upper); // give error message if v out of range

        double tw = vertex.frequencies.b_K2.grid_transf(indices.w);
        double tv1= vertex.frequencies.f_K2.grid_transf(indices.v1);
        if (    fabs(tw) < vertex.frequencies.b_K2.W_upper + this->small
             && fabs(tv1)< vertex.frequencies.f_K2.W_upper + this->small ) {

#if INTERPOLATION == 0
            int index_b = (int) ((tw - vertex.frequencies.b_K2.W_lower) / vertex.frequencies.b_K2.dW);
            int index_f = (int) ((tv1- vertex.frequencies.f_K2.W_lower) / vertex.frequencies.f_K2.dW);
            index_b = min(nBOS2-2, index_b);
            index_f = min(nFER2-2, index_f);
            assert(index_b >= 0 and index_b < nBOS2-1);
            assert(index_f >= 0 and index_f < nFER2-1);


            double x1;
            double x2;
            double y1;
            double y2;
            double xd;
            double yd;
            if (index_b == 0 or index_b == nBOS2-2) {
                x1 = vertex.frequencies.b_K2.Ws[index_b];
                x2 = vertex.frequencies.b_K2.Ws[index_b + 1];
                xd = (tw - x1) / (x2 - x1);
            }
            else {
                x1 = vertex.frequencies.b_K2.Ws[index_b];
                x2 = vertex.frequencies.b_K2.Ws[index_b + 1];
                xd = (indices.w - x1) / (x2 - x1);
            }
            if (index_f == 0 or index_f == nFER2-2) {
                y1 = vertex.frequencies.f_K2.Ws[index_f];
                y2 = vertex.frequencies.f_K2.Ws[index_f + 1];
                yd = (tv1- y1) / (y2 - y1);
            }
            else {
                y1 = vertex.frequencies.f_K2.Ws[index_f];
                y2 = vertex.frequencies.f_K2.Ws[index_f + 1];
                yd = (indices.v1- y1) / (y2 - y1);
            }

#elif INTERPOLATION>=1
            int index_b = (int) ((tw - vertex.frequencies.b_K2.W_lower) / vertex.frequencies.b_K2.dW);
            int index_f = (int) ((tv1- vertex.frequencies.f_K2.W_lower) / vertex.frequencies.f_K2.dW);
            index_b = min(nBOS2-2, index_b);
            index_f = min(nFER2-2, index_f);
            assert(index_b >= 0 and index_b < nBOS2-1);
            assert(index_f >= 0 and index_f < nFER2-1);


            double x1 = vertex.frequencies.b_K2.Ws[index_b];
            double x2 = vertex.frequencies.b_K2.Ws[index_b + 1];
            double y1 = vertex.frequencies.f_K2.Ws[index_f];
            double y2 = vertex.frequencies.f_K2.Ws[index_f + 1];
            double xd = (tw - x1) / (x2 - x1);
            double yd = (tv1- y1) / (y2 - y1);
#elif INTERPOLATION==2
            //TODO: write biquadratic interpolation

#elif INTERPOLATION==3
            //TODO: write bicubic interpolation
#endif


            auto f11 = vertex.K2_val(indices.iK, index_b, index_f, indices.i_in);
            auto f12 = vertex.K2_val(indices.iK, index_b, index_f + 1, indices.i_in);
            auto f21 = vertex.K2_val(indices.iK, index_b + 1, index_f, indices.i_in);
            auto f22 = vertex.K2_val(indices.iK, index_b + 1, index_f + 1, indices.i_in);

            Q result = indices.prefactor * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
            assert(isfinite(result));
            return indices.prefactor * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
        }
        else {
            return 0.;
        }
    }
};

/** Template specialization for K1 */
template <typename Q>
class Interpolate<k1, Q> {
public:
    const double small = 1e-12;
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {
//    assert(vertex.frequencies.b_K1.w_lower <= w && w <= vertex.frequencies.b_K1.w_upper); // give error message if w out of range
        double tw = vertex.frequencies.b_K1.grid_transf(indices.w);
        if (fabs(tw) < vertex.frequencies.b_K1.W_upper + this->small) {
            int index = (int) ((tw - vertex.frequencies.b_K1.W_lower) / vertex.frequencies.b_K1.dW);
            //if (index < 0 or index >= nBOS-1) { // If we get close to the boundaries of the box, make sure to stay within the box.
            //    index = vertex.frequencies.b_K1.fconv(vertex.frequencies.b_K1.w_upper*sign(indices.w), -sign(indices.w)*1e-1);
            //}

#if INTERPOLATION == 0
            index = min(nBOS-2, index);
            assert(index >= 0 and index < nBOS-1);

            double x1;
            double x2;
            double xd;
            if (index == 0 or index == nBOS-2) {
                x1 = vertex.frequencies.b_K1.Ws[index];
                x2 = vertex.frequencies.b_K1.Ws[index + 1];
                xd = (tw - x1) / (x2 - x1);
            }
            else {
                x1 = vertex.frequencies.b_K1.w[index];
                x2 = vertex.frequencies.b_K1.w[index + 1];
                xd = (indices.w - x1) / (x2 - x1);
            }

            auto f1 = vertex.K1_val(indices.iK, index, indices.i_in);
            auto f2 = vertex.K1_val(indices.iK, index + 1, indices.i_in);

            Q result = indices.prefactor * ((1. - xd) * f1 + xd * f2);
#elif INTERPOLATION==1
            index = min(nBOS-2, index);
            assert(index >= 0 and index < nBOS-1);
            double x1 = vertex.frequencies.b_K1.Ws[index];
            double x2 = vertex.frequencies.b_K1.Ws[index + 1];
            double xd = (tw - x1) / (x2 - x1);

            auto f1 = vertex.K1_val(indices.iK, index, indices.i_in);
            auto f2 = vertex.K1_val(indices.iK, index + 1, indices.i_in);

            Q result = indices.prefactor * ((1. - xd) * f1 + xd * f2);
#elif INTERPOLATION==2
            index = min(nBOS-2, index);
            index = max(1, index);
            assert(index >= 0 and index < nBOS-1);

            double x0 = vertex.frequencies.b_K1.Ws[index - 1];
            double x1 = vertex.frequencies.b_K1.Ws[index];
            double x2 = vertex.frequencies.b_K1.Ws[index + 1];

            auto f0 = vertex.K1_val(indices.iK, index - 1, indices.i_in);
            auto f1 = vertex.K1_val(indices.iK, index, indices.i_in);
            auto f2 = vertex.K1_val(indices.iK, index + 1, indices.i_in);

            Q result = indices.prefactor * ((tw - x1)*(tw - x2)/((x0-x1)*(x0-x2)) * f0
                                         +  (tw - x0)*(tw - x2)/((x1-x0)*(x1-x2)) * f1
                                         +  (tw - x0)*(tw - x1)/((x2-x0)*(x2-x1)) * f2);
#elif INTERPOLATION==3

            index = min(nBOS-3, index);
            index = max(1, index);
            assert(index >= 0 and index < nBOS-2);

            double x0 = vertex.frequencies.b_K1.Ws[index - 1];
            double x1 = vertex.frequencies.b_K1.Ws[index    ];
            double x2 = vertex.frequencies.b_K1.Ws[index + 1];
            double x3 = vertex.frequencies.b_K1.Ws[index + 2];

            auto f0 = vertex.K1_val(indices.iK, index - 1, indices.i_in);
            auto f1 = vertex.K1_val(indices.iK, index    , indices.i_in);
            auto f2 = vertex.K1_val(indices.iK, index + 1, indices.i_in);
            auto f3 = vertex.K1_val(indices.iK, index + 2, indices.i_in);

            Q result = indices.prefactor * ((tw - x1)*(tw - x2)*(tw - x3)/((x0-x1)*(x0-x2)*(x0-x3)) * f0
                                           +(tw - x0)*(tw - x2)*(tw - x3)/((x1-x0)*(x1-x2)*(x1-x3)) * f1
                                           +(tw - x0)*(tw - x1)*(tw - x3)/((x2-x0)*(x2-x1)*(x2-x3)) * f2
                                           +(tw - x0)*(tw - x1)*(tw - x2)/((x3-x0)*(x3-x1)*(x3-x2)) * f3);
#endif



            assert(isfinite(result));
            return result;
        }
        else {
            return 0.;
        }
    };
};

/** Template specialization for K3 */
template <typename Q>
class Interpolate<k3, Q> {
public:
    const double small = 1e-12;
    auto operator() (IndicesSymmetryTransformations& indices, const rvert<Q>& vertex) -> Q {

        double tw = vertex.frequencies.b_K3.grid_transf(indices.w);
        double tv1= vertex.frequencies.f_K3.grid_transf(indices.v1);
        double tv2= vertex.frequencies.f_K3.grid_transf(indices.v2);
        if (    fabs(tw) < vertex.frequencies.b_K3.W_upper + this->small
                && fabs(tv1)< vertex.frequencies.f_K3.W_upper + this->small
                && fabs(tv2)< vertex.frequencies.f_K3.W_upper + this->small) {

#if INTERPOLATION == 0
            int index_b = (int) ((tw - vertex.frequencies.b_K3.W_lower) / vertex.frequencies.b_K3.dW);
            int index_f1= (int) ((tv1- vertex.frequencies.f_K3.W_lower) / vertex.frequencies.f_K3.dW);
            int index_f2= (int) ((tv2- vertex.frequencies.f_K3.W_lower) / vertex.frequencies.f_K3.dW);
            index_b = min(nBOS3-2, index_b );
            index_f1= min(nFER3-2, index_f1);
            index_f2= min(nFER3-2, index_f2);


            double x1;
            double x2;
            double y1;
            double y2;
            double z1;
            double z2;
            double xd;
            double yd;
            double zd;
            if (index_b == 0 or index_b == nBOS3-2) {
                x1 = vertex.frequencies.b_K3.Ws[index_b];
                x2 = vertex.frequencies.b_K3.Ws[index_b + 1];
                xd = (tw - x1) / (x2 - x1);
            }
            else {
                x1 = vertex.frequencies.b_K3.Ws[index_b];
                x2 = vertex.frequencies.b_K3.Ws[index_b + 1];
                xd = (indices.w - x1) / (x2 - x1);
            }
            if (index_f1 == 0 or index_f1 == nFER3-2) {
                y1 = vertex.frequencies.f_K3.Ws[index_f1];
                y2 = vertex.frequencies.f_K3.Ws[index_f1 + 1];
                yd = (tv1- y1) / (y2 - y1);
            }
            else {
                y1 = vertex.frequencies.f_K3.Ws[index_f1];
                y2 = vertex.frequencies.f_K3.Ws[index_f1 + 1];
                yd = (indices.v1- y1) / (y2 - y1);
            }
            if (index_f2 == 0 or index_f2 == nFER3-2) {
                z1 = vertex.frequencies.f_K3.Ws[index_f2];
                z2 = vertex.frequencies.f_K3.Ws[index_f2 + 1];
                zd = (tv2- z1) / (z2 - z1);
            }
            else {
                z1 = vertex.frequencies.f_K3.Ws[index_f2];
                z2 = vertex.frequencies.f_K3.Ws[index_f2 + 1];
                zd = (indices.v2- z1) / (z2 - z1);
            }

#elif INTERPOLATION>=1
            int index_b = (int) ((tw - vertex.frequencies.b_K3.W_lower) / vertex.frequencies.b_K3.dW);
            int index_f1= (int) ((tv1- vertex.frequencies.f_K3.W_lower) / vertex.frequencies.f_K3.dW);
            int index_f2= (int) ((tv2- vertex.frequencies.f_K3.W_lower) / vertex.frequencies.f_K3.dW);
            index_b = min(nBOS3-2, index_b );
            index_f1= min(nFER3-2, index_f1);
            index_f2= min(nFER3-2, index_f2);

            assert(index_b >= 0 and index_b <nBOS3-1);
            assert(index_f1>= 0 and index_f1<nFER3-1);
            assert(index_f2>= 0 and index_f2<nFER3-1);


            double x1 = vertex.frequencies.b_K3.Ws[index_b];
            double x2 = vertex.frequencies.b_K3.Ws[index_b + 1];
            double y1 = vertex.frequencies.f_K3.Ws[index_f1];
            double y2 = vertex.frequencies.f_K3.Ws[index_f1 + 1];
            double z1 = vertex.frequencies.f_K3.Ws[index_f2];
            double z2 = vertex.frequencies.f_K3.Ws[index_f2 + 1];

            double xd = (tw - x1) / (x2 - x1);
            double yd = (tv1 - y1) / (y2 - y1);
            double zd = (tv2 - z1) / (z2 - z1);
#elif INTERPOLATION==2
            //TODO: write biquadratic interpolation

#elif INTERPOLATION==3
            //TODO: write bicubic interpolation
#endif

            auto f111 = vertex.K3_val(indices.iK, index_b, index_f1, index_f2, indices.i_in);
            auto f112 = vertex.K3_val(indices.iK, index_b, index_f1, index_f2 + 1, indices.i_in);
            auto f121 = vertex.K3_val(indices.iK, index_b, index_f1 + 1, index_f2, indices.i_in);
            auto f122 = vertex.K3_val(indices.iK, index_b, index_f1 + 1, index_f2 + 1, indices.i_in);
            auto f211 = vertex.K3_val(indices.iK, index_b + 1, index_f1, index_f2, indices.i_in);
            auto f212 = vertex.K3_val(indices.iK, index_b + 1, index_f1, index_f2 + 1, indices.i_in);
            auto f221 = vertex.K3_val(indices.iK, index_b + 1, index_f1 + 1, index_f2, indices.i_in);
            auto f222 = vertex.K3_val(indices.iK, index_b + 1, index_f1 + 1, index_f2 + 1, indices.i_in);

            auto c00 = f111 * (1. - xd) + f211 * xd;
            auto c01 = f112 * (1. - xd) + f212 * xd;
            auto c10 = f121 * (1. - xd) + f221 * xd;
            auto c11 = f122 * (1. - xd) + f222 * xd;
            auto c0 = c00 * (1. - yd) + c10 * yd;
            auto c1 = c01 * (1. - yd) + c11 * yd;

            Q result = indices.prefactor * (c0 * (1. - zd) + c1 * zd);
            assert(isfinite(result));
            return result;
        }
        else {
            return 0.;
        }
    }
};

#endif //KELDYSH_MFRG_INTERPOLATIONS_H