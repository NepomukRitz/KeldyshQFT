//
// Created by Sa.Aguirre on 2/4/20.
//

#ifndef KELDYSH_MFRG_INTERPOLATIONS_H
#define KELDYSH_MFRG_INTERPOLATIONS_H

#include "frequency_grid.h"
#include "parameters.h"
#include "data_structures.h"

//TODO improve to return the edge values
template <typename Q, typename T>
void interpolateK1(Q& ans, double pf, int iK, double w, int i_in, T& vertex){
    assert(w_lower_b<=w && w <=w_upper_b);
    if(fabs(w-w_lower_b)<inter_tol)
        ans = vertex.K1_vval(iK, 0, i_in);
    else if(fabs(w-w_upper_b)<inter_tol)
        ans = vertex.K1_vval(iK, nBOS-1, i_in);
    else {
        int index = fconv_bos(w);
        double x1 = bfreqs[index];
        double x2 = bfreqs[index + 1];
        double xd = (w - x1) / (x2 - x1);

        Q f1 = vertex.K1_vval(iK, index, i_in);
        Q f2 = vertex.K1_vval(iK, index + 1, i_in);

        ans = pf * ((1. - xd) * f1 + xd * f2);
    }
}

template <typename Q, typename T>
void interpolateK2 (Q& ans, double pf, int iK, double w, double v, int i_in, T& vertex){
    assert(w_lower_b<=w && w <=w_upper_b);
    assert(w_lower_f<=v && v <=w_upper_f);

    int index_b = fconv_bos(w);
    int index_f = fconv_fer(v);

    double x1 = bfreqs[index_b];
    double x2 = bfreqs[index_b + 1];
    double y1 = ffreqs[index_f];
    double y2 = ffreqs[index_f + 1];
    double xd = (w - x1) / (x2 - x1);
    double yd = (v - y1) / (y2 - y1);

    Q f11 = vertex.K2_vval(iK, index_b, index_f, i_in);
    Q f12 = vertex.K2_vval(iK, index_b, index_f + 1, i_in);
    Q f21 = vertex.K2_vval(iK, index_b + 1, index_f, i_in);
    Q f22 = vertex.K2_vval(iK, index_b + 1, index_f + 1, i_in);

    ans = pf * ((1. - yd) * ((1. - xd) * f11 + xd * f21) + yd * ((1. - xd) * f12 + xd * f22));
}

template <typename Q, typename T>
void interpolateK3 (Q& ans, double pf, int iK, double w, double v1, double v2, int i_in, T& vertex) {
    assert(w_lower_b<=w && w <=w_upper_b);
    assert(w_lower_f<=v1 && v1 <=w_upper_f);
    assert(w_lower_f<=v2 && v2 <=w_upper_f);

    int index_b = fconv_bos(w);
    int index_f1=fconv_fer(v1);
    int index_f2=fconv_fer(v2);

    double x1 = bfreqs[index_b];
    double x2 = bfreqs[index_b + 1];
    double y1 = ffreqs[index_f1];
    double y2 = ffreqs[index_f1 + 1];
    double z1 = ffreqs[index_f2];
    double z2 = ffreqs[index_f2 + 1];
    double xd = (w - x1) / (x2 - x1);
    double yd = (v1 - y1) / (y2 - y1);
    double zd = (v2- z1) / (z2 - z1);

    Q f111 = vertex.K3_vval(iK, index_b, index_f1, index_f2, i_in);
    Q f112 = vertex.K3_vval(iK, index_b, index_f1, index_f2 + 1, i_in);
    Q f121 = vertex.K3_vval(iK, index_b, index_f1 + 1, index_f2, i_in);
    Q f122 = vertex.K3_vval(iK, index_b, index_f1 + 1, index_f2 + 1, i_in);
    Q f211 = vertex.K3_vval(iK, index_b + 1, index_f1, index_f2, i_in);
    Q f212 = vertex.K3_vval(iK, index_b + 1, index_f1, index_f2 + 1, i_in);
    Q f221 = vertex.K3_vval(iK, index_b + 1, index_f1 + 1, index_f2, i_in);
    Q f222 = vertex.K3_vval(iK, index_b + 1, index_f1 + 1, index_f2 + 1, i_in);

    Q c00 = f111 * (1. - xd) + f211 * xd;
    Q c01 = f112 * (1. - xd) + f212 * xd;
    Q c10 = f121 * (1. - xd) + f221 * xd;
    Q c11 = f122 * (1. - xd) + f222 * xd;
    Q c0 = c00 * (1. - yd) + c10 * yd;
    Q c1 = c01 * (1. - yd) + c11 * yd;

    ans = pf * (c0 * (1. - zd) + c1 * zd);
}

#endif //KELDYSH_MFRG_INTERPOLATIONS_H
