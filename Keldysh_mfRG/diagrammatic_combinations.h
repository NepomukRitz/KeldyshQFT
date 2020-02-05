//
// Created by Elias Walter on 2020-01-15.
//

#ifndef KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
#define KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H

#include "vertex.h"



template <typename Q>
auto left_same_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in, char channel) -> Q
{
    Q gamma0, K1, K2b;
    switch (channel){
        case 'a':
            gamma0 = vertex.densvertex.irred.vval(i1, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.avertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.avertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.tvertex);
#endif
            break;

        case 'p':
            gamma0 = vertex.densvertex.irred.vval(i1, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.pvertex.K1_vvalsmooth(i1, w, i_in);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.pvertex.K2b_vvalsmooth(i1, w, vpp, i_in);
#endif
            break;
        case 't' :
            gamma0 = vertex.densvertex.irred.vval(i1, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.tvertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.tvertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;

    }
    return gamma0 + K1 + K2b;
}
template <typename Q>
auto right_same_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vpp, int i_in, char channel) -> Q
{
    Q gamma0, K1, K2b;
    switch (channel){
        case 'a':
            gamma0 = vertex.densvertex.irred.vval(i3, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.avertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.avertex.K2b_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.tvertex);
#endif
            break;

        case 'p':
            gamma0 = vertex.densvertex.irred.vval(i3, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.pvertex.K1_vvalsmooth(i3, w, i_in);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.pvertex.K2b_vvalsmooth(i3, w, vpp, i_in);
#endif
            break;
        case 't' :
            gamma0 = vertex.densvertex.irred.vval(i3, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.densvertex.tvertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.densvertex.tvertex.K2b_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;

    }
    return gamma0 + K1 + K2b;
}
template <typename Q>
auto left_diff_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in, char channel) -> Q {

    Q K2, K3, gammaRb;

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.avertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.tvertex);
            gammaRb = vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 'a');
#endif
#if DIAG_CLASS >=3
            K3 = vertex.densvertex.avertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.pvertex.K2_vvalsmooth(i1, w, v, i_in);
            gammaRb = vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 'p');
#endif
#if DIAG_CLASS >=3
            K3 = vertex.densvertex.pvertex.K3_vvalsmooth(i1, w, v, vpp, i_in);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex.densvertex.tvertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.avertex);
            gammaRb = vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 't');;
#endif
#if DIAG_CLASS >=3
            K3 = vertex.densvertex.tvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2 + K3 + gammaRb;
}
template <typename Q>
auto right_diff_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in, char channel) -> Q {

    Q K2b, K3, gammaRb;

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex.densvertex.avertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.tvertex);
            gammaRb = vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 'a');
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.densvertex.avertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex.densvertex.pvertex.K2b_vvalsmooth(i3, w, vp, i_in);
            gammaRb = vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 'p');
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.densvertex.pvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex.densvertex.tvertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.avertex);
            gammaRb = vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 't');
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.densvertex.tvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gammaRb;
}



template <typename Q>
auto left_same_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in, int spin, char channel) -> Q
{
    Q gamma0, K1, K2b;
    switch (channel){
        case 'a':
            gamma0 = vertex.spinvertex.irred.vval(i1, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.avertex.K1_vvalsmooth(i1, w, i_in, spin, vertex.spinvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.avertex.K2b_vvalsmooth(i1, w, vpp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;

        case 'p':
            gamma0 = vertex.spinvertex.irred.vval(i1, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.pvertex.K1_vvalsmooth(i1, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.pvertex.K2b_vvalsmooth(i1, w, vpp, i_in, spin);
#endif
            break;
        case 't' :
            gamma0 = vertex.spinvertex.irred.vval(i1, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.tvertex.K1_vvalsmooth(i1, w, i_in, spin, vertex.spinvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.tvertex.K2b_vvalsmooth(i1, w, vpp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;

    }
    return gamma0 + K1 + K2b;
}
template <typename Q>
auto right_same_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vpp, int i_in, int spin, char channel) -> Q
{
    Q gamma0, K1, K2b;
    switch (channel){
        case 'a':
            gamma0 = vertex.spinvertex.irred.vval(i3, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.avertex.K1_vvalsmooth(i3, w, i_in, spin, vertex.spinvertex.tvertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.avertex.K2b_vvalsmooth(i3, w, vpp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;

        case 'p':
            gamma0 = vertex.spinvertex.irred.vval(i3, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.pvertex.K1_vvalsmooth(i3, w, i_in, spin);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.pvertex.K2b_vvalsmooth(i3, w, vpp, i_in, spin);
#endif
            break;
        case 't' :
            gamma0 = vertex.spinvertex.irred.vval(i3, i_in);
#if DIAG_CLASS >=1
            K1 = vertex.spinvertex.tvertex.K1_vvalsmooth(i3, w, i_in, spin, vertex.spinvertex.avertex);
#endif
#if DIAG_CLASS >=2
            K2b = vertex.spinvertex.tvertex.K2b_vvalsmooth(i3, w, vpp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;

    }
    return gamma0 + K1 + K2b;
}
template <typename Q>
auto left_diff_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in, int spin, char channel) -> Q {

    Q K2, K3, gammaRb;

    switch (channel){
        case 'a' :
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.avertex.K2_vvalsmooth(i1, w, v, i_in, spin, vertex.spinvertex.tvertex);
            gammaRb = vertex.spinvertex.gammaRb(i1, w, v, vpp, i_in, spin, 'a');
#endif
#if DIAG_CLASS >=3
            K3 = vertex.spinvertex.avertex.K3_vvalsmooth(i1, w, v, vpp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.pvertex.K2_vvalsmooth(i1, w, v, i_in, spin);
            gammaRb = vertex.spinvertex.gammaRb(i1, w, v, vpp, i_in, spin, 'p');
#endif
#if DIAG_CLASS >=3
            K3 = vertex.spinvertex.pvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, spin);
#endif
            break;
        case 't':
#if DIAG_CLASS >=2
            K2 = vertex.spinvertex.tvertex.K2_vvalsmooth(i1, w, v, i_in, spin, vertex.spinvertex.avertex);
            gammaRb = vertex.spinvertex.gammaRb(i1, w, v, vpp, i_in, spin, 't');;
#endif
#if DIAG_CLASS >=3
            K3 = vertex.spinvertex.tvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2 + K3 + gammaRb;
}
template <typename Q>
auto right_diff_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in, int spin, char channel) -> Q {

    Q K2b, K3, gammaRb;

    switch (channel){
        case 'a' :
#if DIAG_CLASS >= 2
            K2b = vertex.spinvertex.avertex.K2b_vvalsmooth(i3, w, vp, i_in, spin, vertex.spinvertex.tvertex);
            gammaRb = vertex.spinvertex.gammaRb(i3, w, vpp, vp, i_in, spin, 'a');
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.spinvertex.avertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, spin, vertex.spinvertex.tvertex);
#endif
            break;
        case 'p':
#if DIAG_CLASS >= 2
            K2b = vertex.spinvertex.pvertex.K2b_vvalsmooth(i3, w, vp, i_in, spin);
            gammaRb = vertex.spinvertex.gammaRb(i3, w, vpp, vp, i_in, spin, 'p');
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.spinvertex.pvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, spin);
#endif
            break;
        case 't':
#if DIAG_CLASS >= 2
            K2b = vertex.spinvertex.tvertex.K2b_vvalsmooth(i3, w, vp, i_in, spin, vertex.spinvertex.avertex);
            gammaRb = vertex.spinvertex.gammaRb(i3, w, vpp, vp, i_in, spin, 't');
#endif
#if DIAG_CLASS >= 3
            K3 = vertex.spinvertex.tvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, spin, vertex.spinvertex.avertex);
#endif
            break;
        default:
            return 0.;
    }
    return K2b + K3 + gammaRb;
}



#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
