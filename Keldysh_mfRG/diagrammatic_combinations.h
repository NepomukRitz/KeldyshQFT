//
// Created by Elias Walter on 2020-01-15.
//

#ifndef KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
#define KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H

#include "vertex.h"

//template <typename Q, char channel>
//auto left_same_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in) -> Q;
//template <typename Q, char channel>
//auto right_same_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vpp, int i_in) -> Q;
//template <typename Q, char channel>
//auto left_diff_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in) -> Q;
//template <typename Q, char channel>
//auto right_diff_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> Q;



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


//template <>
//auto left_same_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double vpp, int i_in) -> comp {
//    return vertex.densvertex.irred.vval(i1, i_in)
//           + vertex.densvertex.avertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.tvertex)
//           + vertex.densvertex.avertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.tvertex);
//}
//template <>
//auto left_same_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double vpp, int i_in) -> comp {
//    return vertex.densvertex.irred.vval(i1, i_in)
//           + vertex.densvertex.pvertex.K1_vvalsmooth(i1, w, i_in)
//           + vertex.densvertex.pvertex.K2b_vvalsmooth(i1, w, vpp, i_in);
//}
//template <>
//auto left_same_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double vpp, int i_in) -> comp {
//    return vertex.densvertex.irred.vval(i1, i_in)
//           + vertex.densvertex.tvertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.avertex)
//           + vertex.densvertex.tvertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.avertex);
//}
//
//
//template <>
//auto right_same_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vpp, int i_in) -> comp {
//    return vertex.densvertex.irred.vval(i3, i_in)
//           + vertex.densvertex.avertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.tvertex)
//           + vertex.densvertex.avertex.K2_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.tvertex);
//}
//template <>
//auto right_same_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vpp, int i_in) -> comp {
//    return vertex.densvertex.irred.vval(i3, i_in)
//           + vertex.densvertex.pvertex.K1_vvalsmooth(i3, w, i_in)
//           + vertex.densvertex.pvertex.K2_vvalsmooth(i3, w, vpp, i_in);
//}
//template <>
//auto right_same_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vpp, int i_in) -> comp {
//    return vertex.densvertex.irred.vval(i3, i_in)
//           + vertex.densvertex.tvertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.avertex)
//           + vertex.densvertex.tvertex.K2_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.avertex);
//}



template <typename Q>
auto left_diff_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in, char channel) -> Q {
//TODO is gammaRb in K2 and in K3 or only in K3? Or ae there contributions in both cases, but aren't the same?
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
            break;
        default:
            return 0.;
    }
    return K2 + K3 + gammaRb;
}

template <typename Q>
auto right_diff_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in, char channel) -> Q {
    //TODO is gammaRb in K2 and in K3 or only in K3? Or ae there contributions in both cases, but aren't the same?
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



//template <>
//auto left_diff_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double v, double vpp, int i_in) -> comp {
//    return vertex.densvertex.avertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.tvertex)
//           + vertex.densvertex.avertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.tvertex)
//           + vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 'a');
//
//}
//template <>
//auto left_diff_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double v, double vpp, int i_in) -> comp {
//    return vertex.densvertex.pvertex.K2_vvalsmooth(i1, w, v, i_in)
//           + vertex.densvertex.pvertex.K3_vvalsmooth(i1, w, v, vpp, i_in)
//           + vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 'p');
//
//}
//template <>
//auto left_diff_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double v, double vpp, int i_in) -> comp {
//    return vertex.densvertex.tvertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.avertex)
//           + vertex.densvertex.tvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.avertex)
//           + vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 't');
//
//}
//
//
//template <>
//auto right_diff_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> comp {
//    return vertex.densvertex.avertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.tvertex)
//           + vertex.densvertex.avertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.tvertex)
//           + vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 'a');
//
//}
//template <>
//auto right_diff_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> comp {
//    return vertex.densvertex.pvertex.K2b_vvalsmooth(i3, w, vp, i_in)
//           + vertex.densvertex.pvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in)
//           + vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 'p');
//
//}
//template <>
//auto right_diff_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> comp {
//    return vertex.densvertex.tvertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.avertex)
//           + vertex.densvertex.tvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.avertex)
//           + vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 't');
//
//}


#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
