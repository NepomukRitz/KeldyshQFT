//
// Created by Elias Walter on 2020-01-15.
//

#ifndef KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
#define KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H

#include "vertex.h"

template <typename Q, char channel>
auto left_same_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in) -> Q;

template <typename Q, char channel>
auto right_same_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vpp, int i_in) -> Q;

template <typename Q, char channel>
auto left_diff_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double v, double vpp, int i_in) -> Q;

template <typename Q, char channel>
auto right_diff_bare (Vertex<fullvert<Q> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> Q;




template <>
auto left_same_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double vpp, int i_in) -> comp {
    return vertex.densvertex.irred.vval(i1, i_in)
           + vertex.densvertex.avertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.tvertex)
           + vertex.densvertex.avertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.tvertex);
}
template <>
auto left_same_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double vpp, int i_in) -> comp {
    return vertex.densvertex.irred.vval(i1, i_in)
           + vertex.densvertex.pvertex.K1_vvalsmooth(i1, w, i_in)
           + vertex.densvertex.pvertex.K2b_vvalsmooth(i1, w, vpp, i_in);
}
template <>
auto left_same_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double vpp, int i_in) -> comp {
    return vertex.densvertex.irred.vval(i1, i_in)
           + vertex.densvertex.tvertex.K1_vvalsmooth(i1, w, i_in, vertex.densvertex.avertex)
           + vertex.densvertex.tvertex.K2b_vvalsmooth(i1, w, vpp, i_in, vertex.densvertex.avertex);
}


template <>
auto right_same_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vpp, int i_in) -> comp {
    return vertex.densvertex.irred.vval(i3, i_in)
           + vertex.densvertex.avertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.tvertex)
           + vertex.densvertex.avertex.K2_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.tvertex);
}
template <>
auto right_same_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vpp, int i_in) -> comp {
    return vertex.densvertex.irred.vval(i3, i_in)
           + vertex.densvertex.pvertex.K1_vvalsmooth(i3, w, i_in)
           + vertex.densvertex.pvertex.K2_vvalsmooth(i3, w, vpp, i_in);
}
template <>
auto right_same_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vpp, int i_in) -> comp {
    return vertex.densvertex.irred.vval(i3, i_in)
           + vertex.densvertex.tvertex.K1_vvalsmooth(i3, w, i_in, vertex.densvertex.avertex)
           + vertex.densvertex.tvertex.K2_vvalsmooth(i3, w, vpp, i_in, vertex.densvertex.avertex);
}



template <>
auto left_diff_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double v, double vpp, int i_in) -> comp {
    return vertex.densvertex.avertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.tvertex)
           + vertex.densvertex.avertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.tvertex)
           + vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 'a');

}
template <>
auto left_diff_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double v, double vpp, int i_in) -> comp {
    return vertex.densvertex.pvertex.K2_vvalsmooth(i1, w, v, i_in)
           + vertex.densvertex.pvertex.K3_vvalsmooth(i1, w, v, vpp, i_in)
           + vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 'p');

}
template <>
auto left_diff_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i1, double w, double v, double vpp, int i_in) -> comp {
    return vertex.densvertex.tvertex.K2_vvalsmooth(i1, w, v, i_in, vertex.densvertex.avertex)
           + vertex.densvertex.tvertex.K3_vvalsmooth(i1, w, v, vpp, i_in, vertex.densvertex.avertex)
           + vertex.densvertex.gammaRb(i1, w, v, vpp, i_in, 't');

}


template <>
auto right_diff_bare<comp, 'a'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> comp {
    return vertex.densvertex.avertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.tvertex)
           + vertex.densvertex.avertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.tvertex)
           + vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 'a');

}
template <>
auto right_diff_bare<comp, 'p'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> comp {
    return vertex.densvertex.pvertex.K2b_vvalsmooth(i3, w, vp, i_in)
           + vertex.densvertex.pvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in)
           + vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 'p');

}
template <>
auto right_diff_bare<comp, 't'> (Vertex<fullvert<comp> >& vertex, int i3, double w, double vp, double vpp, int i_in) -> comp {
    return vertex.densvertex.tvertex.K2b_vvalsmooth(i3, w, vp, i_in, vertex.densvertex.avertex)
           + vertex.densvertex.tvertex.K3_vvalsmooth(i3, w, vpp, vp, i_in, vertex.densvertex.avertex)
           + vertex.densvertex.gammaRb(i3, w, vpp, vp, i_in, 't');

}


#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
