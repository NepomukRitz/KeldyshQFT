//
// Created by Elias Walter on 2020-01-15.
//

#ifndef KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
#define KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H

#include "vertex.h"

template <typename Q, char channel>
auto left_same_bare (Vertex<fullvert<Q> >& vertex, int i1, double w, double vpp, int i_in) -> Q {}

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



resp1 = vertex1.densvertex.irred.vval(i1, i_in);
resp2 = vertex1.densvertex.avertex.K1_vvalsmooth(i1, wa, i_in, vertex1.densvertex.tvertex);
resp3 = vertex1.densvertex.avertex.K2b_vvalsmooth(i1, wa, vppa, i_in, vertex1.densvertex.tvertex);
resp4 = vertex2.densvertex.irred.vval(i3, i_in);
resp5 = vertex2.densvertex.avertex.K1_vvalsmooth(i3, wa, i_in, vertex2.densvertex.tvertex);
resp6 = vertex2.densvertex.avertex.K2_vvalsmooth(i3, wa, vppa, i_in, vertex2.densvertex.tvertex);


#endif //KELDYSH_MFRG_DIAGRAMMATIC_COMBINATIONS_H
