//
// Created by Sa.Aguirre on 9/16/19.
//

#ifndef KELDYSH_MFRG_LOOP_H
#define KELDYSH_MFRG_LOOP_H

#include "selfenergy.h"
#include "vertex.h"
#include "propagator.h"


//TODO include spinvertex contributions!!

template <typename Q, typename T >
class IntegrandR{
    Vertex<T>& vertex;
    Propagator& propagator;
    double w;

public:
    IntegrandR(Vertex<T>& vertex_in, Propagator& prop_in, double w_in)
        : vertex(vertex_in), propagator(prop_in), w(w_in) {};

    Q operator()(double wp)
    {
        return
        -( (vertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, vertex.densvertex.tvertex) +
            vertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0) +
            vertex.densvertex.tvertex.value(3, 0., wp, w, 0, vertex.densvertex.avertex) +
            vertex.densvertex.irred.vval(3) ) * propagator.pvalsmooth(0, wp) +

           (vertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, vertex.densvertex.tvertex) +
            vertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0) +
            vertex.densvertex.tvertex.value(6, 0., wp, w, 0, vertex.densvertex.avertex) +
            vertex.densvertex.irred.vval(6) ) * conj(propagator.pvalsmooth(0, wp)) +

           (vertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, vertex.densvertex.tvertex) +
            vertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0) +
            vertex.densvertex.tvertex.value(7, 0., wp, w, 0, vertex.densvertex.avertex) +
            vertex.densvertex.irred.vval(7) ) * propagator.pvalsmooth(1, wp) );
    }
};


template <typename Q, typename T >
class IntegrandK{
    Vertex<T>& vertex;
    Propagator& propagator;
    double w;

public:
    IntegrandK(Vertex<T>& vertex_in, Propagator& prop_in, double w_in)
            : vertex(vertex_in), propagator(prop_in), w(w_in) {};

    Q operator()(double wp)
    {
        return
        -( (vertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, vertex.densvertex.tvertex) +
            vertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0) +
            vertex.densvertex.tvertex.value(1, 0., wp, w, 0, vertex.densvertex.avertex) +
            vertex.densvertex.irred.vval(1) ) * propagator.pvalsmooth(0, wp) +

           (vertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, vertex.densvertex.tvertex) +
            vertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0) +
            vertex.densvertex.tvertex.value(4, 0., wp, w, 0, vertex.densvertex.avertex) +
            vertex.densvertex.irred.vval(4) ) * conj(propagator.pvalsmooth(0, wp)) +

           (vertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), 0, vertex.densvertex.tvertex) +
            vertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), 0) +
            vertex.densvertex.tvertex.value(5, 0., wp, w, 0, vertex.densvertex.avertex) +
            vertex.densvertex.irred.vval(5) ) * propagator.pvalsmooth(1, wp) );
    }
};


SelfEnergy<comp> loop(Vertex<fullvert<comp> >& fullvertex, Propagator& prop)
{
    SelfEnergy<comp> resp = SelfEnergy<comp> ();

//#pragma omp parallel for
    for (int i=0; i<nSE; ++i){
        double w = bfreqs[i];

        IntegrandR<comp, fullvert<comp> > integrandR(fullvertex, prop, w);
        IntegrandK<comp, fullvert<comp> > integrandK(fullvertex, prop, w);

        comp integratedR = integrator(integrandR, bfreqs);
        comp integratedK = integrator(integrandK, bfreqs);

        resp.setself(0, i, integratedR);
        resp.setself(1, i, integratedK);
    }

    return resp;
}


#endif //KELDYSH_MFRG_LOOP_H
