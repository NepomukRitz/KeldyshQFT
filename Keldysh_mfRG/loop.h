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
    int i_in;

public:
    IntegrandR(Vertex<T>& vertex_in, Propagator& prop_in, double w_in, int i_in_in)
        : vertex(vertex_in), propagator(prop_in), w(w_in), i_in(i_in_in) {};

    auto operator()(double wp) -> Q
    {
        Q resp1 = vertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) ;   //Result should always be real
        Q resp2 = vertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) ;                              //Should always be 0
        Q resp3 = vertex.densvertex.tvertex.value(3, 0., wp, w, i_in, vertex.densvertex.avertex) ;                      //Equals resp7
        Q resp4 = vertex.densvertex.irred.vval(3, i_in);                                                                      //Always 0.

        Q resp5 = vertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);    //Should always be 0.
        Q resp6 = vertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                               //Result should always be real
        Q resp7 = vertex.densvertex.tvertex.value(6, 0., wp, w, i_in, vertex.densvertex.avertex);                       //Equals resp3
        Q resp8 = vertex.densvertex.irred.vval(6, i_in);                                                                      //Always 0.

        Q resp9  = vertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);   //Something
        Q resp10 = vertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                              //Something
        Q resp11 = vertex.densvertex.tvertex.value(7, 0., wp, w, i_in, vertex.densvertex.avertex);                      //Something
        Q resp12 = vertex.densvertex.irred.vval(7, i_in);  // =0 ?                                                            //Always 0.5*U

//        Q resp1 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(3, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(3, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp2 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(3, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(3, wp+w, i_in, 1);
//        Q resp3 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(3, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(3, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp4 = vertex.spinvertex.irred.vval(3, i_in);
//
//        Q resp5 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(6, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(6, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp6 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(6, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(6, wp+w, i_in, 1);
//        Q resp7 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(6, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(6, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp8 = vertex.spinvertex.irred.vval(6, i_in);
//
//        Q resp9 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(7, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(7, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp10= 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(7, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(7, wp+w, i_in, 1);
//        Q resp11= 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(7, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(7, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp12= vertex.spinvertex.irred.vval(7, i_in);

        Q aid1 = propagator.pvalsmooth(0, wp);
        Q aid2 = conj(propagator.pvalsmooth(0, wp));
        Q aid3 = propagator.pvalsmooth(1, wp);

        Q ans = ((resp1+resp2+resp3+resp4)*aid1 + (resp5+resp6+resp7+resp8)*aid2 + (resp9+resp10+resp11+resp12)*aid3);

        return ans;


//        return
//        -( (vertex.densvertex.avertex.value(3, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(3, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(3, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(3) ) * propagator.pvalsmooth(0, wp) +
//
//           (vertex.densvertex.avertex.value(6, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(6, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(6, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(6) ) * conj(propagator.pvalsmooth(0, wp)) +
//
//           (vertex.densvertex.avertex.value(7, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(7, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(7, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(7) ) * propagator.pvalsmooth(1, wp) );
    }
};


template <typename Q, typename T >
class IntegrandK{
    Vertex<T>& vertex;
    Propagator& propagator;
    double w;
    int i_in;

public:
    IntegrandK(Vertex<T>& vertex_in, Propagator& prop_in, double w_in, int i_in_in)
            : vertex(vertex_in), propagator(prop_in), w(w_in), i_in(i_in_in) {};

    auto operator()(double wp) -> Q
    {
        Q resp1 = vertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) ;
        Q resp2 = vertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) ;
        Q resp3 = vertex.densvertex.tvertex.value(1, 0., wp, w, i_in, vertex.densvertex.avertex) ;
        Q resp4 = vertex.densvertex.irred.vval(1, i_in);      //  =0.?                                                        //Always 0.5*U

        Q resp5 = vertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);    //"Similar" to resp1 (i.e. T3(resp1))
        Q resp6 = vertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                               //"Similar" to resp2 (i.e. TC(resp2))
        Q resp7 = vertex.densvertex.tvertex.value(4, 0., wp, w, i_in, vertex.densvertex.avertex);                       //Equal to resp3
        Q resp8 = vertex.densvertex.irred.vval(4, i_in);      // =0.?                                                         //Always 0.5*U

        Q resp9  = vertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex);   //Note this value is real
        Q resp10 = vertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in);                              //Note this value is real
        Q resp11 = vertex.densvertex.tvertex.value(5, 0., wp, w, i_in, vertex.densvertex.avertex);                      //Should always be 0.
        Q resp12 = vertex.densvertex.irred.vval(5, i_in);                                                                      //Always 0.

        Q aid1 = propagator.pvalsmooth(0, wp);
        Q aid2 = conj(propagator.pvalsmooth(0, wp));
        Q aid3 = propagator.pvalsmooth(1, wp);

        Q ans= ((resp1+resp2+resp3+resp4)*aid1 + (resp5+resp6+resp7+resp8)*aid2 + (resp9+resp10+resp11+resp12)*aid3);

//        Q resp1 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(1, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(1, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp2 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(1, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(1, wp+w, i_in, 1);
//        Q resp3 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(1, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(1, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp4 = vertex.spinvertex.irred.vval(1, i_in);
//
//        Q resp5 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(4, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(4, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp6 = 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(4, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(4, wp+w, i_in, 1);
//        Q resp7 = 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(4, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(4, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp8 = vertex.spinvertex.irred.vval(4, i_in);
//
//        Q resp9 = 2.*vertex.spinvertex.avertex.K1_vvalsmooth(5, wp-w, i_in, 0, vertex.spinvertex.tvertex) + vertex.spinvertex.avertex.K1_vvalsmooth(5, wp-w, i_in, 1, vertex.spinvertex.tvertex);
//        Q resp10= 2.*vertex.spinvertex.pvertex.K1_vvalsmooth(5, wp+w, i_in, 0) + vertex.spinvertex.pvertex.K1_vvalsmooth(5, wp+w, i_in, 1);
//        Q resp11= 2.*vertex.spinvertex.tvertex.K1_vvalsmooth(5, 0.  , i_in, 0, vertex.spinvertex.avertex) + vertex.spinvertex.tvertex.K1_vvalsmooth(5, 0.  , i_in, 1, vertex.spinvertex.avertex);
//        Q resp12= vertex.spinvertex.irred.vval(5, i_in);
//
//        Q aid1 = propagator.pvalsmooth(0, wp);
//        Q aid2 = conj(propagator.pvalsmooth(0, wp));
//        Q aid3 = propagator.pvalsmooth(1, wp);
//
//        Q ans= ((resp1+resp2+resp3+resp4)*aid1 + (resp5+resp6+resp7+resp8)*aid2 + (resp9+resp10+resp11+resp12)*aid3);

        return ans;

//        return
//        -( (vertex.densvertex.avertex.value(1, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(1, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(1, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(1) ) * propagator.pvalsmooth(0, wp) +
//
//           (vertex.densvertex.avertex.value(4, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(4, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(4, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(4) ) * conj(propagator.pvalsmooth(0, wp)) +
//
//           (vertex.densvertex.avertex.value(5, wp-w, 0.5*(w+wp), 0.5*(w+wp), i_in, vertex.densvertex.tvertex) +
//            vertex.densvertex.pvertex.value(5, wp+w, 0.5*(w-wp), 0.5*(w-wp), i_in) +
//            vertex.densvertex.tvertex.value(5, 0., wp, w, i_in, vertex.densvertex.avertex) +
//            vertex.densvertex.irred.vval(5) ) * propagator.pvalsmooth(1, wp) );
    }
};

template <typename Q>
auto loop(Vertex<fullvert<Q> >& fullvertex, Propagator& prop) -> SelfEnergy<comp>
{
    SelfEnergy<comp> resp = SelfEnergy<comp> ();
//TODO complete support for internal structure
#pragma omp parallel for
    for (int iSE=0; iSE<nSE*n_in; ++iSE){
        int i = ( (iSE%(nSE*n_in))/n_in);
        int i_in = iSE%n_in;

        double w = ffreqs[i];

        IntegrandR<Q, fullvert<Q> > integrandR(fullvertex, prop, w, i_in);
        IntegrandK<Q, fullvert<Q> > integrandK(fullvertex, prop, w, i_in);

        comp integratedR = (1./((comp)1.i))*integrator(integrandR, w_lower_f, w_upper_f);
        comp integratedK = (1./((comp)1.i))*integrator(integrandK, w_lower_f, w_upper_f);

        resp.setself(0, i, integratedR);
        resp.setself(1, i, integratedK);
    }

    return resp;
}


#endif //KELDYSH_MFRG_LOOP_H
