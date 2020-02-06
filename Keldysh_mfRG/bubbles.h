//
// Created by Elias Walter on 2020-01-09.
//

#ifndef KELDYSH_MFRG_BUBBLES_H
#define KELDYSH_MFRG_BUBBLES_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"
#include "util.h"
#include "mpi_setup.h"
#include "diagrammatic_combinations.h"



class Bubble{
    Propagator& g;
    Propagator& s;
    bool dot;
public:

    Bubble(Propagator& propagatorG, Propagator& propagatorS, bool dot_in)
        :g(propagatorG), s(propagatorS), dot(dot_in) {};

    auto value(int iK, double v1, double v2) -> comp{
        comp ans;
        if(dot){
            switch (iK) {
                case 3: //AA
                    ans = conj(g.pvalsmooth(0, v1)) * conj(s.pvalsmooth(0, v2)) + conj(s.pvalsmooth(0, v1)) * conj(g.pvalsmooth(0, v2));
                    break;
                case 6: //AR
                    ans = conj(g.pvalsmooth(0, v1)) * s.pvalsmooth(0, v2) + conj(s.pvalsmooth(0, v1)) * g.pvalsmooth(0, v2);
                    break;
                case 7: //AK
                    ans = conj(g.pvalsmooth(0, v1)) * s.pvalsmooth(1, v2) + conj(s.pvalsmooth(0, v1)) * g.pvalsmooth(1, v2);
                    break;
                case 9: //RA
                    ans = g.pvalsmooth(0, v1) * conj(s.pvalsmooth(0, v2)) + s.pvalsmooth(0, v1) * conj(g.pvalsmooth(0, v2));
                    break;
                case 11://KA
                    ans = g.pvalsmooth(1, v1) * conj(s.pvalsmooth(0, v2)) + s.pvalsmooth(1, v1) * conj(g.pvalsmooth(0, v2));
                    break;
                case 12://RR
                    ans = g.pvalsmooth(0, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(0, v1) * g.pvalsmooth(0, v2);
                    break;
                case 13://RK
                    ans = g.pvalsmooth(0, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(0, v1) * g.pvalsmooth(1, v2);
                    break;
                case 14://KR
                    ans = g.pvalsmooth(1, v1) * s.pvalsmooth(0, v2) + s.pvalsmooth(1, v1) *  g.pvalsmooth(0, v2);
                    break;
                case 15://KK
                    ans = g.pvalsmooth(1, v1) * s.pvalsmooth(1, v2) + s.pvalsmooth(1, v1) * g.pvalsmooth(1, v2);
                    break;
                default: ;
            }
        }
        else
        {
            switch (iK){
                case 3: //AA
                    ans = conj(g.pvalsmooth(0, v1)) * conj(g.pvalsmooth(0, v2));
                    break;
                case 6: //AR
                    ans = conj(g.pvalsmooth(0, v1)) * g.pvalsmooth(0, v2);
                    break;
                case 7: //AK
                    ans = conj(g.pvalsmooth(0, v1)) * g.pvalsmooth(1, v2);
                    break;
                case 9: //RA
                    ans = g.pvalsmooth(0, v1) * conj(g.pvalsmooth(0, v2));
                    break;
                case 11://KA
                    ans = g.pvalsmooth(1, v1) * conj(g.pvalsmooth(0, v2));
                    break;
                case 12://RR
                    ans = g.pvalsmooth(0, v1) * g.pvalsmooth(0, v2);
                    break;
                case 13://RK
                    ans = g.pvalsmooth(0, v1) * g.pvalsmooth(1, v2);
                    break;
                case 14://KR
                    ans =  g.pvalsmooth(1, v1) *  g.pvalsmooth(0, v2);
                    break;
                case 15://KK
                    ans =  g.pvalsmooth(1, v1) *  g.pvalsmooth(1, v2);
                    break;
                default: ;
            }
        }
        return ans;
    }
};



template <typename Q> class Integrand_K1 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    char channel;
    double w;
public:
    Integrand_K1(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in,    int i_in_in, char ch_in)
        :                     vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                w(w_in),  i_in(i_in_in), channel(ch_in)
    {
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
    };

    auto operator() (double vpp) -> Q {

        Q res;
        comp Pival;
        vector<int> indices(2);
        int *i1 = &indices[0], *i3 = &indices[1];
        //TODO implement!!! It may be, however, that this implementation will be trivial, since in the higher loops there are no contributions to K1
//        for (auto i2:non_zero_Keldysh_bubble) {
//            switch (channel) {
//                case 'a':
//                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);       //vppa-1/2wa, vppa+1/2wa for the a-channel
//                    tie(i1, i3) = vertex1.spinvertex.avertex.indices_sum(i0, i2);
//                    break;
//                case 'p':
//                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);       //2wp/2+vppp, wp/2-vppp for the p-channel
//                    tie(i1, i3) = vertex1.spinvertex.pvertex.indices_sum(i0, i2);
//                    break;
//                case 't':
//                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);       //vppt-1/2wt, vppt+1/2wt for the t-channel
//                    tie(i1, i3) = vertex1.spinvertex.tvertex.indices_sum(i0, i2);
//                    break;
//                default: ;
//            }
//
//
//        }
        return res;
    }

};
template <typename Q> class Integrand_K2
{
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble&Pi;
    int i0, i_in;
    char channel, part;
    double w, v;
public:
    Integrand_K2(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& Pi_in, int i0_in, double w_in, double v_in,   int i_in_in,     char ch_in, char pt_in)
        :                     vertex1(vertex1_in),              vertex2(vertex2_in),       Pi(Pi_in),                w(w_in),     v(v_in), i_in(i_in_in), channel(ch_in), part(pt_in)
    {
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K2a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K2p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K2t[i0_in]; break;
            default: ;
        }
    };

    auto operator() (double vpp) -> Q {
        if (part != 'L') return 0.;

        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        int *i1 = &indices[0], *i3 = &indices[1];
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Contributions: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    res_l_V = vertex1.spinvertex.gammaRb(*i1, w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, *i3, w,   vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Contributions: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    res_l_V = vertex1.spinvertex.gammaRb(*i1, w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, *i3, w,   vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1.spinvertex.gammaRb(*i1, w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q>(vertex2, *i3, w,   vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    res_l_V = vertex1.spinvertex.gammaRb(*i1, w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, *i3, w,   vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1.spinvertex.gammaRb(*i1, w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q>(vertex2, *i3, w,   vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
                    break;
                default: ;
            }
        }
        return res;
    }
};
template <typename Q> class Integrand_K3
{
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    char channel, part;
    double w, v, vp;
public:
    Integrand_K3(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in, double vp_in, int i_in_in, char ch_in, char pt_in)
        :                     vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                 w(w_in),     v(v_in),    vp(vp_in), i_in(i_in_in), channel(ch_in), part(pt_in)
    {
        i0 = non_zero_Keldysh_K3[i0_in];
    };

    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V,  res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        int *i1 = &indices[0], *i3 = &indices[1];
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                               //Contributions: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    if(part=='L'){
                        res_l_V = vertex1.spinvertex.gammaRb( *i1, w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, *i3, w, vp, vpp, i_in, 0, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, *i1, w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2.spinvertex.gammaRb(*i3, w, vp, vpp, i_in, 0, channel);
                    }
                    else ;
                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                               //Contributions: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    if(part=='L'){
                        res_l_V = vertex1.spinvertex.gammaRb( *i1, w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, *i3, w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = vertex1.spinvertex.gammaRb( *i1, w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = right_diff_bare<Q>(vertex2, *i3, w, vp, vpp, i_in, 1, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, *i1, w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2.spinvertex.gammaRb(*i3, w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = left_diff_bare<Q>(vertex1, *i1, w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = vertex2.spinvertex.gammaRb(*i3, w, vp, vpp, i_in, 1, channel);
                    }
                    else ;
                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                               //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    if(part=='L'){
                        res_l_V = vertex1.spinvertex.gammaRb( *i1, w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, *i3, w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = vertex1.spinvertex.gammaRb( *i1, w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = right_diff_bare<Q>(vertex2, *i3, w, vp, vpp, i_in, 1, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, *i1, w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2.spinvertex.gammaRb(*i3, w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = left_diff_bare<Q>(vertex1, *i1, w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = vertex2.spinvertex.gammaRb(*i3, w, vp, vpp, i_in, 1, channel);
                    }
                    else ;
                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
                    break;
                default: ;
            }
        }
        return res;
    }
};

template <typename Q> class Integrand_K1_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    char channel;
    double w;
public:
    Integrand_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in, int i_in_in, char ch_in)
        :                          vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),               w(w_in), i_in(i_in_in), channel(ch_in)
    {
        assert(i0_in<2);
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
    };

    //This is a call operator
    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        int *i1 = &indices[0], *i3 = &indices[1];
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_same_bare<Q> (vertex1, *i1, w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, *i3, w, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_same_bare<Q> (vertex1, *i1, w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, *i3, w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, *i1, w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, *i3, w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_same_bare<Q> (vertex1, *i1, w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, *i3, w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, *i1, w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, *i3, w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_Vhat+res_l_Vhat) * Pival * res_r_Vhat;

                    break;
                default: ;
            }
        }
        return res;
    }
};
template <typename Q> class Integrand_K2_diff {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble& Pi;
    int i0, i_in;
    char channel;
    double w, v;
public:
    Integrand_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in, int i_in_in, char ch_in)
        :                          vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                w(w_in),     v(v_in), i_in(i_in_in), channel(ch_in)
    {
        assert(i0_in<6);
        switch (channel){
            case 'a': i0 = non_zero_Keldysh_K2a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K2p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K2t[i0_in]; break;
            default: ;
        }
    };

    //This is a call operator
    auto operator() (double vpp) -> Q {

        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        int *i1 = &indices[0], *i3 = &indices[1];
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, *i1, w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, *i3, w,    vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, *i1, w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, *i3, w,    vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, *i1, w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, *i3, w,    vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, *i1, w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, *i3, w,    vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, *i1, w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, *i3, w,    vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_Vhat+res_l_Vhat) * Pival * res_r_Vhat;

                    break;
                default: ;
            }
        }
        return res;
    }
};
template <typename Q> class Integrand_K3_diff {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble&Pi;
    int i0, i_in;
    char channel;
    double w, v, vp;
public:
    Integrand_K3_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in, double vp_in, int i_in_in, char ch_in)
        :                          vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                      w(w_in),     v(v_in),    vp(vp_in), i_in(i_in_in), channel(ch_in)
    {
        assert(i0_in<7);
        i0 = non_zero_Keldysh_K3[i0_in];
    };

    //This is a call operator
    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        int *i1 = &indices[0], *i3 = &indices[1];
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, *i1, w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, *i3, w, vp, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, *i1, w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, *i3, w, vp, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, *i1, w, v,  vpp, i_in, 1, channel);
                    res_r_Vhat = right_diff_bare<Q> (vertex2, *i3, w, vp, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, *i1, w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, *i3, w, vp, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, *i1, w, v,  vpp, i_in, 1, channel);
                    res_r_Vhat = right_diff_bare<Q> (vertex2, *i3, w, vp, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_Vhat+res_l_Vhat) * Pival * res_r_Vhat;

                    break;
                default: ;
            }
        }
        return res;
    }
};


template <typename Q>
void bubble_function(Vertex<fullvert<Q> >& dgamma, Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S, char channel, bool diff, char part)
{
    Bubble Pi(G,S, diff);

    int nw1_w = 0, nw2_w = 0, nw2_v = 0, nw3_w = 0, nw3_v = 0, nw3_vp = 0;
    Q prefactor;

    switch (channel) {
        case 'a':
            nw1_w = nw1_wa;
            nw2_w = nw2_wa;
            nw2_v = nw2_nua;
            nw3_w = nw3_wa;
            nw3_v = nw3_nua;
            nw3_vp = nw3_nuap;
            prefactor = 1.;
            break;
        case 'p':
            nw1_w = nw1_wp;
            nw2_w = nw2_wp;
            nw2_v = nw2_nup;
            nw3_w = nw3_wp;
            nw3_v = nw3_nup;
            nw3_vp = nw3_nupp;
            prefactor = 0.5;
            break;
        case 't':
            nw1_w = nw1_wt;
            nw2_w = nw2_wt;
            nw2_v = nw2_nut;
            nw3_w = nw3_wt;
            nw3_v = nw3_nut;
            nw3_vp = nw3_nutp;
            prefactor = -1.;
            break;
        default: ;
    }

    int mpi_size = mpi_world_size();
    int mpi_rank = mpi_world_rank();

#ifdef DIAG_CLASS
#if DIAG_CLASS>=1
    double tK1 = get_time();
    /*K1 contributions*/
    int n_mpi = nK_K1 * nw1_w;
    int n_omp = n_in;

    vec<Q> K1_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    int iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                int iK1 = i_mpi * n_omp + i_omp;
                int i0 = iK1/(nw1_w*n_in);
                int iw = iK1/(n_in) - i0*nw1_w;
                int i_in = iK1 - i0*nw1_w*n_in - iw*n_in;
                double w = bfreqs[iw];
                Q value;

                if(diff){
                    Integrand_K1_diff<Q> integrand_K1(vertex1, vertex2, Pi, i0, w, i_in, channel);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K1, w_lower_f, w_upper_f);                      //Integration over a fermionic frequency
                }
                else{
                    Integrand_K1<Q> integrand_K1(vertex1, vertex2, Pi, i0, w, i_in, channel);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K1, w_lower_f, w_upper_f);                      //Integration over a fermionic frequency
                }

                K1_buffer[iterator*n_omp + i_omp] = value;
            }
            ++iterator;
        }
    }

    vec<Q> K1_result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(K1_buffer, K1_result, n_mpi, n_omp);
    vec<Q> K1_ordered_result = mpi_reorder_result(K1_result, n_mpi, n_omp);

    switch (channel) {
        case 'a': dgamma.spinvertex.avertex.K1 += K1_ordered_result; break;
        case 'p': dgamma.spinvertex.pvertex.K1 += K1_ordered_result; break;
        case 't': dgamma.spinvertex.tvertex.K1 += K1_ordered_result; break;
        default: ;
    }

    // dgamma.spinvertex.pvertex.K1_addvert(i0, iwp, i_in, value); // old version w/o mpi

    print("K1");  print(channel); print(" done: ");
    get_time(tK1);
#endif

#if DIAG_CLASS>=2
    double tK2 = get_time();
    /*K2 contributions*/
    n_mpi = nK_K2 * nw2_w;
    n_omp = nw2_v * n_in;

    vec<Q> K2_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                int iK2 = i_mpi * n_omp + i_omp;
                int i0 = iK2 /(nw2_w * nw2_v * n_in);
                int iw = iK2 /(nw2_v * n_in) - i0*nw2_w;
                int iv = iK2 / n_in - iw*nw2_v - i0*nw2_w*nw2_v;
                int i_in = iK2 - iv*n_in - iw*nw2_v*n_in - i0*nw2_w * nw2_v * n_in;
                double w = bfreqs[iw];
                double v = ffreqs[iv];
                Q value;

                if(diff){
                    Integrand_K2_diff<Q> integrand_K2 (vertex1, vertex2, Pi, i0, w, v, i_in, channel);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K2, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
                }
                else{
                    Integrand_K2<Q> integrand_K2 (vertex1, vertex2, Pi, i0, w, v, i_in, channel, part);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K2, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
                }

                K2_buffer[iterator*n_omp + i_omp] = value;
            }
            ++iterator;
        }
    }

    vec<Q> K2_result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(K2_buffer, K2_result, n_mpi, n_omp);
    vec<Q> K2_ordered_result = mpi_reorder_result(K2_result, n_mpi, n_omp);

    switch (channel) {
        case 'a': dgamma.spinvertex.avertex.K2 += K2_ordered_result; break;
        case 'p': dgamma.spinvertex.pvertex.K2 += K2_ordered_result; break;
        case 't': dgamma.spinvertex.tvertex.K2 += K2_ordered_result; break;
        default: ;
    }

    print("K2"); print(channel); print(" done: ");
    get_time(tK2);
#endif

#if DIAG_CLASS>=3
    double tK3 = get_time();
    /*K3 contributions*/
    n_mpi = nK_K3 * nw3_w * nw3_v;
    n_omp = nw3_vp * n_in;

    vec<Q> K3_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                int iK3 = i_mpi * n_omp + i_omp;
                int i0 = iK3/(nw3_w * nw3_v * nw3_vp * n_in);
                int iw = iK3/(nw3_v * nw3_vp * n_in) - i0*nw3_w;
                int iv = iK3/(nw3_v * n_in) - i0*nw3_w*nw3_v - iw*nw3_v;
                int ivp =iK3/(n_in) - i0*nw3_w*nw3_v*nw3_vp - iw*nw3_v*nw3_vp - iv*nw3_vp;
                int i_in = iK3 - i0*nw3_w*nw3_v*nw3_vp*n_in - iw*nw3_v*nw3_vp*n_in - iv*nw3_vp*n_in - ivp*n_in;
                double w = bfreqs[iw];
                double v = ffreqs[iv];
                double vp = ffreqs[ivp];
                Q value;


                if (diff){
                    Integrand_K3_diff<Q> integrand_K3 (vertex1, vertex2, Pi, i0, w, v, vp, i_in, channel);

                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K3, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency                    // y
                }
                else{
                    Integrand_K3<Q> integrand_K3 (vertex1, vertex2, Pi, i0, w, v, vp,  i_in, channel, part);

                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K3, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
                }

                K3_buffer[iterator*n_omp + i_omp] = value;
            }
            ++iterator;
        }
    }

    vec<Q> K3_result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(K3_buffer, K3_result, n_mpi, n_omp);
    vec<Q> K3_ordered_result = mpi_reorder_result(K3_result, n_mpi, n_omp);

    switch (channel) {
        case 'a': dgamma.spinvertex.avertex.K3 += K3_ordered_result; break;
        case 'p': dgamma.spinvertex.pvertex.K3 += K3_ordered_result; break;
        case 't': dgamma.spinvertex.tvertex.K3 += K3_ordered_result; break;
        default: ;
    }

    // dgamma.spinvertex.pvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, value); // old version w/o mpi

    print("K3"); print(channel); print(" done: ");
    get_time(tK3);
#endif

#if DIAG_CLASS>=4
    print("Damn son, this is a bad error");
#endif
#endif
}

#endif //KELDYSH_MFRG_BUBBLES_H
