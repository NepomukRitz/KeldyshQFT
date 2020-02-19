/**
 * Classes / functions for computing bubbles.
 *
 * class Bubble            : Object containing two propagators, with call operator
 * class Integrand_Ki      : Object with call operator providing the integrand for bubble frequency integration
 *                           in diagrammatic class Ki:
 *                           Computing Vertex * Bubble * Vertex, and performing the internal Keldysh summation
 * class Integrand_Ki_diff : Same as Integrand_Ki, for differentiated bubble
 * bubble_function()       : Computing the bubble frequency integral, using the integrand classes and the integrator
 *                           from "integrator.h". Using MPI+OMP parallelization in the external arguments.
 */

#ifndef KELDYSH_MFRG_BUBBLES_H
#define KELDYSH_MFRG_BUBBLES_H

#include "vertex.h"
#include "propagator.h"
#include "integrator.h"
#include "selfenergy.h"
#include "util.h"
#include "mpi_setup.h"
#include "diagrammatic_combinations.h"


/**
 * Class combining two propagators, either GG or GS+SG
 */
class Bubble{
    Propagator& g;
    Propagator& s;
    bool dot;
public:

    /**
     * Constructor:
     * @param propagatorG : first propagator (always a standard one)
     * @param propagatorS : second propagator (standard or single-scale/differentiated, depending on "dot_in")
     * @param dot_in      : whether to compute standard (false) or differentiated (true) bubble
     */
    Bubble(Propagator& propagatorG, Propagator& propagatorS, bool dot_in)
        :g(propagatorG), s(propagatorS), dot(dot_in) {};

    /**
     * Call operator:
     * @param iK    : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param v1    : frequency of first propagator
     * @param v2    : frequency of second propagator
     * @return comp : value of the bubble evaluated at (iK, v1, v2)
     */
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
                default:
                    return 0.;
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
                default:
                    return 0.;
            }
        }
        return ans;
    }
};


class IntegrandBubble{
    Propagator& g1;
    Propagator& g2;
    bool diff;
    double w;
    int iK;
    char channel;

public:
    IntegrandBubble(Propagator& g1_in, Propagator& g2_in, bool diff_in, double w_in, int iK_in, char channel_in)
            : g1(g1_in), g2(g2_in), diff(diff_in), w(w_in), iK(iK_in), channel(channel_in) {};

    auto operator() (double vpp) -> comp {
        comp ans;
        double v1, v2;
        Bubble Pi(g1, g2, diff);

        switch(channel){
            case 'p':
                v1 = w/2.+vpp;
                v2 = w/2.-vpp;
                break;

            case 'a': case 't':
                v1 = vpp-w/2.;
                v2 = vpp+w/2.;
                break;
            default:
                v1 = 0.;
                v2 = 0.;
                cout << "Error in IntegrandBubble";
        }

        return Pi.value(iK, v1, v2)/(2.*pi*im_unit);
    }
};

/**
 * Integrand classes for non-differentiated bubble contributing to diagrammatic class K1, K2, K3
 */
template <typename Q> class Integrand_K1 {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    char channel;
    double w;
public:
    /**
     * Constructor:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0 or 1) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     */
    Integrand_K1(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in,    int i_in_in, char ch_in)
        :                     vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                w(w_in),  i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * Pival * res_r_V;

                    break;
                default: ;
            }
        }
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
    /**
     * Constructor:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0,...,4) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param v_in       : external fermionic frequency \nu
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     * @param pt_in      : For multi-loop calculation: specify if one computes left ('L') or right ('R')
     *                     multi-loop contribution.
     */
    Integrand_K2(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in,   Bubble& Pi_in, int i0_in, double w_in, double v_in,   int i_in_in,     char ch_in, char pt_in)
        :                     vertex1(vertex1_in),              vertex2(vertex2_in),       Pi(Pi_in),                w(w_in),     v(v_in), i_in(i_in_in), channel(ch_in), part(pt_in)
    {
        // converting index i0_in (0,...,4) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K2a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K2p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K2t[i0_in]; break;
            default: ;
        }
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) -> Q {
        if (part != 'L') return 0.;  // right part of multi-loop contribution does not contribute to K2 class
        // TODO: attention: central part does contribute, but we treat it as right part of previous loop --> fix this!!

        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Contributions: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    res_l_V = vertex1.spinvertex.gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Contributions: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    res_l_V = vertex1.spinvertex.gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1.spinvertex.gammaRb(indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    res_l_V = vertex1.spinvertex.gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1.spinvertex.gammaRb(indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 1, channel);

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
    /**
     * Constructor:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0,...,5) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param v_in       : external fermionic frequency \nu
     * @param vp_in      : external fermionic frequency \nu'
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     * @param pt_in      : For multi-loop calculation: specify if one computes left ('L') or right ('R')
     *                     multi-loop contribution.
     */
    Integrand_K3(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in, double vp_in, int i_in_in, char ch_in, char pt_in)
        :                     vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                 w(w_in),     v(v_in),    vp(vp_in), i_in(i_in_in), channel(ch_in), part(pt_in)
    {
        i0 = non_zero_Keldysh_K3[i0_in]; // converting index i0_in (0,...,5) into actual Keldysh index i0 (0,...,15)
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V,  res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                               //Contributions: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    if(part=='L'){
                        res_l_V = vertex1.spinvertex.gammaRb( indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 0, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2.spinvertex.gammaRb(indices[1], w, vp, vpp, i_in, 0, channel);
                    }
                    else ;
                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                               //Contributions: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    if(part=='L'){
                        res_l_V = vertex1.spinvertex.gammaRb( indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = vertex1.spinvertex.gammaRb( indices[0], w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 1, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2.spinvertex.gammaRb(indices[1], w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = vertex2.spinvertex.gammaRb(indices[1], w, vp, vpp, i_in, 1, channel);
                    }
                    else ;
                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                               //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    if(part=='L'){
                        res_l_V = vertex1.spinvertex.gammaRb( indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = vertex1.spinvertex.gammaRb( indices[0], w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 1, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2.spinvertex.gammaRb(indices[1], w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = vertex2.spinvertex.gammaRb(indices[1], w, vp, vpp, i_in, 1, channel);
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

/**
 * Integrand classes for differentiated bubble contributing to diagrammatic class K1, K2, K3
 */
template <typename Q> class Integrand_K1_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    char channel;
    double w;
public:
    /**
     * Constructor:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0 or 1) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     */
    Integrand_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in, int i_in_in, char ch_in)
        :                          vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),               w(w_in), i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        assert(i0_in<2);
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * Pival * res_r_V;

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
    /**
     * Constructor:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0,...,4) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param v_in       : external fermionic frequency \nu
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     */
    Integrand_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in, int i_in_in, char ch_in)
        :                          vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                w(w_in),     v(v_in), i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0,...,4) into actual Keldysh index i0 (0,...,15)
        assert(i0_in<6);
        switch (channel){
            case 'a': i0 = non_zero_Keldysh_K2a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K2p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K2t[i0_in]; break;
            default: ;
        }
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) -> Q {

        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * Pival * res_r_V;

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
    /**
     * Constructor:
     * @param vertex1_in : left vertex
     * @param vertex2_in : right vertex
     * @param Pi_in      : Bubble object connecting the left and right vertex
     * @param i0_in      : index (0,...,5) specifying the (external) Keldysh component of integrand object
     *                     (converted into actual Keldysh index i0 within the constructor)
     * @param w_in       : external bosonic frequency \omega
     * @param v_in       : external fermionic frequency \nu
     * @param vp_in      : external fermionic frequency \nu'
     * @param i_in_in    : external index for internal structure
     * @param ch_in      : diagrammatic channel ('a', 'p', 't')
     */
    Integrand_K3_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in, double vp_in, int i_in_in, char ch_in)
        :                          vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                      w(w_in),     v(v_in),    vp(vp_in), i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0,...,5) into actual Keldysh index i0 (0,...,15)
        assert(i0_in<7);
        i0 = non_zero_Keldysh_K3[i0_in];
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1.spinvertex.avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V + V^*Pi*V^
                    vertex1.spinvertex.pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 1, channel);
                    res_r_Vhat = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1.spinvertex.tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 1, channel);
                    res_r_Vhat = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * Pival * res_r_V;

                    break;
                default: ;
            }
        }
        return res;
    }
};


/**
 * Function that computes the bubble frequency integral for all external parameters,
 * including MPI and OMP parallelization. The first (reference) argument contains the result.
 * Distributing the tasks between MPI/OMP is done explicitly before each for-loop for each diagrammatic class
 * K1, K2, K3. (TODO: handle this globally / in parameters.h?)
 *
 * @tparam Q      : data type (comp or double)
 * @param dgamma  : result (full Vertex object)
 * @param vertex1 : left vertex in bubble
 * @param vertex2 : right vertex in bubble
 * @param G       : first propagator of bubble
 * @param S       : second propagator of bubble
 * @param channel : diagrammatic channel ('a', 'p', or 't')
 * @param diff    : whether or not the bubble is a differentiated one
 * @param part    : For multi-loop calculation: specify if one computes left ('L') or right ('R')
 *                  multi-loop contribution. Use '.' for first loop order.
 */
template <typename Q>
void bubble_function(Vertex<fullvert<Q> >& dgamma, Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2,
                     Propagator& G, Propagator& S, char channel, bool diff, char part)
{
    Bubble Pi(G, S, diff); // initialize bubble object

    int nw1_w = 0, nw2_w = 0, nw2_v = 0, nw3_w = 0, nw3_v = 0, nw3_vp = 0;
    Q prefactor;

    // set channel-specific frequency ranges and prefactor (1, 1/2, -1 for a, p, t)
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

    int mpi_size = mpi_world_size(); // number of mpi processes
    int mpi_rank = mpi_world_rank(); // number of the current mpi process

#ifdef DIAG_CLASS
#if DIAG_CLASS>=0
    double tK1 = get_time();
    /*K1 contributions*/
    int n_mpi = nK_K1 * nw1_w; // set external arguments for MPI-parallelization (# of tasks distributed via MPI)
    int n_omp = n_in;          // set external arguments for OMP-parallelization (# of tasks per MPI-task distributed via OMP)

    // initialize buffer into which each MPI process writes their results
    vec<Q> K1_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    // start for-loop over external arguments, using MPI and OMP
    int iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                // converting external MPI/OMP indices to physical indices (TODO: put into extra function(s)?)
                int iK1 = i_mpi * n_omp + i_omp;
                int i0 = iK1/(nw1_w*n_in);
                int iw = iK1/(n_in) - i0*nw1_w;
                int i_in = iK1 - i0*nw1_w*n_in - iw*n_in;
                double w = bfreqs[iw];
                Q value;

                // initialize the integrand object and perform frequency integration
                // (distinguishing between differentiated and non-differentiated bubble)
                if(diff){
                    Integrand_K1_diff<Q> integrand_K1(vertex1, vertex2, Pi, i0, w, i_in, channel);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K1, w_lower_f, w_upper_f);                      //Integration over a fermionic frequency
                } //TODO: prefactor -1./(2.*pi*im_unit) into Integrand classes?
                else{
                    Integrand_K1<Q> integrand_K1(vertex1, vertex2, Pi, i0, w, i_in, channel);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K1, w_lower_f, w_upper_f);                      //Integration over a fermionic frequency
                }

                K1_buffer[iterator*n_omp + i_omp] = value; // write result of integration into MPI buffer
            }
            ++iterator;
        }
    }

    // collect+combine results from different MPI processes, reorder them appropriately
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

    // initialize buffer into which each MPI process writes their results
    vec<Q> K2_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    // start for-loop over external arguments, using MPI and OMP
    iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                // converting external MPI/OMP indices to physical indices
                int iK2 = i_mpi * n_omp + i_omp;
                int i0 = iK2 /(nw2_w * nw2_v * n_in);
                int iw = iK2 /(nw2_v * n_in) - i0*nw2_w;
                int iv = iK2 / n_in - iw*nw2_v - i0*nw2_w*nw2_v;
                int i_in = iK2 - iv*n_in - iw*nw2_v*n_in - i0*nw2_w * nw2_v * n_in;
                double w = bfreqs[iw];
                double v = ffreqs[iv];
                Q value;

                // initialize the integrand object and perform frequency integration
                // (distinguishing between differentiated and non-differentiated bubble)
                if(diff){
                    Integrand_K2_diff<Q> integrand_K2 (vertex1, vertex2, Pi, i0, w, v, i_in, channel);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K2, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
                }
                else{
                    Integrand_K2<Q> integrand_K2 (vertex1, vertex2, Pi, i0, w, v, i_in, channel, part);
                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K2, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
                }

                K2_buffer[iterator*n_omp + i_omp] = value; // write result of integration into MPI buffer
            }
            ++iterator;
        }
    }

    // collect+combine results from different MPI processes, reorder them appropriately
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

    // initialize buffer into which each MPI process writes their results
    vec<Q> K3_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    // start for-loop over external arguments, using MPI and OMP
    iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                // converting external MPI/OMP indices to physical indices
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

                // initialize the integrand object and perform frequency integration
                // (distinguishing between differentiated and non-differentiated bubble)
                if (diff){
                    Integrand_K3_diff<Q> integrand_K3 (vertex1, vertex2, Pi, i0, w, v, vp, i_in, channel);

                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K3, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency                    // y
                }
                else{
                    Integrand_K3<Q> integrand_K3 (vertex1, vertex2, Pi, i0, w, v, vp,  i_in, channel, part);

                    value = prefactor*(-1./(2.*pi*im_unit))*integrator(integrand_K3, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency
                }

                K3_buffer[iterator*n_omp + i_omp] = value; // write result of integration into MPI buffer
            }
            ++iterator;
        }
    }

    // collect+combine results from different MPI processes, reorder them appropriately
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
