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

#include <cmath>                        // for using the macro M_PI as pi

#include "vertex.h"                     // vertex class
#include "selfenergy.h"                 // self-energy class
#include "propagator.h"                 // propagator class
#include "integrator.h"                 // integration routines
#include "util.h"                       // measuring time, printing text output
#include "mpi_setup.h"                  // mpi parallelization routines
#include "diagrammatic_combinations.h"  // combinations of diagrammatic classes that go into the left/right vertex
                                        // in the bubble
#include "correctionFunctions.h"        // correction terms due to finite integration range

/// Class combining two propagators, either GG or GS+SG
class Bubble{
    const Propagator& g;
    const Propagator& s;
    const bool dot;
public:

    /**
     * Constructor:
     * @param propagatorG : first propagator (always a standard one)
     * @param propagatorS : second propagator (standard or single-scale/differentiated, depending on "dot_in")
     * @param dot_in      : whether to compute standard (false) or differentiated (true) bubble
     */
    Bubble(const Propagator& propagatorG, const Propagator& propagatorS, const bool dot_in)
        :g(propagatorG), s(propagatorS), dot(dot_in) {};

    /**
     * Call operator:
     * @param iK    : Keldysh index of combined bubble object (0 <= iK <= 15)
     * @param v1    : frequency of first propagator
     * @param v2    : frequency of second propagator
     * @return comp : value of the bubble evaluated at (iK, v1, v2)
     */
    auto value(int iK, double v1, double v2, int i_in) const -> comp{
        comp ans;
        if(dot){
            switch (iK) {
                case 3: //AA
                    ans = conj(g.valsmooth(0, v1, i_in)) * conj(s.valsmooth(0, v2, i_in)) + conj(s.valsmooth(0, v1, i_in)) * conj(g.valsmooth(0, v2, i_in));
                    break;
                case 6: //AR
                    ans = conj(g.valsmooth(0, v1, i_in)) * s.valsmooth(0, v2, i_in) + conj(s.valsmooth(0, v1, i_in)) * g.valsmooth(0, v2, i_in);
                    break;
                case 7: //AK
                    ans = conj(g.valsmooth(0, v1, i_in)) * s.valsmooth(1, v2, i_in) + conj(s.valsmooth(0, v1, i_in)) * g.valsmooth(1, v2, i_in);
                    break;
                case 9: //RA
                    ans = g.valsmooth(0, v1, i_in) * conj(s.valsmooth(0, v2, i_in)) + s.valsmooth(0, v1, i_in) * conj(g.valsmooth(0, v2, i_in));
                    break;
                case 11://KA
                    ans = g.valsmooth(1, v1, i_in) * conj(s.valsmooth(0, v2, i_in)) + s.valsmooth(1, v1, i_in) * conj(g.valsmooth(0, v2, i_in));
                    break;
                case 12://RR
                    ans = g.valsmooth(0, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
                    break;
                case 13://RK
                    ans = g.valsmooth(0, v1, i_in) * s.valsmooth(1, v2, i_in) + s.valsmooth(0, v1, i_in) * g.valsmooth(1, v2, i_in);
                    break;
                case 14://KR
                    ans = g.valsmooth(1, v1, i_in) * s.valsmooth(0, v2, i_in) + s.valsmooth(1, v1, i_in) *  g.valsmooth(0, v2, i_in);
                    break;
                case 15://KK
                    ans = g.valsmooth(1, v1, i_in) * s.valsmooth(1, v2, i_in) + s.valsmooth(1, v1, i_in) * g.valsmooth(1, v2, i_in);
                    break;
                default:
                    return 0.;
            }
        }
        else {
            switch (iK){ // labelling propagators from top (t: left) to bottom (t: right); a,t: G(v+w/2)G(v-w/2), p: G(w/2-v)G(w/2+v)
                case 3: //AA
                    ans = conj(g.valsmooth(0, v1, i_in)) * conj(g.valsmooth(0, v2, i_in));
                    break;
                case 6: //AR
                    ans = conj(g.valsmooth(0, v1, i_in)) * g.valsmooth(0, v2, i_in);
                    break;
                case 7: //AK
                    ans = conj(g.valsmooth(0, v1, i_in)) * g.valsmooth(1, v2, i_in);
                    break;
                case 9: //RA
                    ans = g.valsmooth(0, v1, i_in) * conj(g.valsmooth(0, v2, i_in));
                    break;
                case 11://KA
                    ans = g.valsmooth(1, v1, i_in) * conj(g.valsmooth(0, v2, i_in));
                    break;
                case 12://RR
                    ans = g.valsmooth(0, v1, i_in) * g.valsmooth(0, v2, i_in);
                    break;
                case 13://RK
                    ans = g.valsmooth(0, v1, i_in) * g.valsmooth(1, v2, i_in);
                    break;
                case 14://KR
                    ans =  g.valsmooth(1, v1, i_in) *  g.valsmooth(0, v2, i_in);
                    break;
                case 15://KK
                    ans =  g.valsmooth(1, v1, i_in) *  g.valsmooth(1, v2, i_in);
                    break;
                default:
                    return 0.;
            }
        }
        return ans;
    }
};

//Class created for debugging of the Bubbles
class IntegrandBubble{
    const Propagator& g1;
    const Propagator& g2;
    bool diff;
    double w;
    int iK;
    char channel;

public:
    /**
     * Constructor for the IntegrandBubble
     * @param g1_in     : Propagator object for the lower (right) leg
     * @param g2_in     : Propagator object for the upper (left) leg
     * @param diff_in   : Boolean defining whether bubble is differentiated or not
     * @param w_in      : Transfer frequency at which the bubble should be integrated
     * @param iK_in     : Keldysh index to be taken
     * @param channel_in: Char indicating the channel in which the bubble should be calculated and which determines the frequency transformations
     */
    IntegrandBubble(const Propagator& g1_in, const Propagator& g2_in, bool diff_in, double w_in, int iK_in, char channel_in)
            : g1(g1_in), g2(g2_in), diff(diff_in), w(w_in), iK(iK_in), channel(channel_in) {};

    /**
     * Call operator
     * @param vpp : v'', the frequency over which is being integrated
     * @return The value g1(v1)*g2(v2), where v1 and v2 are calculated according to the channel. The components of the
     * propagators taken depend on the Keldysh component
     */
    auto operator() (double vpp) const -> comp {
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
        //Make reference to the Bubble object of the actual code, making this into a useful test of code correctnes and compliance
        return Pi.value(iK, v1, v2, 0)/(2.*M_PI*glb_i);
    }
};

/// Integrand classes for non-differentiated bubble contributing to diagrammatic class K1, K2, K3
template <typename Q> class Integrand_K1 {
    const Vertex<Q>& vertex1;
    const Vertex<Q>& vertex2;
    const Bubble& Pi;
    int i0;
    const int i_in;
    const char channel;
    const double w;
#if DIAG_CLASS <= 1
    Q res_l_V[16], res_r_V[16], res_l_Vhat[16], res_r_Vhat[16];
#endif
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
    Integrand_K1(const Vertex<Q>& vertex1_in, const Vertex<Q>& vertex2_in, const Bubble& Pi_in, int i0_in, const double w_in,
                 const int i_in_in, const char ch_in)
            :            vertex1(vertex1_in),         vertex2(vertex2_in),           Pi(Pi_in),                w(w_in),
                    i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
#if DIAG_CLASS <= 1
        // For K1 class, left and right vertices do not depend on integration frequency -> precompute them to save time
        vector<int> indices(2);
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    res_l_V[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 0, channel);
                    res_r_V[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 0, channel);
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    res_l_V[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 0, channel);
                    res_r_V[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 0, channel);
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    res_l_V[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 0, channel);
                    res_r_V[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 0, channel);

                    res_l_Vhat[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 1, channel);
                    res_r_Vhat[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 1, channel);
                    break;
                default: ;
            }
        }
#endif
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q {
        Q Pival;
        Q res;
#if DIAG_CLASS >= 2
        Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        vector<int> indices(2);
#endif

        //Iterates over all Keldysh components of the bubble which are nonzero
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Flow eq: V*Pi*V
#if DIAG_CLASS <= 1 // only nonzero combinations of \int dvpp Gamma_0 Pi(vpp) Gamma_0
                    if(i0==1 && (i2 != 11 && i2 != 13)) continue;
                    if(i0==3 && (i2 != 6 && i2 != 7 && i2 != 9 && i2 != 11 && i2 != 13 && i2 != 14 && i2 != 15)) continue;
#endif
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel
#if DIAG_CLASS <= 1
                    res += res_l_V[i2] * Pival * res_r_V[i2];
#else
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
#endif
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V //+ V^*Pi*V^
#if DIAG_CLASS <= 1 // only nonzero combinations of \int dvpp Gamma_0 Pi(vpp) Gamma_0
                    if(i0==1 && (i2 != 7 && i2 != 11)) continue;
                    if(i0==5 && (i2 != 3 && i2 != 7 && i2 != 11 && i2 != 12 && i2 != 13 && i2 != 14 && i2 != 15)) continue;
#endif
                    Pival = Pi.value(i2, w/2. + vpp, w/2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
#if DIAG_CLASS <= 1
                    res += res_l_V[i2] * Pival * res_r_V[i2];
#else
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    //This is commented out on the ground of p-channel contributions being cross-symmetric
                    //Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                    //the bubble_function(...) must be changed.*/
//                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
//                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;
#endif
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
#if DIAG_CLASS <= 1 // only nonzero combinations of \int dvpp Gamma_0 Pi(vpp) Gamma_0
                    if(i0==1 && (i2 != 11 && i2 != 13)) continue;
                    if(i0==3 && (i2 != 6 && i2 != 7 && i2 != 9 && i2 != 11 && i2 != 13 && i2 != 14 && i2 != 15)) continue;
#endif
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
#if DIAG_CLASS <= 1
                    res += res_l_V[i2] * Pival * (res_r_V[i2]+res_r_Vhat[i2])
                            + (res_l_V[i2]+res_l_Vhat[i2]) * Pival * res_r_V[i2];
#else
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * Pival * res_r_V;
#endif
                    break;
                default: ;
            }
        }
        return res;
    }

};
template <typename Q> class Integrand_K2
{
    const Vertex<Q>& vertex1;
    const Vertex<Q>& vertex2;
    const Bubble&Pi;
    int i0;
    const int i_in;
    const char channel, part;
    const double w, v;
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
    Integrand_K2(const Vertex<Q>& vertex1_in, const Vertex<Q>& vertex2_in, const  Bubble& Pi_in, int i0_in,
            const double w_in, double v_in,   const int i_in_in,    const char ch_in, const char pt_in)
              :          vertex1(vertex1_in),         vertex2(vertex2_in),            Pi(Pi_in),                w(w_in),
              v(v_in), i_in(i_in_in), channel(ch_in), part(pt_in)
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
    auto operator() (double vpp) const -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        //Iterates over all Keldysh components of the bubble which are nonzero
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Contributions: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Contributions: V*Pi*V// + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, w/2. + vpp, w/2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 0, channel);

                    //This is commented out on the ground of p-channel contributions being cross-symmetric
                    //Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                    //the bubble_function(...) must be changed.*/
//                    res_l_Vhat = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 1, channel);
//                    res_r_Vhat = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w,   vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 1, channel);
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
    const Vertex<Q>& vertex1;
    const Vertex<Q>& vertex2;
    const Bubble& Pi;
    int i0;
    const int i_in;
    const char channel, part;
    const double w, v, vp;
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
    Integrand_K3(const Vertex<Q>& vertex1_in, const Vertex<Q>& vertex2_in, const Bubble& Pi_in, int i0_in,
            const double w_in, const double v_in, const double vp_in, const int i_in_in, const char ch_in, const char pt_in)
            :            vertex1(vertex1_in),         vertex2(vertex2_in),           Pi(Pi_in),                 w(w_in),
            v(v_in),    vp(vp_in), i_in(i_in_in), channel(ch_in), part(pt_in)
    {
        i0 = non_zero_Keldysh_K3[i0_in]; // converting index i0_in (0,...,5) into actual Keldysh index i0 (0,...,15)
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q {
        Q res, res_l_V, res_r_V,  res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        //Iterates over all Keldysh components of the bubble which are nonzero
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                               //Contributions: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel
                    if(part=='L'){
                        res_l_V = vertex1[0].gammaRb( indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 0, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2[0].gammaRb(indices[1], w, vp, vpp, i_in, 0, channel);
                    }
                    else ;
                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                               //Contributions: V*Pi*V; + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, w/2. + vpp, w/2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
                    if(part=='L'){
                        res_l_V = vertex1[0].gammaRb( indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 0, channel);
                        /*This is commented out on the ground of p-channel contributions being cross-symmetric
                         *Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                         * the bubble_function(...) must be changed.*/
//                        res_l_Vhat = vertex1[0].gammaRb( indices[0], w,  v, vpp, i_in, 1, channel);
//                        res_r_Vhat = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 1, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2[0].gammaRb(indices[1], w, vp, vpp, i_in, 0, channel);

//                        res_l_Vhat = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 1, channel);
//                        res_r_Vhat = vertex2[0].gammaRb(indices[1], w, vp, vpp, i_in, 1, channel);
                    }
                    else ;
                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                               //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
                    if(part=='L'){
                        res_l_V = vertex1[0].gammaRb( indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = vertex1[0].gammaRb( indices[0], w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = right_diff_bare<Q>(vertex2, indices[1], w, vp, vpp, i_in, 1, channel);
                    }
                    else if(part=='R'){
                        res_l_V = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 0, channel);
                        res_r_V = vertex2[0].gammaRb(indices[1], w, vp, vpp, i_in, 0, channel);

                        res_l_Vhat = left_diff_bare<Q>(vertex1, indices[0], w,  v, vpp, i_in, 1, channel);
                        res_r_Vhat = vertex2[0].gammaRb(indices[1], w, vp, vpp, i_in, 1, channel);
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

/// Integrand classes for differentiated bubble contributing to diagrammatic class K1, K2, K3
template <typename Q> class Integrand_K1_diff {
    const Vertex<Q>& vertex1;
    const Vertex<Q>& vertex2;
    const Bubble& Pi;
    int i0;
    const int i_in;
    const char channel;
    const double w;
#if DIAG_CLASS <= 1
    Q res_l_V[16], res_r_V[16], res_l_Vhat[16], res_r_Vhat[16];
#endif
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
    Integrand_K1_diff(const Vertex<Q>& vertex1_in, const Vertex<Q>& vertex2_in, const Bubble& Pi_in, int i0_in,
                      const double w_in, const int i_in_in, const char ch_in)
            :                 vertex1(vertex1_in),         vertex2(vertex2_in),           Pi(Pi_in),
            w(w_in), i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0 or 1) into actual Keldysh index i0 (0,...,15)
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
#if DIAG_CLASS <= 1
        // For K1 class, left and right vertices do not depend on integration frequency -> precompute them to save time
        vector<int> indices(2);
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    res_l_V[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 0, channel);
                    res_r_V[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 0, channel);
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    res_l_V[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 0, channel);
                    res_r_V[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 0, channel);
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    res_l_V[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 0, channel);
                    res_r_V[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 0, channel);

                    res_l_Vhat[i2] =  left_same_bare<Q> (vertex1, indices[0], w, 0, i_in, 1, channel);
                    res_r_Vhat[i2] = right_same_bare<Q> (vertex2, indices[1], w, 0, i_in, 1, channel);
                    break;
                default: ;
            }
        }
#endif
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q {
        Q Pival;
        Q res;
#if DIAG_CLASS >= 2
        Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        vector<int> indices(2);
#endif
        //Iterates over all Keldysh components of the bubble which are nonzero
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Flow eq: V*Pi*V
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel
#if DIAG_CLASS <= 1
                    res += res_l_V[i2] * Pival * res_r_V[i2];
#else
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
#endif
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                    Pival = Pi.value(i2, w/2. + vpp, w/2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
#if DIAG_CLASS <= 1
                    res += res_l_V[i2] * Pival * res_r_V[i2];
#else
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    //This is commented out on the ground of p-channel contributions being cross-symmetric
                    //Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                    //the bubble_function(...) must be changed.
//                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
//                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;
#endif
                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
#if DIAG_CLASS <= 1
                    res += res_l_V[i2] * Pival * (res_r_V[i2]+res_r_Vhat[i2])
                            + (res_l_V[i2]+res_l_Vhat[i2]) * Pival * res_r_V[i2];
#else
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    res_l_V =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat =  left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * (res_r_V+res_r_Vhat) + (res_l_V+res_l_Vhat) * Pival * res_r_V;
#endif
                    break;
                default: ;
            }
        }
        return res;
    }
};
template <typename Q> class Integrand_K2_diff {
    const Vertex<Q> &vertex1;
    const Vertex<Q> &vertex2;
    const Bubble& Pi;
    int i0;
    const int i_in;
    const char channel;
    const double w, v;
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
    Integrand_K2_diff(const Vertex<Q> &vertex1_in, const Vertex<Q> &vertex2_in, const Bubble& Pi_in, int i0_in,
                      const double w_in, const double v_in, const int i_in_in, const char ch_in)
            :                 vertex1(vertex1_in),         vertex2(vertex2_in),           Pi(Pi_in),
            w(w_in),     v(v_in), i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0,...,4) into actual Keldysh index i0 (0,...,15)
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
    auto operator() (double vpp) const -> Q {

        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        //Iterates over all Keldysh components of the bubble which are nonzero
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V; + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, w/2. + vpp, w/2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 0, channel);

                    //This is commented out on the ground of p-channel contributions being cross-symmetric
                    //Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                    //the bubble_function(...) must be changed.
                    res_l_Vhat =  left_diff_bare<Q> (vertex1, indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q> (vertex2, indices[1], w,    vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
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
    const Vertex<Q> & vertex1;
    const Vertex<Q> & vertex2;
    const Bubble& Pi;
    int i0;
    const int i_in;
    const char channel;
    const double w, v, vp;
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
    Integrand_K3_diff(const Vertex<Q> &vertex1_in, const Vertex<Q> &vertex2_in, const Bubble& Pi_in, int i0_in,
                      const double w_in, const double v_in, const double vp_in, const int i_in_in, const char ch_in)
            :                 vertex1(vertex1_in),         vertex2(vertex2_in),           Pi(Pi_in),
            w(w_in),     v(v_in),    vp(vp_in), i_in(i_in_in), channel(ch_in)
    {
        // converting index i0_in (0,...,5) into actual Keldysh index i0 (0,...,15)
        i0 = non_zero_Keldysh_K3[i0_in];
    };

    /**
     * Call operator:
     * @param vpp : frequency at which to evaluate integrand (to be integrated over)
     * @return Q  : value of the integrand object evaluated at frequency vpp (comp or double)
     */
    auto operator() (double vpp) const -> Q {
        Q res, res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        Q Pival;
        vector<int> indices(2);
        //Iterates over all Keldysh components of the bubble which are nonzero
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Flow eq: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                    res += res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V// + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, w/2. + vpp, w/2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 0, channel);
                    res_r_V = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 0, channel);

                    //This is commented out on the ground of p-channel contributions being cross-symmetric
                     //Should this not hold, must return to calculating this too, bearing in mind that the prefactor in
                     //the bubble_function(...) must be changed.
//                    res_l_Vhat =  left_diff_bare<Q> (vertex1, indices[0], w, v,  vpp, i_in, 1, channel);
//                    res_r_Vhat = right_diff_bare<Q> (vertex2, indices[1], w, vp, vpp, i_in, 1, channel);

                    res += res_l_V * Pival * res_r_V + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Flow V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w/2., vpp + w/2., i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
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
void bubble_function(Vertex<Q>& dgamma, const Vertex<Q>& vertex1, const Vertex<Q>& vertex2,
                     const Propagator& G, const Propagator& S, const char channel, const bool diff, const char part)
{
    Bubble Pi(G, S, diff); // initialize bubble object

    int nw1_w = 0, nw2_w = 0, nw2_v = 0, nw3_w = 0, nw3_v = 0, nw3_v_p = 0;
    Q prefactor = 1.;

    // set channel-specific frequency ranges and prefactor (1, 1, -1 for a, p, t) for sum over spins.
    switch (channel) {
        case 'a':
            nw1_w = nw1_wa;
            nw2_w = nw2_wa;
            nw2_v = nw2_va;
            nw3_w = nw3_wa;
            nw3_v = nw3_va;
            nw3_v_p = nw3_vap;
            prefactor *= 1.;
            break;
        case 'p':
            nw1_w = nw1_wp;
            nw2_w = nw2_wp;
            nw2_v = nw2_vp;
            nw3_w = nw3_wp;
            nw3_v = nw3_vp;
            nw3_v_p = nw3_vpp;
            prefactor *= 1.;
            break;
        case 't':
            nw1_w = nw1_wt;
            nw2_w = nw2_wt;
            nw2_v = nw2_vt;
            nw3_w = nw3_wt;
            nw3_v = nw3_vt;
            nw3_v_p = nw3_vtp;
            prefactor *= -1.;
            break;
        default: ;
    }

    int mpi_size = mpi_world_size(); // number of mpi processes
    int mpi_rank = mpi_world_rank(); // number of the current mpi process

#ifdef DIAG_CLASS
#if DIAG_CLASS>=0
//    double tK1 = get_time();
    /*K1 contributions*/
    int n_mpi = 1;                      // set external arguments for MPI-parallelization (# of tasks distributed via MPI)
    int n_omp = nK_K1 * nw1_w * n_in;   // set external arguments for OMP-parallelization (# of tasks per MPI-task distributed via OMP)

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
                    value = prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_K1, w_lower_f, w_upper_f, -w/2., w/2.);                      //Integration over a fermionic frequency
                }
                else{
                    Integrand_K1<Q> integrand_K1(vertex1, vertex2, Pi, i0, w, i_in, channel);
                    value  = prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_K1, w_lower_f, w_upper_f, -w/2., w/2.);                      //Integration over a fermionic frequency
                    value += prefactor*(1./(2.*M_PI*glb_i))*asymp_corrections_K1(vertex1, vertex2, w_upper_b, w_upper_b, w, i0, i_in, channel); //Correction needed for the K1 class
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
        case 'a': dgamma[0].avertex.K1 += K1_ordered_result; break;
        case 'p': dgamma[0].pvertex.K1 += K1_ordered_result; break;
        case 't': dgamma[0].tvertex.K1 += K1_ordered_result; break;
        default: ;
    }

    // dgamma[0].pvertex.K1_addvert(i0, iwp, i_in, value); // old version w/o mpi

//    print("K1", channel, " done: ");
//    get_time(tK1);
#endif

#if DIAG_CLASS>=2
//    double tK2 = get_time();
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
                double w = bfreqs2[iw];
                double v = ffreqs2[iv];
                Q value;

                // initialize the integrand object and perform frequency integration
                // (distinguishing between differentiated and non-differentiated bubble)
                if(diff){
                    Integrand_K2_diff<Q> integrand_K2 (vertex1, vertex2, Pi, i0, w, v, i_in, channel);
                    value = prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_K2, w_lower_f, w_upper_f, -w/2., w/2.);                      //Integration over vppp, a fermionic frequency
                }
                else{
                    if (part != 'L') value = 0.;  // right part of multi-loop contribution does not contribute to K2 class
                    // TODO: attention: central part does contribute, but we treat it as right part of previous loop --> fix this!! --> ?
                    else {
                        Integrand_K2<Q> integrand_K2(vertex1, vertex2, Pi, i0, w, v, i_in, channel, part);
                        value = prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_K2, w_lower_f, w_upper_f, -w/2., w/2.);                      //Integration over vppp, a fermionic frequency
                    }
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
        case 'a': dgamma[0].avertex.K2 += K2_ordered_result; break;
        case 'p': dgamma[0].pvertex.K2 += K2_ordered_result; break;
        case 't': dgamma[0].tvertex.K2 += K2_ordered_result; break;
        default: ;
    }

//    print("K2", channel, " done: ");
//    get_time(tK2);
#endif

#if DIAG_CLASS>=3
    double tK3 = get_time();
    /*K3 contributions*/
    n_mpi = nK_K3 * nw3_w * nw3_v;
    n_omp = nw3_v_p * n_in;

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
                int i0 = iK3/(nw3_w * nw3_v * nw3_v_p * n_in);
                int iw = iK3/(nw3_v * nw3_v_p * n_in) - i0*nw3_w;
                int iv = iK3/(nw3_v * n_in) - i0*nw3_w*nw3_v - iw*nw3_v;
                int ivp =iK3/(n_in) - i0*nw3_w*nw3_v*nw3_v_p - iw*nw3_v*nw3_v_p - iv*nw3_v_p;
                int i_in = iK3 - i0*nw3_w*nw3_v*nw3_v_p*n_in - iw*nw3_v*nw3_v_p*n_in - iv*nw3_v_p*n_in - ivp*n_in;
                double w = bfreqs3[iw];
                double v = ffreqs3[iv];
                double vp = ffreqs3[ivp];
                Q value;

                // initialize the integrand object and perform frequency integration
                // (distinguishing between differentiated and non-differentiated bubble)
                if (diff){
                    Integrand_K3_diff<Q> integrand_K3 (vertex1, vertex2, Pi, i0, w, v, vp, i_in, channel);

                    value = prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_K3, w_lower_f, w_upper_f, -w/2., w/2.);                      //Integration over vppp, a fermionic frequency                    // y
                }
                else{
                    Integrand_K3<Q> integrand_K3 (vertex1, vertex2, Pi, i0, w, v, vp,  i_in, channel, part);

                    value = prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_K3, w_lower_f, w_upper_f, -w/2., w/2.);                      //Integration over vppp, a fermionic frequency
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
        case 'a': dgamma[0].avertex.K3 += K3_ordered_result; break;
        case 'p': dgamma[0].pvertex.K3 += K3_ordered_result; break;
        case 't': dgamma[0].tvertex.K3 += K3_ordered_result; break;
        default: ;
    }

    // dgamma[0].pvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, value); // old version w/o mpi

    print("K3", channel, " done: ");
    get_time(tK3);
#endif

#if DIAG_CLASS>=4
    print("Damn son, this is a bad error");
#endif
#endif
}

#endif //KELDYSH_MFRG_BUBBLES_H
