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


//TODO If structure of ALL bubbles is equal, why not have a master-bubble class, takes in g and s and a parameter to choose if regular or differentiated and calculates value accordingly?
/*Class defining the bubble object (two propagators) with a Keldysh structure*/
class Bubble{
    Propagator& g;
public:
    explicit Bubble(Propagator& propagator)
            : g(propagator){};

    /*This function returns the value of the bubble (two propagators) for the Keldysh index iK and propagators frequencies v1 and v2*/
    auto value(int iK, double v1, double v2) -> comp
    {
        comp ans;
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
        return ans;
    }
};

/*Class defining the differentiated bubble object with a Keldysh structure*/
class Diff_Bubble{
    Propagator& g;
    Propagator& s;
public:
    Diff_Bubble(Propagator& propagatorG, Propagator& propagatorS)
            : g(propagatorG), s(propagatorS)    {};

    /*This function returns the value of the differentiated bubble (two propagators) for the Keldysh index iK and propagators frequencies v1 and v2*/
    auto value(int iK, double v1, double v2) -> comp
    {
        comp ans;
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
        return ans;
    }
};

template <typename Q, char channel> class Integrand_K1 { // this class is only for 2OPT
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    double w;
public:
    Integrand_K1(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in,   int i_in_in)
             :                vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                w(w_in),  i_in(i_in_in)
                {
                    switch (channel) {
                        case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
                        case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
                        case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
                        default: ;
                    }
                };

    auto operator() (double vpp) -> Q {
        int i1, i3;
        Q resp;
        comp Pival;
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);       //vppa-1/2wa, vppa+1/2wa for the a-channel
                    tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
                    break;
                case 'p':
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);       //2wp/2+vppp, wp/2-vppp for the p-channel
                    tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);
                    break;
                case 't':
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);       //vppt-1/2wt, vppt+1/2wt for the t-channel
                    tie(i1, i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);
                    break;
                default: ;
            }

            Q add1 = vertex1.densvertex.irred.vval(i1, i_in) * Pival * vertex2.densvertex.irred.vval(i3, i_in);

            resp += add1;
            //Augments to RPA
//            resp += vertex1.densvertex.irred.vval(i1, i_in) * PiPval * vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
//            resp += vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) * PiPval * vertex2.densvertex.irred.vval(i3, i_in);
//            resp += vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) * PiPval * vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
        }
        return resp;
    }

};


template <typename Q, char channel, char part> class Integrand_K2
{
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble &Pi;
    int i0, i_in;
    double w, v;
public:
    Integrand_K2(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in, double v_in,   int i_in_in)
             :                vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),                w(w_in),     v(v_in), i_in(i_in_in)
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

        int i1=0, i3=0;
        Q res, res_l, res_r;
        Q Pival;
        for (auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':
                    tie(i1, i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    break;
                case 'p':
                    tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    break;
                case 't':
                    tie(i1, i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);
                    Pival = Pi.value(i2, vpp - 0.5 * w, vpp + 0.5 * w);
                    break;
                default: ;
            }
            res_l = vertex1.densvertex.gammRb(i1, w, v, vpp, i_in, channel);
            res_r = right_same_bare<Q, channel> (vertex2, i3, w, vpp, i_in);
            res += res_l * Pival * res_r;
        }
        return res;
    }
};
template <typename Q, typename Bubble> class Integrand_p_K3
{
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp, vp, vpp;
public:
    Integrand_p_K3(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiP_in, int i0_in, double wp_in, double vp_in, double vpp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in),   vpp(vpp_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    auto operator() (double vppp) -> Q {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K3: (gammaP)Pi(gammaP)
            resp += (vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

};

template <typename Q, char channel> class Integrand_K1_diff { // TODO: finish
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& Pi;
    int i0, i_in;
    double w;
public:
    Integrand_K1_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& Pi_in, int i0_in, double w_in, int i_in_in)
             :                     vertex1(vertex1_in),              vertex2(vertex2_in),     Pi(Pi_in),               w(w_in), i_in(i_in_in)
    {
        switch (channel) {
            case 'a': i0 = non_zero_Keldysh_K1a[i0_in]; break;
            case 'p': i0 = non_zero_Keldysh_K1p[i0_in]; break;
            case 't': i0 = non_zero_Keldysh_K1t[i0_in]; break;
            default: ;
        }
    };

    //This is a call operator
    auto operator() (double vpp) -> Q {
        int i1=0, i3=0;
        Q res, res_l, res_r;
        Q Pival;
        for(auto i2:non_zero_Keldysh_bubble) {
            switch (channel) {
                case 'a':
                    tie(i1,i3) = vertex1.densvertex.avertex.indices_sum(i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppa-1/2wa, vppa+1/2wa for the a-channel
                case 'p':
                    tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);
                    Pival = Pi.value(i2, 0.5*w+vpp, 0.5*w-vpp);                                //wp/2+vppp, wp/2-vppp for the p-channel
                case 't':
                    tie(i1,i3) = vertex1.densvertex.tvertex.indices_sum(i0, i2);
                    Pival = Pi.value(i2, vpp-0.5*w, vpp+0.5*w);                                //vppt-1/2wt, vppt+1/2wt for the t-channel
                default: ;
            }
            res_l = left_same_bare<comp, channel> (vertex1, i1, w, vpp, i_in);
            res_r = right_same_bare<comp, channel> (vertex2, i3, w, vpp, i_in);

            res += res_l * Pival * res_r;

//            resp += vertex1.densvertex.irred.vval(i1, i_in) * PiPval * vertex2.densvertex.irred.vval(i3, i_in);

//            //Contributions to K1: (u+K1+K2b)Pi(u+K1+K2)
//            resp += (vertex1.densvertex.irred.vval(i1, i_in) +
//                     vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in)  +
//                     vertex1.densvertex.pvertex.K2b_vvalsmooth(i1, wp, vppp, i_in)) *
//                    PiP.value(i2, 0.5*wp+vppp, 0.5*wp-vppp) *                                       //wp/2+vppp, wp/2-vppp for the p-channel
//                    (vertex2.densvertex.irred.vval(i3, i_in) +
//                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
//                     vertex2.densvertex.pvertex.K2_vvalsmooth (i3, wp, vppp, i_in) );

            //TODO: put this into functions
            resp1 = vertex1.densvertex.irred.vval(i1, i_in);
            resp2 = vertex1.densvertex.pvertex.K1_vvalsmooth(i1, wp, i_in) ;
            resp3 = vertex1.densvertex.pvertex.K2b_vvalsmooth(i1, wp, vppp, i_in);
            resp4 = vertex2.densvertex.irred.vval(i3, i_in);
            resp5 = vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
            resp6 = vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in);

            resp += (resp1 + resp2 + resp3) * PiPval * (resp4 + resp5 + resp6);
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_p_K2_diff {
    Vertex<fullvert<Q> > &vertex1;
    Vertex<fullvert<Q> > &vertex2;
    Bubble &PiP;
    int i0, i_in;
    double wp, vp;
public:
    Integrand_p_K2_diff(Vertex<fullvert<Q> > &vertex1_in, Vertex<fullvert<Q> > &vertex2_in, Bubble &PiP_in, int i0_in, double wp_in, double vp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in), i_in(non_zero_Keldysh_K2p[i_in_in]) {};

    //This is a call operator
    auto operator() (double vppp) -> Q {
        int i1, i3;
        Q resp;
        Q resp1, resp2, resp3, resp4, resp5, resp6;
        for (auto i2:non_zero_Keldysh_bubble) {
            tie(i1, i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            auto PiPvalue = PiP.value(i2, 0.5*wp+vppp, 0.5*wp - vppp);

            //Contributions to K2: (K2 +K3 + gammaP)Pi(K1+K2)
//            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) +
//                     vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
//                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
//                    PiP.value(i2, 0.5*wp+vppp, 0.5*wp-vppp) *                                       //wp/2+vppp, wp/2-vppp for the p-channel
//                    (vertex2.densvertex.irred.vval(i3, i_in) +
//                     vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in) +
//                     vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in));

            resp1 = vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in);
            resp2 = vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in);
            resp3 = vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p');
            resp4 = vertex2.densvertex.irred.vval(i3, i_in);
            resp5 = vertex2.densvertex.pvertex.K1_vvalsmooth(i3, wp, i_in);
            resp6 = vertex2.densvertex.pvertex.K2_vvalsmooth(i3, wp, vppp, i_in);

            resp += (resp1+resp2+resp3)*PiPvalue*(resp4+resp5+resp6);
        }
        return resp;
    }
};
template <typename Q, typename Bubble> class Integrand_p_K3_diff {
    Vertex<fullvert<Q> >& vertex1;
    Vertex<fullvert<Q> >& vertex2;
    Bubble& PiP;
    int i0, i_in;
    double wp, vp, vpp;
public:
    Integrand_p_K3_diff(Vertex<fullvert<Q> >& vertex1_in, Vertex<fullvert<Q> >& vertex2_in, Bubble& PiP_in, int i0_in, double wp_in, double vp_in, double vpp_in, int i_in_in)
            :                   vertex1(vertex1_in),              vertex2(vertex2_in),    PiP(PiP_in), i0(i0_in),    wp(wp_in),    vp(vp_in),   vpp(vpp_in), i_in(non_zero_Keldysh_K3[i_in_in]) {};

    //This is a call operator
    auto operator() (double vppp) -> Q {
        int i1, i3;
        Q resp;
        for(auto i2:non_zero_Keldysh_bubble) {
            tie(i1,i3) = vertex1.densvertex.pvertex.indices_sum(i0, i2);

            //Contributions to K3: (K2 +K3 + gammaP)Pi(K2b + K3 + gammaP)
            resp += (vertex1.densvertex.pvertex.K2_vvalsmooth(i1, wp, vp, i_in) +
                     vertex1.densvertex.pvertex.K3_vvalsmooth(i1, wp, vp, vppp, i_in) +
                     vertex1.densvertex.gammaRb(i1, wp, vp, vppp, i_in, 'p')) *
                    PiP.value(i2, vppp-0.5*wp, vppp+0.5*wp) *
                    (vertex2.densvertex.pvertex.K2b_vvalsmooth(i3, wp, vp, i_in) +
                     vertex2.densvertex.pvertex.K3_vvalsmooth(i3, wp, vppp, vpp, i_in) +
                     vertex2.densvertex.gammaRb(i3, wp, vppp, vpp, i_in, 'p'));
        }
        return resp;
    }

};



template <typename Q, char channel, bool diff>
void bubble_function(Vertex<fullvert<Q> >& dgamma, Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, Propagator& S)
{
    Bubble Pi(G);
    if (diff) Diff_Bubble Pi(G,S);

    int nw1_w = 0, nw2_w = 0, nw2_nu = 0, nw3_w = 0, nw3_nu = 0, nw3_nu_p = 0;
    Q prefactor;

    switch (channel) {
        case 'a':
            nw1_w = nw1_wa;
            nw2_w = nw2_wa;
            nw2_nu = nw2_nua;
            nw3_w = nw3_wa;
            nw3_nu = nw3_nua;
            nw3_nu_p = nw3_nuap;
            prefactor = 1.;
            break;
        case 'p':
            nw1_w = nw1_wp;
            nw2_w = nw2_wp;
            nw2_nu = nw2_nup;
            nw3_w = nw3_wp;
            nw3_nu = nw3_nup;
            nw3_nu_p = nw3_nupp;
            prefactor = 0.5;
            break;
        case 't':
            nw1_w = nw1_wt;
            nw2_w = nw2_wt;
            nw2_nu = nw2_nut;
            nw3_w = nw3_wt;
            nw3_nu = nw3_nut;
            nw3_nu_p = nw3_nutp;
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
                int i0 = (iK1 % (nK_K1 * nw1_w * n_in)) / (nw1_w * n_in);
                int iwp = (iK1 % (nw1_w * n_in)) / n_in;
                int i_in = iK1 % n_in;
                double wp = bfreqs[iwp];

                Integrand_K1<Q, channel> integrand_K1 (vertex1, vertex2, Pi, i0, wp, i_in);
                if (diff) Integrand_K1_diff<Q, channel> integrand_K1 (vertex1, vertex2, Pi, i0, wp, i_in);

                Q value = prefactor*(-1./(2.*pi*(comp)1.i))*integrator(integrand_K1, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency

                K1_buffer[iterator*n_omp + i_omp] = value;
            }
            ++iterator;
        }
    }

    vec<Q> K1_result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(K1_buffer, K1_result, n_mpi, n_omp);
    vec<Q> K1_ordered_result = mpi_reorder_result(K1_result, n_mpi, n_omp);

    switch (channel) {
        case 'a': dgamma.densvertex.avertex.K1 += K1_ordered_result; break;
        case 'p': dgamma.densvertex.pvertex.K1 += K1_ordered_result; break;
        case 't': dgamma.densvertex.tvertex.K1 += K1_ordered_result; break;
        default: ;
    }

    // dgamma.densvertex.pvertex.K1_addvert(i0, iwp, i_in, value); // old version w/o mpi

    cout << "K1p done: ";
    get_time(tK1);
#endif

#if DIAG_CLASS>=2
    double tK2 = get_time();
    /*K2 contributions*/
    n_mpi = nK_K2 * nw2_w;
    n_omp = nw2_nu * n_in;

    vec<Q> K2_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                int iK2 = i_mpi * n_omp + i_omp;
                int i0 = (iK2 % (nK_K2 * nw2_w * nw2_nu * n_in)) / (nw2_w * nw2_nu * n_in);
                int iwp = (iK2 % (nw2_w * nw2_nu * n_in)) / (nw2_nu * n_in);
                int ivp = (iK2 % (nw2_nu * n_in)) / n_in;
                int i_in = iK2 % n_in;
                double wp = bfreqs[iwp];
                double vp = ffreqs[ivp];

                Integrand_K2<Q, channel> integrand_K2 (vertex1, vertex2, Pi, i0, wp, vp, i_in);
                if (diff) Integrand_K2_diff<Q, channel> integrand_K2 (vertex1, vertex2, Pi, i0, wp, vp, i_in);

                Q value = prefactor*(-1./(2.*pi*(comp)1.i))*integrator(integrand_K2, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency

                K2_buffer[iterator*n_omp + i_omp] = value;
            }
            ++iterator;
        }
    }

    vec<Q> K2_result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(K2_buffer, K2_result, n_mpi, n_omp);
    vec<Q> K2_ordered_result = mpi_reorder_result(K2_result, n_mpi, n_omp);

    switch (channel) {
        case 'a': dgamma.densvertex.avertex.K2 += K2_ordered_result; break;
        case 'p': dgamma.densvertex.pvertex.K2 += K2_ordered_result; break;
        case 't': dgamma.densvertex.tvertex.K2 += K2_ordered_result; break;
        default: ;
    }

    // dgamma.densvertex.pvertex.K2_addvert(i0, iwp, ivp, i_in, value); // old version w/o mpi

    cout << "K2p done: ";
    get_time(tK2);
#endif

#if DIAG_CLASS>=3
    double tK3 = get_time();
    /*K3 contributions*/
    n_mpi = nK_K3 * nw3_w * nw3_nu;
    n_omp = nw3_nu_p * n_in;

    vec<Q> K3_buffer = mpi_initialize_buffer<Q>(n_mpi, n_omp);

    iterator = 0;
    for (int i_mpi=0; i_mpi<n_mpi; ++i_mpi) {
        if (i_mpi % mpi_size == mpi_rank) {
#pragma omp parallel for
            for (int i_omp=0; i_omp<n_omp; ++i_omp) {
                int iK3 = i_mpi * n_omp + i_omp;
                int i0 = (iK3 % (nK_K3 * nw3_w * nw3_nu * nw3_nu_p * n_in)) / (nw3_w * nw3_nu * nw3_nu_p * n_in);
                int iwp = (iK3 % (nw3_wp * nw3_nu * nw3_nu_p * n_in)) / (nw3_nu * nw3_nu_p * n_in);
                int ivp = (iK3 % (nw3_nu * nw3_nu_p * n_in)) / (nw3_nu_p * n_in);
                int ivpp = (iK3 % (nw3_nu_p * n_in))/ n_in;
                int i_in = iK3 % n_in;
                double wp = bfreqs[iwp];
                double vp = ffreqs[ivp];
                double vpp = ffreqs[ivpp];

                Integrand_K3<Q, channel> integrand_K3 (vertex1, vertex2, Pi, i0, wp, vp, vpp,  i_in);
                if (diff) Integrand_K3_diff<Q, channel> integrand_K3 (vertex1, vertex2, Pi, i0, wp, vp, vpp,  i_in);

                Q value = prefactor*(-1./(2.*pi*(comp)1.i))*integrator(integrand_K3, w_lower_f, w_upper_f);                      //Integration over vppp, a fermionic frequency

                K3_buffer[iterator*n_omp + i_omp] = value;
            }
            ++iterator;
        }
    }

    vec<Q> K3_result = mpi_initialize_result<Q> (n_mpi, n_omp);
    mpi_collect(K3_buffer, K3_result, n_mpi, n_omp);
    vec<Q> K3_ordered_result = mpi_reorder_result(K3_result, n_mpi, n_omp);

    switch (channel) {
        case 'a': dgamma.densvertex.avertex.K3 += K3_ordered_result; break;
        case 'p': dgamma.densvertex.pvertex.K3 += K3_ordered_result; break;
        case 't': dgamma.densvertex.tvertex.K3 += K3_ordered_result; break;
        default: ;
    }

    // dgamma.densvertex.pvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, value); // old version w/o mpi

    cout << "K3p done: ";
    get_time(tK3);
#endif

#if DIAG_CLASS>=4
    cout << "Damn son, this is a bad error";
#endif
#endif
}

template <typename Q> void p_bubble_function(Vertex<fullvert<Q> >& gamma, Vertex<fullvert<Q> >& vertex1, Vertex<fullvert<Q> >& vertex2, Propagator& G, char side)
{
    P_Bubble PiP(G);

    //These lines are to test the SOPT results - there are no K1 contributions in the corrections!
    double t0 = get_time();
    /*K1 contributions*/
#pragma omp parallel for
    for (int iK1=0; iK1<nK_K1*nw1_wp*n_in; ++iK1) {

        int i0 = (iK1 % (nK_K1 * nw1_wp * n_in)) / (nw1_wp * n_in);
        int iwp = (iK1 % (nw1_wp * n_in)) / n_in;
        int i_in = iK1 % n_in;
        double wp = bfreqs[iwp];

        Integrand_p_K1 <Q, P_Bubble> integrand_p_K1 (vertex1, vertex2, PiP, i0, wp, i_in);

        Q value = (0.5)*(-1./(2.*pi*(comp)1.i))*integrator(integrand_p_K1, w_lower_f, w_upper_f);                   //Integration over vppp, a fermionic frequency

        gamma.densvertex.pvertex.K1_addvert(i0, iwp, i_in, value);
    }
    cout << "K1p done:" << endl;
    get_time(t0);

//    if(side == 'L')
//    {
//        //In this case, there are only contributions to K2 and K3
//        /*K2 contributions*/
//#pragma omp parallel for
//        for(int iK2=0; iK2<nK_K2*nw2_wp*nw2_nup*n_in; iK2++)
//        {
//            int i0 = (iK2 % (nK_K2 * nw2_wp * nw2_nup * n_in)) / (nw2_wp * nw2_nup * n_in);
//            int iwp = (iK2 % (nw2_wp * nw2_nup * n_in)) / (nw2_nup * n_in);
//            int ivp = (iK2 % (nw2_nup * n_in)) / n_in;
//            int i_in = iK2 % n_in;
//            double wp = bfreqs[iwp];
//            double vp = ffreqs[ivp];
//
//            Integrand_p_K2<Q, P_Bubble> integrand_p_K2 (vertex1, vertex2, PiP, i0, wp, vp, i_in);
//
//            resp.densvertex.K2_addvert(i0, iwp, ivp, i_in, integrator(integrand_p_K2, ffreqs)); //
//        }
//    }
//    else if(side == 'R')
//    {
//        //In this case, there are only contributions to K2b and K3
//        /*K2b contributions*/
//        //TODO implement
//    }
//    else
//    {
//        return resp;
//    }

    /*K3 contributions*/
//#pragma omp parallel for
//    for(int iK3=0; iK3<nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in; iK3++)
//    {
//        int i0 = (iK3 % (nK_K3 * nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_wp * nw3_nup * nw3_nupp * n_in);
//        int iwp = (iK3 % (nw3_wp * nw3_nup * nw3_nupp * n_in)) / (nw3_nup * nw3_nupp * n_in);
//        int ivp = (iK3 % (nw3_nup * nw3_nupp * n_in)) / (nw3_nupp * n_in);
//        int ivpp = (iK3 % (nw3_nupp * n_in))/ n_in;
//        int i_in = iK3 % n_in;
//        double wp = bfreqs[iwp];
//        double vp = ffreqs[ivp];
//        double vpp = ffreqs[ivpp];
//
//        Integrand_p_K3<Q, P_Bubble> integrand_p_K3 (vertex1, vertex2, PiP, i0, wp, vp, vpp, i_in);
//
//        resp.densvertex.K3_addvert(i0, iwp, ivp, ivpp, i_in, integrator(integrand_p_K3, ffreqs));
//    }

}

#endif //KELDYSH_MFRG_BUBBLES_H
