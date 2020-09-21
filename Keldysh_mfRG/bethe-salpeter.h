//
// Created by SAguirre on 15/08/2020.
//

#ifndef KELDYSH_MFRG_BETHE_SALPETER_H
#define KELDYSH_MFRG_BETHE_SALPETER_H

#include "parameters.h"
#include "state.h"
#include "diagrammatic_combinations.h"
#include "integrator.h"
#include "bubbles.h"
#include "solvers.h"
#include "hdf5_routines.h"
#include "data_structures.h"
#include "write_data2file.h"
#include "loop.h"
#include <iostream>

rvec reconstruct_grid(){
    const double X_ini = log_substitution(Lambda_ini), X_fin = log_substitution(Lambda_fin); // substitute limits
    const double dX = (X_fin-X_ini)/((double)nODE);         // equidistant grid in substituted variable X

    // create non-linear integration grid using substitution
    rvec x_vals (nODE+1);               // integration values
    x_vals[0] = Lambda_ini;                          // start with initial value
    for (int i=1; i<=nODE; ++i) {
        x_vals[i] = log_resubstitution(X_ini + i*dX);      // value i
    }
    add_points_to_Lambda_grid(x_vals);
    return x_vals;
}

/// Integrand classes for non-differentiated bubble contributing to diagrammatic class K1, K2, K3
template <typename Q> class Integrand_BS_K1 {
    const Vertex<Q>& vertex1;
    const Vertex<Q>& vertex2;
    const Bubble& Pi;
    int i0;
    int i2;
    const double w;
    const int i_in;
    const char channel;
    const char side;

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
    Integrand_BS_K1(const Vertex<Q>& vertex1_in, const Vertex<Q>& vertex2_in, const Bubble& Pi_in,
                 int i0_in, int i2_in, const double w_in, const int i_in_in, const char ch_in, const char side_in)
            : vertex1(vertex1_in),         vertex2(vertex2_in),           Pi(Pi_in),
              i2(i2_in), w(w_in), i_in(i_in_in), channel(ch_in), side(side_in)
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
    auto operator() (double vpp) const -> Q {
        Q Pival;
        Q res;
        Q res_l_V, res_r_V, res_l_Vhat, res_r_Vhat;
        vector<int> indices(2);


        if(side == 'L') { //In this case, we ant to have I_r on the left side
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Flow eq: V*Pi*V
                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2., i_in);             //vppa-1/2wa, vppa+1/2wa for the a-channel
                    vertex1[0].avertex.indices_sum(indices, i0, i2);

                    res_l_V = vertex1[0].irred.val(indices[0], i_in, 0);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res = res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V //+ V^*Pi*V^

                    Pival = Pi.value(i2, w / 2. + vpp, w / 2. - vpp, i_in);              //wp/2+vppp, wp/2-vppp for the p-channel
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    res_l_V = vertex1[0].irred.val(indices[0], i_in, 0);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res = res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;

                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V

                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2., i_in);              //vppt-1/2wt, vppt+1/2wt for the t-channel

                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    res_l_V = vertex1[0].irred.val(indices[0], i_in, 0);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1[0].irred.val(indices[0], i_in, 1);
                    res_r_Vhat = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
                    break;
                default:;
            }
        }
        else{
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Flow eq: V*Pi*V
                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2., i_in);             //vppa-1/2wa, vppa+1/2wa for the a-channel
                    vertex1[0].avertex.indices_sum(indices, i0, i2);

                    res_l_V = left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = vertex2[0].irred.val(indices[1], i_in, 0);

                    res = res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Flow eq: V*Pi*V //+ V^*Pi*V^

                    Pival = Pi.value(i2, w / 2. + vpp, w / 2. - vpp, i_in);              //wp/2+vppp, wp/2-vppp for the p-channel
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);

                    res_l_V = left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = vertex2[0].irred.val(indices[1], i_in, 0);

                    res = res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;

                    break;
                case 't':                                                                       //Flow eq: V*Pi*(V+V^) + (V+V^)*Pi*V

                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2., i_in);              //vppt-1/2wt, vppt+1/2wt for the t-channel
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);

                    res_l_V = left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 0, channel);
                    res_r_V = vertex2[0].irred.val(indices[1], i_in, 0);

                    res_l_Vhat = left_same_bare<Q> (vertex1, indices[0], w, vpp, i_in, 1, channel);
                    res_r_Vhat = vertex2[0].irred.val(indices[1], i_in, 1);

                    res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
                    break;
                default:;
            }
        }

        return res;
    }
};

template <typename Q> class Integrand_BS_K2
{
    const Vertex<Q>& vertex1;
    const Vertex<Q>& vertex2;
    const Bubble&Pi;
    int i0;
    int i2;
    const int i_in;
    const char channel, side;
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
    Integrand_BS_K2(const Vertex<Q>& vertex1_in, const Vertex<Q>& vertex2_in, const  Bubble& Pi_in,
                 int i0_in, int i2_in, const double w_in, double v_in, const int i_in_in, const char ch_in, const char side_in)
            : vertex1(vertex1_in), vertex2(vertex2_in), Pi(Pi_in),
              i2(i2_in), w(w_in), v(v_in), i_in(i_in_in), channel(ch_in), side(side_in)

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

        if(side == 'L') {
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Contributions: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2.,
                                     i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel

                    res_l_V = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res = res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Contributions: V*Pi*V// + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, w / 2. + vpp, w / 2. - vpp,
                                     i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel
                    res_l_V = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res = res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2.,
                                     i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel
                    res_l_V = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 0, channel);

                    res_l_Vhat = vertex1[0].gammaRb(indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = right_same_bare<Q>(vertex2, indices[1], w, vpp, i_in, 1, channel);

                    res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
                    break;
                default:;
            }
        }
        else{   //if side == 'R', i.e. I_r is on the right side
            switch (channel) {
                //According to channel, indices of the left and right vertices are determined.
                //Then, the value of the multiplication of the two propagators is calculated.
                //Left and right values of the vertices are determined. Keep in mind which spin components of the vertex
                //contribute to the relevant spin components
                //Add contribution to the result.
                case 'a':                                                                       //Contributions: V*Pi*V
                    vertex1[0].avertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2., i_in);                         //vppa-1/2wa, vppa+1/2wa for the a-channel

                    res_l_V = left_diff_bare<Q>(vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = vertex2[0].irred.val(indices[1], i_in, 0);

                    res = res_l_V * Pival * res_r_V;
                    break;
                case 'p':                                                                       //Contributions: V*Pi*V// + V^*Pi*V^
                    vertex1[0].pvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, w / 2. + vpp, w / 2. - vpp, i_in);                         //wp/2+vppp, wp/2-vppp for the p-channel

                    res_l_V = left_diff_bare<Q>(vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = vertex2[0].irred.val(indices[1], i_in, 0);

                    res = res_l_V * Pival * res_r_V;// + res_l_Vhat * Pival * res_r_Vhat;
                    break;
                case 't':                                                                       //Contributions: V*Pi*(V+V^) + (V+V^)*Pi*V
                    vertex1[0].tvertex.indices_sum(indices, i0, i2);
                    Pival = Pi.value(i2, vpp - w / 2., vpp + w / 2.,i_in);                         //vppt-1/2wt, vppt+1/2wt for the t-channel

                    res_l_V = left_diff_bare<Q>(vertex1, indices[0], w, v, vpp, i_in, 0, channel);
                    res_r_V = vertex2[0].irred.val(indices[1], i_in, 0);

                    res_l_V = left_diff_bare<Q>(vertex1, indices[0], w, v, vpp, i_in, 1, channel);
                    res_r_Vhat = vertex2[0].irred.val(indices[1], i_in, 1);

                    res = res_l_V * Pival * (res_r_V + res_r_Vhat) + (res_l_V + res_l_Vhat) * Pival * res_r_V;
                    break;
                default:;
            }
        }
        return res;
    }
};

template <typename Q>
void Bethe_Salpeter_bubble(Vertex<Q>& lhs_vertex, const Vertex<Q>& vertex1, const Vertex<Q>& vertex2,
                           const Propagator& G, const Propagator& S, const char channel, const char side){

    Bubble Pi(G, S, false); // initialize bubble object

    int nw1_w = 0, nw2_w = 0, nw2_v = 0, nw3_w = 0, nw3_v = 0, nw3_v_p = 0;
    Q prefactor = 1.;

    // set channel-specific frequency ranges and prefactor (1, 1, -1 for a, p, t) for sum over spins.
    switch (channel) {
        case 'a':
            nw1_w = nw1_a;
            nw2_w = nw2_a;
            nw2_v = nv2_a;
            nw3_w = nw3_a;
            nw3_v = nv3_a;
            nw3_v_p = nv3_a;
            prefactor *= 1.;
            break;
        case 'p':
            nw1_w = nw1_p;
            nw2_w = nw2_p;
            nw2_v = nv2_p;
            nw3_w = nw3_p;
            nw3_v = nv3_p;
            nw3_v_p = nv3_p;
            prefactor *= 1.;
            break;
        case 't':
            nw1_w = nw1_t;
            nw2_w = nw2_t;
            nw2_v = nv2_t;
            nw3_w = nw3_t;
            nw3_v = nv3_t;
            nw3_v_p = nv3_t;
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
                // converting external MPI/OMP indices to physical indices
                int iK1 = i_mpi * n_omp + i_omp;
                int i0 = iK1/(nw1_w*n_in);
                int iw = iK1/(n_in) - i0*nw1_w;
                int i_in = iK1 - i0*nw1_w*n_in - iw*n_in;
                double w = bfreqs[iw];
                Q value;

                for(auto i2: non_zero_Keldysh_bubble){
                    Integrand_BS_K1<Q> integrand_BS_K1 (vertex1, vertex2, Pi, i0, i2, w, i_in, channel, side);

                    value += prefactor * (1. / (2. * M_PI * glb_i))*integrator(integrand_BS_K1, glb_v_lower, glb_v_upper);   //Integration over a fermionic frequency
//                    value += prefactor * (1. / (2. * M_PI * glb_i))*asymp_corrections_K1(vertex1, vertex2, -glb_v_lower, glb_v_upper, w, i0, i2, i_in, channel);
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
        case 'a': lhs_vertex[0].avertex.K1 += K1_ordered_result; break;
        case 'p': lhs_vertex[0].pvertex.K1 += K1_ordered_result; break;
        case 't': lhs_vertex[0].tvertex.K1 += K1_ordered_result; break;
        default: ;
    }
#endif


#if DIAG_CLASS>=2
//    double tK2 = get_time();
    /*K2 contributions*/
    n_mpi = nK_K2;
    n_omp = nw2_w * nw2_v * n_in;

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

                for(auto i2:non_zero_Keldysh_bubble){
                    Integrand_BS_K2<Q> integrand_BS_K2 (vertex1, vertex2, Pi, i0, i2, w, v, i_in, channel, side);

                    value += prefactor*(1./(2.*M_PI*glb_i))*integrator(integrand_BS_K2, glb_v_lower, glb_v_upper);
                    value += prefactor*(1./(2.*M_PI*glb_i))*asymp_corrections_K2(vertex1, vertex2, -glb_v_lower, glb_v_upper, w, v, i0, i2, i_in, channel);
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
        case 'a': lhs_vertex[0].avertex.K2 += K2_ordered_result; break;
        case 'p': lhs_vertex[0].pvertex.K2 += K2_ordered_result; break;
        case 't': lhs_vertex[0].tvertex.K2 += K2_ordered_result; break;
        default: ;
    }
#endif

#endif
}



Vertex<comp> calculate_Bethe_Salpeter_equation(const Vertex<comp>& Gamma, const Propagator& G, const char side){

    Vertex<comp> bethe_salpeter(n_spin);
    bethe_salpeter[0].irred.initialize(-glb_U/2.);

    Bethe_Salpeter_bubble(bethe_salpeter, Gamma, Gamma, G, G, 'a', side);
    Bethe_Salpeter_bubble(bethe_salpeter, Gamma, Gamma, G, G, 'p', side);
    Bethe_Salpeter_bubble(bethe_salpeter, Gamma, Gamma, G, G, 't', side);

    return bethe_salpeter;

}

void check_BSE_and_SDE(const string& dir, const H5std_string& filename){

    rvec lambdas = reconstruct_grid();
    rvec actual_gird;       //lambdas has length nODE + U_NRG.size() which, for Gamma>0 is larger than needed!
    rvec norm_K1_fRG;
    rvec norm_K1_BSE;
    rvec norm_K2_fRG;
    rvec norm_K2_BSE;
    rvec K1_diff_rel;
    rvec K1_diff_abs;
    rvec K2_diff_rel;
    rvec K2_diff_abs;
    rvec norm_S_fRG;
    rvec norm_S_SDE;
    rvec Sigma_diff_rel;
    rvec Sigma_diff_abs;

    bool done = false;

    for(int i=0; i<lambdas.size(); i++){

        if(done){
            break;
        }

        if(lambdas[i] == 0){
            done = true;
        }

        actual_gird.emplace_back(lambdas[i]);

        print("Current nLambda = " + to_string(i)+". Represents Lambda = " + to_string(lambdas[i]), true);

        State<comp> state = read_hdf(dir + filename, i, lambdas.size());

        Propagator G (lambdas[i], state.selfenergy, 'g');

        //Bethe-Salpeter equation for gamma_r. Due to channel-crossing symmetries, only consider full vertex
        Vertex<comp> bethe_salpeter_L = calculate_Bethe_Salpeter_equation(state.vertex, G, 'L');
        Vertex<comp> bethe_salpeter_R = calculate_Bethe_Salpeter_equation(state.vertex, G, 'R');
        Vertex<comp> bethe_salpeter_vertex = (bethe_salpeter_L + bethe_salpeter_R) * 0.5;

        Vertex<comp> fRG_vertex = state.vertex;
        Vertex<comp> vertex_diff = fRG_vertex - bethe_salpeter_vertex;

        //p of the p-norm
        int p = 2;

#if DIAG_CLASS >= 0
        norm_K1_fRG.emplace_back(fRG_vertex[0].norm_K1(p));
        norm_K1_BSE.emplace_back(bethe_salpeter_vertex[0].norm_K1(p));
        K1_diff_rel.emplace_back(vertex_diff[0].norm_K1(p)/fRG_vertex[0].norm_K1(p));
        K1_diff_abs.emplace_back(vertex_diff[0].norm_K1(p));
#endif
#if DIAG_CLASS >= 2
        norm_K2_fRG.emplace_back(fRG_vertex[0].norm_K2(p));
        norm_K2_BSE.emplace_back(bethe_salpeter_vertex[0].norm_K2(p));
        K2_diff_rel.emplace_back(vertex_diff[0].norm_K2(p)/fRG_vertex[0].norm_K2(p));
        K2_diff_abs.emplace_back(vertex_diff[0].norm_K2(p));
#endif

        //Schwinger-Dyson comparison
        Vertex<comp> bare_vertex (n_spin);
        bare_vertex[0].initialize(-glb_U/2.);

        Vertex<comp> temp_bubble_a(n_spin);
        //Calculate K1-contributions. part='.' yields only K1
        bubble_function(temp_bubble_a, bare_vertex, fRG_vertex, G, G, 'a', false, '.');
        bubble_function(temp_bubble_a, fRG_vertex, bare_vertex, G, G, 'a', false, '.');

        //Calculare K2-contributions. part='L' only one to yield non-zero result
//        bubble_function(temp_bubble_a, bare_vertex, fRG_vertex, G, G, 'a', false, 'R');
        bubble_function(temp_bubble_a, fRG_vertex, bare_vertex, G, G, 'a', false, 'L');
        temp_bubble_a *= 0.5;

        Vertex<comp> temp_bubble_p(n_spin);
        //Calculate K1-contributions. part='.' yields only K1
        bubble_function(temp_bubble_p, bare_vertex, fRG_vertex, G, G, 'p', false, '.');
        bubble_function(temp_bubble_p, fRG_vertex, bare_vertex, G, G, 'p', false, '.');

        //Calculare K2-contributions. part='L' only one to yield non-zero result
//        bubble_function(temp_bubble_p, bare_vertex, fRG_vertex, G, G, 'p', false, 'R');
        bubble_function(temp_bubble_p, fRG_vertex, bare_vertex, G, G, 'p', false, 'L');
        temp_bubble_p *= 0.5;


        SelfEnergy<comp> schwinger_dyson_p;
        SelfEnergy<comp> schwinger_dyson_a;
        SelfEnergy<comp> schwinger_dyson_self_energy;
        SelfEnergy<comp> fRG_self_energy = state.selfenergy;


        schwinger_dyson_p.initialize(glb_U/2., 0.);   //Bare-vertex contribution known analytically for ph-symmetric
        loop(schwinger_dyson_p, temp_bubble_p, G, false);

        schwinger_dyson_a.initialize(glb_U/2., 0.);   //Bare-vertex contribution known analytically for ph-symmetric
        loop(schwinger_dyson_a, temp_bubble_a, G, false);

        schwinger_dyson_self_energy = (schwinger_dyson_a + schwinger_dyson_p)*0.5;

        SelfEnergy<comp> self_energy_diff = fRG_self_energy - schwinger_dyson_self_energy;
        SelfEnergy<comp> hartree;
        fRG_self_energy.asymp_val_R = glb_U/2.;     //Info in sot read off of HDF5 data, must be manually input
        hartree.initialize(fRG_self_energy.asymp_val_R, 0.);

        norm_S_fRG.emplace_back( (fRG_self_energy-hartree).norm(p));
        norm_S_SDE.emplace_back( (schwinger_dyson_self_energy-hartree).norm(p));
        Sigma_diff_rel.emplace_back(self_energy_diff.norm(p)/((fRG_self_energy-hartree).norm(p)));
        Sigma_diff_abs.emplace_back(self_energy_diff.norm(p));

        State<comp> parquet;
        parquet.vertex = bethe_salpeter_vertex;
        parquet.selfenergy = schwinger_dyson_self_energy;

        save_to_hdf(dir + "parquet_check_" + filename, i, lambdas.size(), parquet, lambdas, (bool)i);

        write_h5_rvecs(dir + "error_analysis_" + filename,
                       {"lambdas", "norm_K1_fRG", "norm_K1_BSE", "K1_diff_rel", "K1_diff_abs",
                        "norm_K2_fRG", "norm_K2_BSE", "K2_diff_rel", "K2_diff_abs",
                        "norm_S_fRG",  "norm_S_SDE", "Sigma_diff_rel", "Sigma_diff_abs"},
                       {actual_gird, norm_K1_fRG, norm_K1_BSE, K1_diff_rel, K1_diff_abs,
                        norm_K2_fRG, norm_K2_BSE, K2_diff_rel, K2_diff_abs,
                        norm_S_fRG,  norm_S_SDE,  Sigma_diff_rel, Sigma_diff_abs});

    }


    write_h5_rvecs(dir + "error_analysis_" + filename,
                   {"lambdas", "norm_K1_fRG", "norm_K1_BSE", "K1_diff_rel", "K1_diff_abs",
                                      "norm_K2_fRG", "norm_K2_BSE", "K2_diff_rel", "K2_diff_abs",
                                      "norm_S_fRG",  "norm_S_SDE", "Sigma_diff_rel", "Sigma_diff_abs"},
                   {actual_gird, norm_K1_fRG, norm_K1_BSE, K1_diff_rel, K1_diff_abs,
                                         norm_K2_fRG, norm_K2_BSE, K2_diff_rel, K2_diff_abs,
                                         norm_S_fRG,  norm_S_SDE,  Sigma_diff_rel, Sigma_diff_abs});

}



#endif //KELDYSH_MFRG_BETHE_SALPETER_H