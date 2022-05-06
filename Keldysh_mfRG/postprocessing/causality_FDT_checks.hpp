#ifndef KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
#define KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H

#include "../data_structures.hpp"        // vector class
#include "../correlation_functions/state.hpp"                  // State class
#include "../correlation_functions/two_point/selfenergy.hpp"             // Self-energy class
#include "../utilities/util.hpp"                   // print output
#include "../utilities/hdf5_routines.hpp"

/**
* Function that checks causality of self-energy: Im(Sigma^R)<=0.
*/
template <typename Q>
void check_SE_causality(const SelfEnergy<Q>& selfEnergy) {
    if (KELDYSH) {
        utils::print("Causality check of self-energy: Im(Sigma^R)<=0.", true);
        using SE_buffertype = typename SelfEnergy<Q>::buffer_type;
        SE_buffertype Sigma = selfEnergy.Sigma;                        // take self-energy
        vec<Q> Sigma_R(Sigma.get_vec().begin(), Sigma.get_vec().begin() + (nFER));     // take first half of self-energy (retarded comp.)
        assert(Sigma_R.size() == nFER);


        // check if Im(Sigma^R) is positive for every data point
        int cnt = 0;
        double sum = 0.;
        double max = 0.;
        for (unsigned int i = 0; i < Sigma_R.size(); ++i) {
            double val = myimag(Sigma_R[i]);
            if (val > 0.) {
                cnt += 1;
                sum += val;
                max = std::max(max, val);
            }
        }
        selfEnergy.Sigma.initInterpolator();
        if (cnt > 0) {
            utils::print("Selfenergy is non-causal: ", true);
            utils::print(cnt, " values of Im(Sigma^R) are positive, with a maximum of ", max, true);
            utils::print("value at v=0: ", selfEnergy.valsmooth(0, 0., 0), true);
        } else {
            utils::print("Selfenergy is causal.", true);
            utils::print("value at v=0: ", selfEnergy.valsmooth(0, 0., 0), true);
        }

    }
    else {
        utils::print("Causality check of self-energy: Im[Sigma(w)]*w<=0.", true);

        auto Sigma = selfEnergy.Sigma;                        // take self-energy

        // check if Im(Sigma^R) is positive for every data point
        int cnt = 0;
        double sum = 0.;
        for (int i = 0; i < nFER; ++i) {

            double val = myimag(Sigma.acc(i)) * sign(selfEnergy.Sigma.frequencies.  primary_grid.get_frequency(i));

            if (val > 0.) {
                //cout << "i: " << i << "\t for w = " << selfenergy.Sigma.frequencies.  primary_grid.get_frequency(i) << "; \t Sigma[i] = " << Sigma[i] << "\n";
                cnt += 1;
                sum += val;
            }
        }
        if (cnt > 0) {
            utils::print("Im[Selfenergy] is not negative for positive w (vice versa): ", true);
            utils::print(cnt, " values of Im(Sigma) have the wrong sign, with a sum of ", sum, true);
        } else
            utils::print("Selfenergy has the right sign.", true);
    }
}

// wrapper for the function above, taking a State instead of a SelfEnergy
template <typename Q>
void check_SE_causality(const State<Q>& state) {
    check_SE_causality(state.selfenergy);
}

template <typename Q>
void check_SE_causality(const Q& selfEnergy) {}

/**
 * Function that checks FDTs for self-energy and K1 in all channels for given input state:
 */
template <typename Q>
void check_FDTs_selfenergy(const SelfEnergy<Q>& selfenergy, const bool verbose) {
    std::array<size_t,SE_config.rank> dims_K = SE_config.dims;
    dims_K[my_defs::SE::keldysh] = 1;
    using buffer_t = multidimensional::multiarray<Q,SE_config.rank>;
    buffer_t SE_K_from_data(dims_K);
    buffer_t SE_K_from_FDTs(dims_K);

    const size_t size_flat = getFlatSize<SE_config.rank>(dims_K);
    for (unsigned int i = 0; i < size_flat; i++) {
        std::array<my_index_t,SE_config.rank> idx_R;
        getMultIndex<SE_config.rank>(idx_R,i,dims_K);

        std::array<my_index_t,SE_config.rank> idx_K = idx_R;
        idx_K[my_defs::SE::keldysh] = 1;

        const Q val_K_from_data = selfenergy.Sigma.val(idx_K);
        const Q val_R_from_data = selfenergy.Sigma.val(idx_R);

        // Compute Keldysh component from FDTs:
        int iv = idx_K[my_defs::SE::nu];
        double v;
        selfenergy.Sigma.frequencies.get_freqs_w(v, iv);

        const Q val_K_from_FDTs = glb_i * 2. * Eff_fac(v) * myimag(val_R_from_data);

        SE_K_from_data.at(idx_R) = val_K_from_data;
        SE_K_from_FDTs.at(idx_R) = val_K_from_FDTs;
    }

    buffer_t difference = SE_K_from_data - SE_K_from_FDTs;
    double max_deviation = difference.max_norm();
    double max_norm_Kdat = SE_K_from_data.max_norm();
    if (verbose) {
        utils::print("Check FDT for selfenergy: \n");
        utils::print("Maximal deviation from FDTs (absolute): \t", max_deviation, "\n");
        utils::print("Maximal absolute value of Sigma^K: \t", max_norm_Kdat, "\n");
    }
}


template <typename Q>
void check_FDTs_K1(const rvert<Q>& rvert4K1, const bool verbose) {
    if constexpr(!CONTOUR_BASIS) {
        std::array<size_t, K1at_config.rank> dims_K = rvert4K1.K1.get_dims();
        dims_K[my_defs::K1::keldysh] = 1;
        using buffer_t = multidimensional::multiarray<Q, K1at_config.rank>;
        buffer_t Im_K1_R_from_data(dims_K);
        buffer_t Im_K1_R_from_FDTs(dims_K);

        const size_t size_flat = getFlatSize<K1at_config.rank>(dims_K);
        for (unsigned int i = 0; i < size_flat; i++) {
            std::array<my_index_t, K1at_config.rank> idx_R;
            getMultIndex<K1at_config.rank>(idx_R, i, dims_K);

            std::array<my_index_t, K1at_config.rank> idx_K = idx_R;
            idx_K[my_defs::K1::keldysh] = 1;

            const Q val_K_from_data = rvert4K1.K1.val(idx_K);
            const Q val_Im_R_from_data = myimag(rvert4K1.K1.val(idx_R));

            // Compute Keldysh component from FDTs:
            int iv = idx_K[my_defs::K1::omega];
            double v;
            rvert4K1.K1.frequencies.get_freqs_w(v, iv);

            const Q val_Im_R_from_FDTs = glb_i * 0.5 * Eff_fac(v) * val_K_from_data;

            Im_K1_R_from_data.at(idx_R) = val_Im_R_from_data;
            Im_K1_R_from_FDTs.at(idx_R) = val_Im_R_from_FDTs;
        }

        buffer_t difference = Im_K1_R_from_data - Im_K1_R_from_FDTs;
        double max_deviation = difference.max_norm();
        double max_norm_Kdat = Im_K1_R_from_data.max_norm();
        if (verbose) {
            utils::print("Check FDT for K1", rvert4K1.channel, ": \n");
            utils::print("Maximal deviation from FDTs (absolute): \t", max_deviation, "\n");
            utils::print("Maximal absolute value of Im(K1^R): \t", max_norm_Kdat, "\n");
        }
    }
}


template <typename Q>
void check_FDTs(const State<Q>& state, bool verbose) {
    check_FDTs_selfenergy(state.selfenergy, verbose);
    for (char ch : {'a', 'p', 't'}) {
        check_FDTs_K1(state.vertex.get_rvertex(ch), verbose);
    }
}


/**
 * Function that computes vertex components with the help of fluctuation-dissipation relations.
 * Components obtainable by FDTs:
 * K1r[1] (with r = a/p/t)
 * K2a, K2a[2], K2a[3]
 * K2p[0], K2p[1], K2p[3]
 * K2t[0], K2t[2], K2t[3]
 * K3r[0], K3r[1] (with r = a/p/t)
 *
 * These components are written into state_out.
 * The other components are copied to state_out.
 * CAUTION: For the K2-class the FDTs involve numerically diverging prefactors for zero bosonic frequency
 */
template <typename Q, char channel>
void compute_components_through_FDTs(fullvert<Q>& vertex_out, const fullvert<Q>& vertex_in, const fullvert<Q>& vertex_half2_in) {
    if constexpr(CONTOUR_BASIS) {
        vertex_out = vertex_in;
    } else {
        /// FDTs for K1
        if constexpr(channel == 'a') {
#pragma omp parallel for collapse(3)
            for (int itw = 0; itw < nw1; itw++) {
                for (int itspin = 0; itspin < n_in; itspin++) {
                    for (int itin = 0; itin < n_in; itin++) {
                        my_defs::K1::index_type idx;
                        idx[my_defs::K1::keldysh] = 1;
                        idx[my_defs::K1::spin] = itspin;
                        idx[my_defs::K1::omega] = itw;
                        idx[my_defs::K1::internal] = itin;

                        double w, N1;
                        Q G1, G2;

                        vertex_in.avertex.K1.frequencies.get_freqs_w(w, itw);
                        if (std::abs(w) > glb_T * 25.) {
                            N1 = 1. / Fermi_fac(w, glb_mu);
                            VertexInput inputG1(1, itspin, w, 0, 0, itin, 'a');  // advanced component
                            G1 = vertex_in.avertex.template valsmooth<k1>(inputG1, vertex_in.tvertex,
                                                                          vertex_half2_in.avertex,
                                                                          vertex_half2_in.tvertex);
                            G2 = N1 * (myconj(G1) - G1);
                            vertex_out.avertex.K1.setvert(G2, idx);
                        }
                    }
                }

            }
        } else if constexpr(channel == 'p') {
#pragma omp parallel for collapse(3)
            for (int itw = 0; itw < nw1; itw++) {
                for (int itspin = 0; itspin < n_in; itspin++) {
                    for (int itin = 0; itin < n_in; itin++) {
                        my_defs::K1::index_type idx;
                        idx[my_defs::K1::keldysh] = 1;
                        idx[my_defs::K1::spin] = itspin;
                        idx[my_defs::K1::omega] = itw;
                        idx[my_defs::K1::internal] = itin;

                        double w, N1;
                        Q G1, G2;

                        vertex_in.pvertex.K1.frequencies.get_freqs_w(w, itw);
                        if (std::abs(w) > glb_T * 25.) {
                            N1 = 1. / Fermi_fac(w, glb_mu);
                            VertexInput inputG1(1, itspin, w, 0, 0, itin, 'a');  // advanced component
                            G1 = vertex_in.pvertex.template valsmooth<k1>(inputG1, vertex_in.pvertex,
                                                                          vertex_half2_in.pvertex,
                                                                          vertex_half2_in.pvertex);
                            G2 = N1 * (myconj(G1) - G1);
                            vertex_out.pvertex.K1.setvert(G2, idx);

                        }
                    }
                }

            }
        } else {
            static_assert(channel == 't', "Invalid channel");
#pragma omp parallel for collapse(3)
            for (int itw = 0; itw < nw1; itw++) {
                for (int itspin = 0; itspin < n_in; itspin++) {
                    for (int itin = 0; itin < n_in; itin++) {
                        my_defs::K1::index_type idx;
                        idx[my_defs::K1::keldysh] = 1;
                        idx[my_defs::K1::spin] = itspin;
                        idx[my_defs::K1::omega] = itw;
                        idx[my_defs::K1::internal] = itin;

                        double w, N1;
                        Q G1, G2;

                        vertex_in.tvertex.K1.frequencies.get_freqs_w(w, itw);
                        if (std::abs(w) > glb_T * 25.) {
                            N1 = 1. / Fermi_fac(w, glb_mu);
                            VertexInput inputG1(1, itspin, w, 0, 0, itin, 'a');  // advanced component
                            G1 = vertex_in.tvertex.template valsmooth<k1>(inputG1, vertex_in.avertex,
                                                                          vertex_half2_in.tvertex,
                                                                          vertex_half2_in.avertex);
                            G2 = N1 * (myconj(G1) - G1);
                            vertex_out.tvertex.K1.setvert(G2, idx);

                        }
                    }
                }

            }

        }


#if MAX_DIAG_CLASS > 1 /// FDTs for K2 objects:
        if constexpr(channel == 'a') {
#pragma omp parallel for collapse(4)
            for (int itw = 0; itw < nw2; itw++) {
                for (int itv = 0; itv < nv2; itv++) {
                    for (int itspin = 0; itspin < n_in; itspin++) {
                        for (int itin = 0; itin < n_in; itin++) {
                            my_defs::K2::index_type idx;
                            //idx[my_defs::K2::keldysh] = ;
                            idx[my_defs::K2::spin] = itspin;
                            idx[my_defs::K2::omega] = itw;
                            idx[my_defs::K2::nu] = itv;
                            idx[my_defs::K2::internal] = itin;

                            double w, v, N1, N2, N3;
                            Q G1, G2, G3, G12, G23, G123;

                            vertex_in.avertex.K2.frequencies.get_freqs_w(w, v, itw, itv);
                            if (std::abs(w) > glb_T) {
                                N1 = Fermi_fac(-v - w / 2, glb_mu);
                                N2 = Fermi_fac(v - w / 2, glb_mu);
                                N3 = 1. / Fermi_fac(w, glb_mu);
                                VertexInput inputG1(8, itspin, w, v, 0, itin, 'a');
                                VertexInput inputG2(1, itspin, w, v, 0, itin, 'a');
                                VertexInput inputG3(11, itspin, w, v, 0, itin, 'a');
                                G1 = vertex_in.avertex.template valsmooth<k2>(inputG1, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);
                                G2 = vertex_in.avertex.template valsmooth<k2>(inputG2, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);
                                G3 = vertex_in.avertex.template valsmooth<k2>(inputG3, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);

                                G123 = myconj(G1 + G2 + G3)
                                       + myreal(G1 * N2 * N3
                                                + N1 * G2 * N3
                                                + N1 * N2 * G3) * 2.;
                                G12 = N2 * (myconj(G3) - G1) + N1 * (myconj(G3) - G2);
                                //G13 = N3 * (myconj(G2) - G1) + N1 * (myconj(G2) - G3);
                                G23 = N3 * (myconj(G1) - G2) + N2 * (myconj(G1) - G3);

                                idx[my_defs::K2::keldysh] = 0;
                                vertex_out.avertex.K2.setvert(G12, idx);
                                idx[my_defs::K2::keldysh] = 2;
                                vertex_out.avertex.K2.setvert(G123, idx);
                                idx[my_defs::K2::keldysh] = 3;
                                vertex_out.avertex.K2.setvert(G23, idx);

                            }

                        }
                    }
                }
            }
        } else if constexpr (channel == 'p') {
#pragma omp parallel for collapse(4)
            for (int itw = 0; itw < nw2; itw++) {
                for (int itv = 0; itv < nv2; itv++) {
                    for (int itspin = 0; itspin < n_in; itspin++) {
                        for (int itin = 0; itin < n_in; itin++) {
                            my_defs::K2::index_type idx;
                            //idx[my_defs::K2::keldysh] = ;
                            idx[my_defs::K2::spin] = itspin;
                            idx[my_defs::K2::omega] = itw;
                            idx[my_defs::K2::nu] = itv;
                            idx[my_defs::K2::internal] = itin;

                            double w, v, N1, N2, N3;
                            Q G1, G2, G3, G12, G13, G123;

                            vertex_in.pvertex.K2.frequencies.get_freqs_w(w, v, itw, itv);
                            if (std::abs(w) > glb_T * 1.) {//
                                N1 = Fermi_fac(v + w / 2, glb_mu);
                                N2 = Fermi_fac(-v + w / 2, glb_mu);
                                N3 = 1. / Fermi_fac(-w, glb_mu);
                                VertexInput inputG1(4, itspin, w, v, 0, itin, 'p');
                                VertexInput inputG2(8, itspin, w, v, 0, itin, 'p');
                                VertexInput inputG3(13, itspin, w, v, 0, itin, 'p');
                                G1 = vertex_in.pvertex.template valsmooth<k2>(inputG1, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);
                                G2 = vertex_in.pvertex.template valsmooth<k2>(inputG2, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);
                                G3 = vertex_in.pvertex.template valsmooth<k2>(inputG3, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);

                                G123 = myconj(G1 + G2 + G3)
                                       + myreal(G1 * N2 * N3
                                                + N1 * G2 * N3
                                                + N1 * N2 * G3) * 2.;
                                G12 = N2 * (myconj(G3) - G1) + N1 * (myconj(G3) - G2);
                                G13 = N3 * (myconj(G2) - G1) + N1 * (myconj(G2) - G3);
                                //G23 = N3 * (myconj(G1) - G2) + N2 * (myconj(G1) - G3);
                                idx[my_defs::K2::keldysh] = 0;
                                vertex_out.pvertex.K2.setvert(G12, idx);
                                idx[my_defs::K2::keldysh] = 1;
                                vertex_out.pvertex.K2.setvert(G123, idx);
                                idx[my_defs::K2::keldysh] = 3;
                                vertex_out.pvertex.K2.setvert(G13, idx);

                            }

                        }
                    }
                }
            }

        } else {
            static_assert(channel == 't', "Invalid channel");
#pragma omp parallel for collapse(4)
            for (int itw = 0; itw < nw2; itw++) {
                for (int itv = 0; itv < nv2; itv++) {
                    for (int itspin = 0; itspin < n_in; itspin++) {
                        for (int itin = 0; itin < n_in; itin++) {
                            my_defs::K2::index_type idx;
                            idx[my_defs::K2::spin] = itspin;
                            idx[my_defs::K2::omega] = itw;
                            idx[my_defs::K2::nu] = itv;
                            idx[my_defs::K2::internal] = itin;

                            double w, v, N1, N2, N3;
                            Q G1, G2, G3, G12, G23, G123;

                            vertex_in.tvertex.K2.frequencies.get_freqs_w(w, v, itw, itv);
                            if (std::abs(w) > glb_T * 1.) {
                                N1 = Fermi_fac(-v - w / 2, glb_mu);
                                N2 = Fermi_fac(v - w / 2, glb_mu);
                                N3 = 1. / Fermi_fac(w, glb_mu);
                                VertexInput inputG1(4, itspin, w, v, 0, itin, 't');
                                VertexInput inputG2(1, itspin, w, v, 0, itin, 't');
                                VertexInput inputG3(7, itspin, w, v, 0, itin, 't');
                                G1 = vertex_in.tvertex.template valsmooth<k2>(inputG1, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);
                                G2 = vertex_in.tvertex.template valsmooth<k2>(inputG2, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);
                                G3 = vertex_in.tvertex.template valsmooth<k2>(inputG3, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);

                                G123 = myconj(G1 + G2 + G3)
                                       + myreal(G1 * N2 * N3
                                                + N1 * G2 * N3
                                                + N1 * N2 * G3) * 2.;
                                G12 = N2 * (myconj(G3) - G1) + N1 * (myconj(G3) - G2);
                                //G13 = N3 * (myconj(G2) - G1) + N1 * (myconj(G2) - G3);
                                G23 = N3 * (myconj(G1) - G2) + N2 * (myconj(G1) - G3);
                                idx[my_defs::K2::keldysh] = 0;
                                vertex_out.tvertex.K2.setvert(G12, idx);
                                idx[my_defs::K2::keldysh] = 2;
                                vertex_out.tvertex.K2.setvert(G123, idx);
                                idx[my_defs::K2::keldysh] = 3;
                                vertex_out.tvertex.K2.setvert(G23, idx);
                            }
                        }
                    }
                }
            }
        }


#endif

#if MAX_DIAG_CLASS > 2
        if constexpr(channel == 'a') {
#pragma omp parallel for collapse(5)
            for (int itw = 0; itw < nw3; itw++) {
                for (int itv = 0; itv < nv3; itv++) {
                    for (int itvp = 0; itvp < (GRID != 2 ? nv3 : (nv3 - 1) / 2 + 1); itvp++) {
                        for (int itspin = 0; itspin < n_spin; itspin++) {
                            for (int itin = 0; itin < n_in; itin++) {
                                my_defs::K3::index_type idx;
                                idx[my_defs::K3::spin] = itspin;
                                idx[my_defs::K3::omega] = itw;
                                idx[my_defs::K3::nu] = itv;
                                idx[my_defs::K3::nup] = itvp;
                                idx[my_defs::K3::internal] = itin;

                                double w, v, vp, N1, N2, N3, N4;
                                Q G1, G2, G3, G4, G12, G13, G14, G23, G24, G34, G1234, G123;

                                vertex_in.avertex.K3.frequencies.get_freqs_w(w, v, vp, itw, itv, itvp);
                                N1 = Fermi_fac(v - w / 2, glb_mu);
                                N2 = Fermi_fac(vp + w / 2, glb_mu);
                                N3 = Fermi_fac(-vp + w / 2, glb_mu);
                                N4 = Fermi_fac(-v - w / 2, glb_mu);
                                VertexInput inputG1(7, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG2(11, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG3(13, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG4(14, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG12(3, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG13(5, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG14(6, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG23(9, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG24(10, itspin, w, v, vp, itin, 'a');
                                VertexInput inputG34(12, itspin, w, v, vp, itin, 'a');
                                G1 = vertex_in.avertex.template valsmooth<k3>(inputG1, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);
                                G2 = vertex_in.avertex.template valsmooth<k3>(inputG2, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);
                                G3 = vertex_in.avertex.template valsmooth<k3>(inputG3, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);
                                G4 = vertex_in.avertex.template valsmooth<k3>(inputG4, vertex_in.tvertex,
                                                                              vertex_half2_in.avertex,
                                                                              vertex_half2_in.tvertex);
                                G12 = vertex_in.avertex.template valsmooth<k3>(inputG12, vertex_in.tvertex,
                                                                               vertex_half2_in.avertex,
                                                                               vertex_half2_in.tvertex);
                                G13 = vertex_in.avertex.template valsmooth<k3>(inputG13, vertex_in.tvertex,
                                                                               vertex_half2_in.avertex,
                                                                               vertex_half2_in.tvertex);
                                G14 = vertex_in.avertex.template valsmooth<k3>(inputG14, vertex_in.tvertex,
                                                                               vertex_half2_in.avertex,
                                                                               vertex_half2_in.tvertex);
                                G23 = vertex_in.avertex.template valsmooth<k3>(inputG23, vertex_in.tvertex,
                                                                               vertex_half2_in.avertex,
                                                                               vertex_half2_in.tvertex);
                                G24 = vertex_in.avertex.template valsmooth<k3>(inputG24, vertex_in.tvertex,
                                                                               vertex_half2_in.avertex,
                                                                               vertex_half2_in.tvertex);
                                G34 = vertex_in.avertex.template valsmooth<k3>(inputG34, vertex_in.tvertex,
                                                                               vertex_half2_in.avertex,
                                                                               vertex_half2_in.tvertex);

                                G123 = (1 + N1 * N2 + N1 * N3 + N2 * N3) * myconj(G4)
                                       - (G1 * N2 * N3
                                          + N1 * G2 * N3
                                          + N1 * N2 * G3
                                          + N1 * G23
                                          + G12 * N3
                                          + G13 * N2);
                                G1234 = G1 * N2 * N3 * N4 * 2.
                                        + N1 * G2 * N3 * N4 * 2.
                                        + N1 * N2 * G3 * N4 * 2.
                                        + N1 * N2 * N3 * G4 * 2.
                                        + (N2 * N3 * N4 + N2 + N3 + N4) * myconj(G1)
                                        + (N1 * N3 * N4 + N1 + N3 + N4) * myconj(G2)
                                        + (N1 * N2 * N4 + N1 + N2 + N4) * myconj(G3)
                                        + (N1 * N2 * N3 + N1 + N2 + N3) * myconj(G4)
                                        + N3 * N4 * G12
                                        + N2 * N4 * G13
                                        + N2 * N3 * G14
                                        + N1 * N4 * G23
                                        + N1 * N3 * G24
                                        + N1 * N2 * G34;
                                idx[my_defs::K3::keldysh] = 0;
                                vertex_out.avertex.K3.setvert(G1234, idx);
                                idx[my_defs::K3::keldysh] = 1;
                                vertex_out.avertex.K3.setvert(G123, idx);
                            }
                        }
                    }
                }
            }
        } else if constexpr (channel == 'p') {
#pragma omp parallel for collapse(5)
            for (int itw = 0; itw < nw3; itw++) {
                for (int itv = 0; itv < nv3; itv++) {
                    for (int itvp = 0; itvp < (GRID != 2 ? nv3 : (nv3 - 1) / 2 + 1); itvp++) {
                        for (int itspin = 0; itspin < n_spin; itspin++) {
                            for (int itin = 0; itin < n_in; itin++) {
                                my_defs::K3::index_type idx;
                                idx[my_defs::K3::spin] = itspin;
                                idx[my_defs::K3::omega] = itw;
                                idx[my_defs::K3::nu] = itv;
                                idx[my_defs::K3::nup] = itvp;
                                idx[my_defs::K3::internal] = itin;

                                double w, v, vp, N1, N2, N3, N4;
                                Q G1, G2, G3, G4, G12, G13, G14, G23, G24, G34, G1234, G123;

                                vertex_in.pvertex.K3.frequencies.get_freqs_w(w, v, vp, itw, itv, itvp);
                                N1 = Fermi_fac(v + w / 2, glb_mu);
                                N2 = Fermi_fac(-v + w / 2, glb_mu);
                                N3 = Fermi_fac(-vp - w / 2, glb_mu);
                                N4 = Fermi_fac(vp - w / 2, glb_mu);
                                VertexInput inputG1(7, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG2(11, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG3(13, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG4(14, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG12(3, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG13(5, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG14(6, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG23(9, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG24(10, itspin, w, v, vp, itin, 'p');
                                VertexInput inputG34(12, itspin, w, v, vp, itin, 'p');
                                G1 = vertex_in.pvertex.template valsmooth<k3>(inputG1, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);
                                G2 = vertex_in.pvertex.template valsmooth<k3>(inputG2, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);
                                G3 = vertex_in.pvertex.template valsmooth<k3>(inputG3, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);
                                G4 = vertex_in.pvertex.template valsmooth<k3>(inputG4, vertex_in.pvertex,
                                                                              vertex_half2_in.pvertex,
                                                                              vertex_half2_in.pvertex);
                                G12 = vertex_in.pvertex.template valsmooth<k3>(inputG12, vertex_in.pvertex,
                                                                               vertex_half2_in.pvertex,
                                                                               vertex_half2_in.pvertex);
                                G13 = vertex_in.pvertex.template valsmooth<k3>(inputG13, vertex_in.pvertex,
                                                                               vertex_half2_in.pvertex,
                                                                               vertex_half2_in.pvertex);
                                G14 = vertex_in.pvertex.template valsmooth<k3>(inputG14, vertex_in.pvertex,
                                                                               vertex_half2_in.pvertex,
                                                                               vertex_half2_in.pvertex);
                                G23 = vertex_in.pvertex.template valsmooth<k3>(inputG23, vertex_in.pvertex,
                                                                               vertex_half2_in.pvertex,
                                                                               vertex_half2_in.pvertex);
                                G24 = vertex_in.pvertex.template valsmooth<k3>(inputG24, vertex_in.pvertex,
                                                                               vertex_half2_in.pvertex,
                                                                               vertex_half2_in.pvertex);
                                G34 = vertex_in.pvertex.template valsmooth<k3>(inputG34, vertex_in.pvertex,
                                                                               vertex_half2_in.pvertex,
                                                                               vertex_half2_in.pvertex);

                                G123 = (1 + N1 * N2 + N1 * N3 + N2 * N3) * myconj(G4)
                                       - (G1 * N2 * N3
                                          + N1 * G2 * N3
                                          + N1 * N2 * G3
                                          + N1 * G23
                                          + G12 * N3
                                          + G13 * N2);
                                G1234 = G1 * N2 * N3 * N4 * 2.
                                        + N1 * G2 * N3 * N4 * 2.
                                        + N1 * N2 * G3 * N4 * 2.
                                        + N1 * N2 * N3 * G4 * 2.
                                        + (N2 * N3 * N4 + N2 + N3 + N4) * myconj(G1)
                                        + (N1 * N3 * N4 + N1 + N3 + N4) * myconj(G2)
                                        + (N1 * N2 * N4 + N1 + N2 + N4) * myconj(G3)
                                        + (N1 * N2 * N3 + N1 + N2 + N3) * myconj(G4)
                                        + N3 * N4 * G12
                                        + N2 * N4 * G13
                                        + N2 * N3 * G14
                                        + N1 * N4 * G23
                                        + N1 * N3 * G24
                                        + N1 * N2 * G34;

                                idx[my_defs::K3::keldysh] = 0;
                                vertex_out.pvertex.K3.setvert(G1234, itspin, itw, itv, itvp, 0, itin);
                                idx[my_defs::K3::keldysh] = 1;
                                vertex_out.pvertex.K3.setvert(G123, itspin, itw, itv, itvp, 1, itin);

                            }
                        }
                    }
                }
            }
        } else {
            static_assert(channel == 't', "Invalid channel");

#pragma omp parallel for collapse(5)
            for (int itw = 0; itw < nw3; itw++) {
                for (int itv = 0; itv < nv3; itv++) {
                    for (int itvp = 0; itvp < (GRID != 2 ? nv3 : (nv3 - 1) / 2 + 1); itvp++) {
                        for (int itspin = 0; itspin < n_spin; itspin++) {
                            for (int itin = 0; itin < n_in; itin++) {
                                my_defs::K3::index_type idx;
                                idx[my_defs::K3::spin] = itspin;
                                idx[my_defs::K3::omega] = itw;
                                idx[my_defs::K3::nu] = itv;
                                idx[my_defs::K3::nup] = itvp;
                                idx[my_defs::K3::internal] = itin;

                                double w, v, vp, N1, N2, N3, N4;
                                Q G1, G2, G3, G4, G12, G13, G14, G23, G24, G34, G1234, G123;

                                vertex_in.tvertex.K3.frequencies.get_freqs_w(w, v, vp, itw, itv, itvp);
                                N1 = Fermi_fac(vp + w / 2, glb_mu);
                                N2 = Fermi_fac(v - w / 2, glb_mu);
                                N3 = Fermi_fac(-vp + w / 2, glb_mu);
                                N4 = Fermi_fac(-v - w / 2, glb_mu);

                                VertexInput inputG1(7, itspin, w, v, vp, itin, 't');
                                VertexInput inputG2(11, itspin, w, v, vp, itin, 't');
                                VertexInput inputG3(13, itspin, w, v, vp, itin, 't');
                                VertexInput inputG4(14, itspin, w, v, vp, itin, 't');
                                VertexInput inputG12(3, itspin, w, v, vp, itin, 't');
                                VertexInput inputG13(5, itspin, w, v, vp, itin, 't');
                                VertexInput inputG14(6, itspin, w, v, vp, itin, 't');
                                VertexInput inputG23(9, itspin, w, v, vp, itin, 't');
                                VertexInput inputG24(10, itspin, w, v, vp, itin, 't');
                                VertexInput inputG34(12, itspin, w, v, vp, itin, 't');
                                G1 = vertex_in.tvertex.template valsmooth<k3>(inputG1, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);
                                G2 = vertex_in.tvertex.template valsmooth<k3>(inputG2, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);
                                G3 = vertex_in.tvertex.template valsmooth<k3>(inputG3, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);
                                G4 = vertex_in.tvertex.template valsmooth<k3>(inputG4, vertex_in.avertex,
                                                                              vertex_half2_in.tvertex,
                                                                              vertex_half2_in.avertex);
                                G12 = vertex_in.tvertex.template valsmooth<k3>(inputG12, vertex_in.avertex,
                                                                               vertex_half2_in.tvertex,
                                                                               vertex_half2_in.avertex);
                                G13 = vertex_in.tvertex.template valsmooth<k3>(inputG13, vertex_in.avertex,
                                                                               vertex_half2_in.tvertex,
                                                                               vertex_half2_in.avertex);
                                G14 = vertex_in.tvertex.template valsmooth<k3>(inputG14, vertex_in.avertex,
                                                                               vertex_half2_in.tvertex,
                                                                               vertex_half2_in.avertex);
                                G23 = vertex_in.tvertex.template valsmooth<k3>(inputG23, vertex_in.avertex,
                                                                               vertex_half2_in.tvertex,
                                                                               vertex_half2_in.avertex);
                                G24 = vertex_in.tvertex.template valsmooth<k3>(inputG24, vertex_in.avertex,
                                                                               vertex_half2_in.tvertex,
                                                                               vertex_half2_in.avertex);
                                G34 = vertex_in.tvertex.template valsmooth<k3>(inputG34, vertex_in.avertex,
                                                                               vertex_half2_in.tvertex,
                                                                               vertex_half2_in.avertex);

                                G123 = (1 + N1 * N2 + N1 * N3 + N2 * N3) * myconj(G4)
                                       - (G1 * N2 * N3
                                          + N1 * G2 * N3
                                          + N1 * N2 * G3
                                          + N1 * G23
                                          + G12 * N3
                                          + G13 * N2);
                                G1234 = G1 * N2 * N3 * N4 * 2.
                                        + N1 * G2 * N3 * N4 * 2.
                                        + N1 * N2 * G3 * N4 * 2.
                                        + N1 * N2 * N3 * G4 * 2.
                                        + (N2 * N3 * N4 + N2 + N3 + N4) * myconj(G1)
                                        + (N1 * N3 * N4 + N1 + N3 + N4) * myconj(G2)
                                        + (N1 * N2 * N4 + N1 + N2 + N4) * myconj(G3)
                                        + (N1 * N2 * N3 + N1 + N2 + N3) * myconj(G4)
                                        + N3 * N4 * G12
                                        + N2 * N4 * G13
                                        + N2 * N3 * G14
                                        + N1 * N4 * G23
                                        + N1 * N3 * G24
                                        + N1 * N2 * G34;
                                idx[my_defs::K3::keldysh] = 0;
                                vertex_out.tvertex.K3.setvert(G1234, itspin, itw, itv, itvp, 0, itin);
                                idx[my_defs::K3::keldysh] = 1;
                                vertex_out.tvertex.K3.setvert(G123, itspin, itw, itv, itvp, 1, itin);

                            }
                        }
                    }
                }
            }
        }
#endif
    }
}

/*
 * Wrapper for above function (here defined for GeneralVertex)
 */
template <typename Q, vertexType symmetry_type>
void compute_components_through_FDTs(GeneralVertex<Q,symmetry_type>& vertex_out, const GeneralVertex<Q,symmetry_type>& vertex_in) {
    vertex_in.initializeInterpol();
    compute_components_through_FDTs<Q,'a'>(vertex_out.half1(), vertex_in.half1(), vertex_in.half1());
    compute_components_through_FDTs<Q,'p'>(vertex_out.half1(), vertex_in.half1(), vertex_in.half1());
    compute_components_through_FDTs<Q,'t'>(vertex_out.half1(), vertex_in.half1(), vertex_in.half1());
    vertex_in.set_initializedInterpol(false);
}

/*
 *
 */
template <typename Q, vertexType symmetry_type>
void compare_with_FDTs(const GeneralVertex<Q,symmetry_type>& vertex_in, double Lambda, int Lambda_it, std::string filename_prefix, bool write_flag = false, int nLambda = 1) {
    if(KELDYSH) {
        if (CONTOUR_BASIS and not ZERO_T) {} //assert(false);
        if (CONTOUR_BASIS and ZERO_T) {
            utils::print("Checking the FDTs for Lambda_it", Lambda_it, true);

            State<Q> state_diff(Lambda);

            double max_deviation_K1 = 0, max_deviation_K2 = 0, max_deviation_K3 = 0;
            for (char r: {'a', 'p', 't'}) {
                auto K1 = vertex_in.get_rvertex(r).K1;
                for (int iflat = 0; iflat < getFlatSize(K1.get_dims()); iflat++) {
                    my_defs::K1::index_type idx;
                    getMultIndex<rank_K1>(idx, iflat, K1.get_dims());
                    int itK = (int) idx[my_defs::K1::keldysh];
                    my_index_t itw = idx[my_defs::K1::omega];
                    double w;
                    K1.frequencies.get_freqs_w(w, itw);

                    if (is_zero_due_to_FDTs<k1>(itK, w, 0., 0., r)) {
                        Q value = K1.val(idx);
                        state_diff.vertex.get_rvertex(r).K1.setvert(value, idx);
                        max_deviation_K1 = std::max(max_deviation_K1, std::abs(value));

                    }
                }
            }

            if (MAX_DIAG_CLASS > 1) {
                for (char r: {'a', 'p', 't'}) {
                    auto K2 = vertex_in.get_rvertex(r).K2;
                    for (int iflat = 0; iflat < getFlatSize(K2.get_dims()); iflat++) {
                        my_defs::K2::index_type idx;
                        getMultIndex<rank_K2>(idx, iflat, K2.get_dims());
                        int itK = (int) idx[my_defs::K2::keldysh];
                        my_index_t itw = idx[my_defs::K2::omega];
                        my_index_t itv = idx[my_defs::K2::nu];
                        double w, v;
                        K2.frequencies.get_freqs_w(w, v, itw, itv);

                        if (is_zero_due_to_FDTs<k2>(itK, w, v, 0., r)) {
                            Q value = K2.val(idx);
                            max_deviation_K2 = std::max(max_deviation_K2, std::abs(value));
                            state_diff.vertex.get_rvertex(r).K2.setvert(value, idx);
                        }
                    }
                }
            }

            if (MAX_DIAG_CLASS > 2) {
                for (char r: {'a', 'p', 't'}) {
                    auto K3 = vertex_in.get_rvertex(r).K3;
                    for (int iflat = 0; iflat < getFlatSize(K3.get_dims()); iflat++) {
                        my_defs::K3::index_type idx;
                        getMultIndex<rank_K3>(idx, iflat, K3.get_dims());
                        int itK = (int) idx[my_defs::K3::keldysh];
                        my_index_t itw = idx[my_defs::K3::omega];
                        my_index_t itv = idx[my_defs::K3::nu];
                        my_index_t itvp= idx[my_defs::K3::nup];
                        double w, v, vp;
                        K3.frequencies.get_freqs_w(w, v, vp, itw, itv, itvp);

                        if (is_zero_due_to_FDTs<k3>(itK, w, v, vp, r)) {
                            Q value = K3.val(idx);
                            max_deviation_K3 = std::max(max_deviation_K3, std::abs(value));
                            state_diff.vertex.get_rvertex(r).K3.setvert(value, idx);
                        }
                    }
                }
            }

            if (mpi_world_rank() == 0) {
                utils::print("K1: max-norm of deviation = ", false);
                std::cout << max_deviation_K1 << std::scientific << '\n';
                utils::print("K1: relative deviation = ", false);
                std::cout << max_deviation_K1 / vertex_in.half1().norm_K1(0) << std::scientific << '\n';
            }


            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 1) {
                utils::print("K2: max-norm of deviation = ", false);
                std::cout << max_deviation_K2 << std::scientific << '\n';
                utils::print("K2: relative deviation = ", false);
                std::cout << max_deviation_K2 / vertex_in.half1().norm_K2(0) << std::scientific << '\n';
            }


            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 2) {
                utils::print("K3: max-norm of deviation = ", false);
                std::cout << max_deviation_K3 << std::scientific << '\n';
                utils::print("K3: relative deviation = ", false);
                std::cout << max_deviation_K3 / vertex_in.half1().norm_K3(0) << std::scientific << '\n';
            }

            write_state_to_hdf(data_dir + filename_prefix + "_FDTdiff", Lambda, nLambda, state_diff);

        }
        else {
            GeneralVertex<Q, symmetry_type> vertex_out = vertex_in;
            compute_components_through_FDTs(vertex_out, vertex_in);

            utils::print("Checking the FDTs for Lambda_it", Lambda_it, true);
            GeneralVertex<Q, symmetry_type> vertex_diff = vertex_in - vertex_out;
            utils::print("K2: max-norm of deviation = ", false);
            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 1)
                std::cout << vertex_diff.half1().norm_K2(0) << std::scientific << '\n';
            utils::print("K2: relative deviation = ", false);
            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 1)
                std::cout << vertex_diff.half1().norm_K2(0) / vertex_out.half1().norm_K2(0) << std::scientific << '\n';
            utils::print("K2: max-norm = ", false);
            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 1)
                std::cout << vertex_out.half1().norm_K2(0) << std::scientific << '\n';
            utils::print("K3: max-norm of deviation = ", false);
            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 2)
                std::cout << vertex_diff.half1().norm_K3(0) << std::scientific << '\n';
            utils::print("K3: relative deviation = ", false);
            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 2)
                std::cout << vertex_diff.half1().norm_K3(0) / vertex_out.half1().norm_K3(0) << std::scientific << '\n';
            utils::print("K3: max-norm ", false);
            if (mpi_world_rank() == 0 and MAX_DIAG_CLASS > 2)
                std::cout << vertex_out.half1().norm_K3(0) << std::scientific << '\n';
            //

            if (write_flag) {
                SelfEnergy<Q> SE_empty(Lambda);
                Vertex<Q> temp_diff(Lambda);
                temp_diff.half1() = vertex_diff.half1();
                Vertex<Q> temp_out(Lambda);
                temp_out.half1() = vertex_out.half1();
                State<Q> state_out(temp_out, SE_empty, Lambda);
                State<Q> state_diff(temp_diff, SE_empty, Lambda);
                if (Lambda_it == 0) {
                    write_state_to_hdf(data_dir + filename_prefix + "_FDTresult", Lambda, nLambda, state_out);
                    write_state_to_hdf(data_dir + filename_prefix + "_FDTdiff", Lambda, nLambda, state_diff);
                } else {
                    add_state_to_hdf(data_dir + filename_prefix + "_FDTresult", Lambda_it, state_out);
                    add_state_to_hdf(data_dir + filename_prefix + "_FDTdiff", Lambda_it, state_diff);
                }
            }
        }

    }

}


/*
 *
 */
//template <typename Q>
void compare_flow_with_FDTs(const std::string filename, bool write_flag = false);


#endif //KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
