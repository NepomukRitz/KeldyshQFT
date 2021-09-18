#ifndef KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
#define KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H

#include "../data_structures.h"        // vector class
#include "../state.h"                  // State class
#include "../selfenergy.h"             // Self-energy class
#include "../utilities/util.h"                   // print output

/**
* Function that checks causality of self-energy: Im(Sigma^R)<=0.
*/
template <typename Q>
void check_SE_causality(const SelfEnergy<Q>& selfEnergy) {
    if (KELDYSH) {
        print("Causality check of self-energy: Im(Sigma^R)<=0.", true);

        vec<Q> Sigma = selfEnergy.Sigma;                        // take self-energy
        vec<Q> Sigma_R(&Sigma[0], &Sigma[Sigma.size() / 2]);     // take first half of self-energy (retarded comp.)

        // check if Im(Sigma^R) is positive for every data point
        int cnt = 0;
        double sum = 0.;
        for (int i = 0; i < Sigma_R.size(); ++i) {
            double val = myimag(Sigma_R[i]);
            if (val > 0.) {
                cnt += 1;
                sum += val;
            }
        }
        if (cnt > 0) {
            print("Selfenergy is non-causal: ", true);
            print(cnt, " values of Im(Sigma^R) are positive, with a sum of ", sum, true);
        } else
            print("Selfenergy is causal.", true);
    }
    else {
        print("Causality check of self-energy: Im[Sigma(w)]*w<=0.", true);

        vec<Q> Sigma = selfEnergy.Sigma;                        // take self-energy

        // check if Im(Sigma^R) is positive for every data point
        int cnt = 0;
        double sum = 0.;
        for (int i = 1; i < nFER-1; ++i) {

            double val = myimag(Sigma[i]) * sign(selfEnergy.frequencies.ws[i]);

            if (val > 0.) {
                //cout << "i: " << i << "\t for w = " << selfEnergy.frequencies.ws[i] << "; \t Sigma[i] = " << Sigma[i] << "\n";
                cnt += 1;
                sum += val;
            }
        }
        if (cnt > 0) {
            print("Im[Selfenergy] is not negative for positive w (vice versa): ", true);
            print(cnt, " values of Im(Sigma) have the wrong sign, with a sum of ", sum, true);
        } else
            print("Selfenergy has the right sign.", true);
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
 * Function that checks FDTs for self-energy and K1 in all channels for given input state: Re(Sigma^K)=0, Re(K1r^K)=0.
 * If verbose is true, maximum values of Re(Sigma^K) and Re(K1r^K) are always printed. If verbose is false (default),
 * output is only printed if checks fail.
 */
template <typename Q>
void check_FDTs(const State<Q>& state, bool verbose=false) {
    if (verbose)
        print("Check of FDTs for self-energy and K1: Re(Sigma^K)=0, Re(K1r^K)=0.", true);

    const double EPS = std::numeric_limits<double>::epsilon(); // double precision used as error estimate

    /** 1st check: real part of Keldysh component of the selfenergy has to be zero */

    vec<Q> Sigma = state.selfenergy.Sigma;                          // take self-energy
    vec<Q> Sigma_K (&Sigma[Sigma.size()/2], &Sigma[Sigma.size()]);  // take second half of self-energy (Keldysh comp.)
    double max_Sigma_K = Sigma_K.real().max_norm();                 // take max. value of Re(Sigma^K)

    if (verbose) {
        print("Maximal value of Re(Sigma^K): ", max_Sigma_K, false);
        if (max_Sigma_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(Sigma^K) is non-zero!", true);
    }
    else {
        if (max_Sigma_K > 10 * EPS)
            print("Maximal value of Re(Sigma^K): ", max_Sigma_K,
                  "  --> CHECK FAILED: Re(Sigma^K) is non-zero!", true);
    }

    /** 2nd check: real part of Keldysh component of K1 in all channels has to be zero */

    // take K1 vertices in all channels
    vec<Q> K1a = state.vertex[0].avertex().get_K1();
    vec<Q> K1p = state.vertex[0].pvertex().get_K1();
    vec<Q> K1t = state.vertex[0].tvertex().get_K1();
    // take second half of K1 vertices (Keldysh comp.)
    vec<Q> K1a_K (&K1a[K1a.size()/2], &K1a[K1a.size()]);
    vec<Q> K1p_K (&K1p[K1p.size()/2], &K1p[K1p.size()]);
    vec<Q> K1t_K (&K1t[K1t.size()/2], &K1t[K1t.size()]);
    // take max. value of Re(K1r^K)
    double max_K1a_K = K1a_K.real().max_norm();
    double max_K1p_K = K1p_K.real().max_norm();
    double max_K1t_K = K1t_K.real().max_norm();

    if (verbose) {
        print("Maximal value of Re(K1a^K):   ", max_K1a_K, false);
        if (max_K1a_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(K1a^K) is non-zero!", true);

        print("Maximal value of Re(K1p^K):   ", max_K1p_K, false);
        if (max_K1p_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(K1p^K) is non-zero!", true);

        print("Maximal value of Re(K1t^K):   ", max_K1t_K, false);
        if (max_K1t_K < 10 * EPS) print_add("  --> Check passed.", true);
        else print_add("  --> CHECK FAILED: Re(K1t^K) is non-zero!", true);
    }
    else {
        if (max_K1a_K > 10 * EPS)
            print("Maximal value of Re(K1a^K):   ", max_K1a_K,
                  "  --> CHECK FAILED: Re(K1a^K) is non-zero!", true);
        if (max_K1p_K > 10 * EPS)
            print("Maximal value of Re(K1p^K):   ", max_K1p_K,
                  "  --> CHECK FAILED: Re(K1p^K) is non-zero!", true);
        if (max_K1t_K > 10 * EPS)
            print("Maximal value of Re(K1t^K):   ", max_K1a_K,
                  "  --> CHECK FAILED: Re(K1t^K) is non-zero!", true);
    }
}

template <typename Q>
void check_FDTs(const Q& state, bool verbose=false) {}

#endif //KELDYSH_MFRG_TESTING_CAUSALITY_FDT_CHECKS_H
