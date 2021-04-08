//
// Created by nepomuk on 08.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
#define KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H

#include "../parameters.h"
#include "../data_structures.h"
#include "../state.h"
#include "../perturbation_theory.h"
#include "../propagator.h"
#include "../selfenergy.h"
#include "../bubbles.h"

template<typename Q>
class test_PrecalculateBubble{
public:
    Propagator g = Propagator (Lambda_ini, 'g');
    Propagator s = Propagator (Lambda_ini, 's');
#ifdef KELDYSH_FORMALISM
    int number_of_Keldysh_components = 9;
#else
    int number_of_Keldysh_components = 1;
#endif
public:
    test_PrecalculateBubble(){
        State<comp> testing_state (Lambda_ini);
        testing_state.initialize();
        sopt_state(testing_state, Lambda_ini);

        g = Propagator(Lambda_ini, testing_state.selfenergy, 's');
        s = Propagator(Lambda_ini, testing_state.selfenergy, 'g');
    }

    void perform_test(){
        double max_error = iterate_through_channels();
        if (max_error != 0){
            std::cout << "The largest deviation is " << max_error;
        }
        else {
            std::cout << "Perfect agreement!";
        }
    }

    double iterate_through_channels(){ // currently only 'a' channel
        double max_error = 0;
        const double error_undot =  find_largest_deviation_from_bubble(false, 'a');
        const double error_dot =  find_largest_deviation_from_bubble(true, 'a');
        double larger_error = max(abs(error_undot), abs(error_dot));
        if (larger_error > max_error){
            max_error = larger_error;
        }
        }
        return max_error;
    }

    double find_largest_deviation_from_bubble(const bool dot_in, const char channel){
        PrecalculateBubble<Q> Pre_Bubble(g, s, dot_in, channel);
        Bubble Usual_Bubble (g, s, dot_in);
        double largest_deviation = 0;
        for (int iK = 0; iK < number_of_Keldysh_components; ++iK) {
            for (int iw = 0; iw < nBOS; ++iw) {
                const double w = g.selfenergy.frequencies.w[iw];
                for (int ivpp = 0; ivpp < nFER; ++ivpp) {
                    const double vpp = g.selfenergy.frequencies.w[ivpp];
                    for (int i_in = 0; i_in < n_in; ++i_in) {
                        Q deviation = Pre_Bubble.value(iK, w, vpp, i_in)
                                - Usual_Bubble.value(iK, w, vpp, i_in, channel);
                        if (abs(deviation) > largest_deviation) {
                            largest_deviation = abs(deviation);
                        }
                    }
                }
            }
        }
        return largest_deviation;
    }
};

#endif //KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
