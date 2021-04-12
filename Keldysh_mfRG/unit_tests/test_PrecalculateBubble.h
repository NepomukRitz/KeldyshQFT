//
// Created by nepomuk on 08.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
#define KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H

#include <iostream>
#include <string>
#include "../parameters.h"
#include "../data_structures.h"
#include "../state.h"
#include "../perturbation_theory.h"
#include "../bubbles.h"
#include "../propagator.h"
#include "../write_data2file.h"            // write vectors into hdf5 file




template<typename Q>
class test_PrecalculateBubble{
public:
    Propagator g;
    Propagator s;
#ifdef KELDYSH_FORMALISM
    int number_of_Keldysh_components = 9;
#else
    int number_of_Keldysh_components = 1;
#endif
public:

    test_PrecalculateBubble(): g(Lambda_ini, 'g'), s(Lambda_ini, 's'){
        State<comp> testing_state (Lambda_ini);
        testing_state.initialize();
        sopt_state(testing_state, Lambda_ini);

        g = Propagator(Lambda_ini, testing_state.selfenergy, 's');
        s = Propagator(Lambda_ini, testing_state.selfenergy, 'g');
    }

    void perform_test();
    double iterate_through_channels();
    double find_largest_deviation_from_bubble(bool dot_in, char channel);
    void save_data(string& filename,
                   vec<Q>& ValuesOfPreBubble, vec<Q>& ValuesOfUsualBubble,
                   vec<double>& AbsoluteDeviations);
    string build_filename(bool dot_in, char channel);
    void test_for_zero_value(Q& value, int& number_of_zero_values);

};

template<typename Q> void test_PrecalculateBubble<Q>::perform_test(){
    double max_error = iterate_through_channels();
    if (max_error != 0.){
        string deviation = to_string(max_error);
        std::cout << "The largest deviation is "  << deviation;
    }
    else {
        std::cout << "Perfect agreement!";
    }
}

template<typename Q> double test_PrecalculateBubble<Q>::iterate_through_channels(){ // currently only 'a' channel
    double max_error = 0;
    const double error_undot =  find_largest_deviation_from_bubble(false, 'a');
    const double error_dot =  find_largest_deviation_from_bubble(true, 'a');
    double larger_error = max(abs(error_undot), abs(error_dot));
    if (larger_error > max_error){
        max_error = larger_error;
    }
    return max_error;
}

template<typename Q> double test_PrecalculateBubble<Q>::find_largest_deviation_from_bubble(const bool dot_in, const char channel){
    PrecalculateBubble<Q> Pre_Bubble(g, s, dot_in, channel);
    Bubble Usual_Bubble (g, s, dot_in);

    vec<Q> ValuesOfPreBubble (number_of_Keldysh_components*nBOS*nFER*n_in);
    vec<Q> ValuesOfUsualBubble (number_of_Keldysh_components*nBOS*nFER*n_in);

    vec<double> AbsoluteDeviations (number_of_Keldysh_components*nBOS*nFER*n_in);

    double largest_deviation = 0;
    int number_of_zero_pre_results = 0;

    for (int iK = 0; iK < number_of_Keldysh_components; ++iK) {
        for (int iw = 0; iw < nBOS; ++iw) {
            const double w = g.selfenergy.frequencies.w[iw];
            for (int ivpp = 0; ivpp < nFER; ++ivpp) {
                const double vpp = g.selfenergy.frequencies.w[ivpp];
                for (int i_in = 0; i_in < n_in; ++i_in) {
                    Q Bubble_Value = Usual_Bubble.value(iK, w, vpp, i_in, channel);
                    Q Pre_Bubble_Value = Pre_Bubble.value(iK, w, vpp, i_in);
                    //test_for_zero_value(Pre_Bubble_Value, number_of_zero_pre_results);

                    ValuesOfPreBubble[Pre_Bubble.composite_index(iK, iw, ivpp, i_in)] = Pre_Bubble_Value;
                    ValuesOfUsualBubble[Pre_Bubble.composite_index(iK, iw, ivpp, i_in)] = Bubble_Value;

                    Q deviation = Pre_Bubble_Value - Bubble_Value;
                    AbsoluteDeviations[Pre_Bubble.composite_index(iK, iw, ivpp, i_in)] = abs(deviation);
                    if (abs(deviation) > largest_deviation) {
                        largest_deviation = abs(deviation);
                        // string dev_string = to_string(largest_deviation);
                        // std::cout << dev_string << "\n";
                    }
                }
            }
        }
    }
    std::cout << "Number of zero results encountered for precalculated bubble = "
    + to_string(number_of_zero_pre_results) + "\n";
    string filename = build_filename(dot_in, channel);
    save_data(filename, ValuesOfPreBubble, ValuesOfUsualBubble, AbsoluteDeviations);
    return largest_deviation;
}

template<typename Q>
void
test_PrecalculateBubble<Q>::save_data(string& filename,
                                      vec<Q>& ValuesOfPreBubble,
                                      vec<Q>& ValuesOfUsualBubble,
                                      vec<double>& AbsoluteDeviations) {
write_h5_rvecs(filename,
               {"Frequencies",
                "RealValuesOfPrecalculatedBubble", "RealValuesOfUsualBubble",
                "ImaginaryValuesOfPrecalculatedBubble", "ImaginaryValuesOfUsualBubble",
                "AbsoluteDeviations"},
               {g.selfenergy.frequencies.w,
                ValuesOfPreBubble.real(), ValuesOfUsualBubble.real(),
                ValuesOfPreBubble.imag(), ValuesOfUsualBubble.imag(),
                AbsoluteDeviations});
}

template<typename Q>
string test_PrecalculateBubble<Q>::build_filename(bool dot_in, char channel) {
    string filename = "/scratch-local/Nepomuk.Ritz/testing_data/test_PrecalculateBubble_";
    //string filename = "../../Data/test_PrecalculateBubble_";
    filename += channel;
    if (dot_in) {filename += "_dot";}
    filename += "_nK=" + to_string(number_of_Keldysh_components)
                + "_nBOS=" + to_string(nBOS)
                + "_nFER=" + to_string(nFER)
                + "_n_in=" + to_string(n_in)
                + ".h5";
    return filename;
}

template<typename Q>
void test_PrecalculateBubble<Q>::test_for_zero_value(Q& value, int& number_of_zero_values) {
    if ((value.real() == 0)
        && (value.imag() == 0)){
        std::cout << "WARNING! Precalculated value exactly zero! \n";
        number_of_zero_values += 1;
    }
    else{
        std::cout << "Re Pre-Bubble = "
                     + to_string(value.real()) + "\n";
        std::cout << "Im Pre-Bubble = "
                     + to_string(value.imag()) + "\n";
    }
}


#endif //KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
