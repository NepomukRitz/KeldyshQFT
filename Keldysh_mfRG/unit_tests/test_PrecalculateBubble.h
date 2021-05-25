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
#include "../util.h"




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

        g = Propagator(Lambda_ini, testing_state.selfenergy, 'g');
        s = Propagator(Lambda_ini, testing_state.selfenergy, 's');
        // g = Propagator(Lambda_ini, 'g');
        // s = Propagator(Lambda_ini, 's');
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
    string filename = "/scratch-local/Nepomuk.Ritz/testing_data/true_test_PrecalculateBubble_";
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

template<typename Q>
class Runtime_comparison{
#ifdef KELDYSH_FORMALISM
    int number_of_Keldysh_components = 9;
#else
    int number_of_Keldysh_components = 1;
#endif
    Propagator g;
    Propagator s;
    PrecalculateBubble<Q> Pre_Bubble; // free versions of the bubbles
    Bubble Usual_Bubble;
public:
    Runtime_comparison(): g (Lambda_ini, 'g'), s (Lambda_ini, 's'),
                       Pre_Bubble (g, s, 0, 'a'), Usual_Bubble (g, s, 0){
     State<comp> testing_state (Lambda_ini);
     testing_state.initialize();
     sopt_state(testing_state, Lambda_ini);

     g = Propagator (Lambda_ini, testing_state.selfenergy,'g'); // this is done to obtain the frequency grid
     s = Propagator (Lambda_ini, testing_state.selfenergy,'s');
    }

    void test_runtimes(int max_number_of_iterations);
    double run_iterations(int iterations, bool precalculated);
};

template<typename Q>
void Runtime_comparison<Q>::test_runtimes(int max_number_of_iterations) {
    vec<double> times_usual (max_number_of_iterations);
    vec<double> times_precalculated (max_number_of_iterations);
    for (int t = 0; t < max_number_of_iterations; ++t) {
        times_usual[t] = run_iterations(t, 1);
        times_precalculated[t] = run_iterations(t, 0);
    }
    string filename = "/scratch-local/Nepomuk.Ritz/testing_data/runtime_comparisons.h5";
    write_h5_rvecs(filename,
                   {"runtimes_usual", "runtimes_precalculated"},
                   {times_usual, times_precalculated});
}

template<typename Q>
double Runtime_comparison<Q>::run_iterations(int iterations, bool precalculated) {
    double starting_time = get_time();
    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (int iK = 0; iK < number_of_Keldysh_components; ++iK) {
            for (int iw = 0; iw < nBOS; ++iw) {
                const double w = g.selfenergy.frequencies.w[iw];
                for (int ivpp = 0; ivpp < nFER; ++ivpp) {
                    const double vpp = g.selfenergy.frequencies.w[ivpp];
                    for (int i_in = 0; i_in < n_in; ++i_in) {
                        if (precalculated){
                            Q Pre_Bubble_Value = Pre_Bubble.value(iK, w, vpp, i_in);
                        }
                        else {
                            Q Bubble_Value = Usual_Bubble.value(iK, w, vpp, i_in, 'a');
                        }
                    }
                }
            }
        }
    }
    double end_time = get_time();
    return end_time - starting_time; // time given in milliseconds
}

void test_Bubble_in_Momentum_Space(){
    Propagator g (0, 'g');
    Propagator s (0, 's');

    double starting_time = get_time();
    PrecalculateBubble<comp> DotBubble (g, s, true);
    double end_time = get_time();
    double diff = (end_time - starting_time); // time given in seconds
    std::cout << "Time for differentiated Bubble = " << diff << " s." << "\n";

    starting_time =  get_time();
    PrecalculateBubble<comp> Bubble (g, s, false);
    end_time = get_time();
    diff = (end_time - starting_time);
    std::cout << "Time for undifferentiated Bubble = " << diff << " s." << "\n";

#ifdef KELDYSH_FORMALISM
    vec<comp> g_R (nFER * n_in);
    vec<comp> g_K (nFER * n_in);
    vec<comp> s_R (nFER * n_in);
    vec<comp> s_K (nFER * n_in);
    for (int iv = 0; iv < nFER; ++iv) {
        for (int i_in = 0; i_in < n_in; ++i_in) {
            double v = g.selfenergy.frequencies.w[iv];
            g_R[iv * n_in + i_in] = g.valsmooth(0, v, i_in);
            g_K[iv * n_in + i_in] = g.valsmooth(1, v, i_in);
            s_R[iv * n_in + i_in] = s.valsmooth(0, v, i_in);
            s_K[iv * n_in + i_in] = s.valsmooth(1, v, i_in);
        }
    }

    string filename = "/scratch-local/Nepomuk.Ritz/testing_data/KELDYSH_bubble_in_mom_space_Nq_"
                      + to_string(glb_N_q) + ".h5";
    write_h5_rvecs(filename,
                   {"propagator_frequencies", "bubble_frequencies",
                    "RealValuesOfRetardedPropagator", "ImaginaryValuesOfRetardedPropagator",
                    "RealValuesOfKeldyshPropagator", "ImaginaryValuesOfKeldyshPropagator",
                    "RealValuesOfRetardedSingleScale", "ImaginaryValuesOfRetardedSingleScale",
                    "RealValuesOfKeldyshSingleScale", "ImaginaryValuesOfKeldyshSingleScale",
                    "RealValuesOfBubble", "ImaginaryValuesOfBubble",
                    "RealValuesOfDottedBubble", "ImaginaryValuesOfDottedBubble"},
                   {g.selfenergy.frequencies.w, Bubble.fermionic_grid.w,
                    g_R.real(), g_R.imag(), g_K.real(), g_K.imag(),
                    s_R.real(), s_R.imag(), s_K.real(), s_K.imag(),
                    Bubble.FermionicBubble.real(), Bubble.FermionicBubble.imag(),
                    DotBubble.FermionicBubble.real(), DotBubble.FermionicBubble.imag()});
#else
    vec<comp> prop (nFER * n_in);
    vec<comp> single_scale (nFER * n_in);
    for (int iv = 0; iv < nFER; ++iv) {
        for (int i_in = 0; i_in < n_in; ++i_in) {
            double v = g.selfenergy.frequencies.w[iv];
            prop[iv * n_in + i_in]         = g.valsmooth(0, v, i_in);
            single_scale[iv * n_in + i_in] = s.valsmooth(0, v, i_in);
        }
    }

    string filename = "/scratch-local/Nepomuk.Ritz/testing_data/FFT_parallelized_full_bubble_in_mom_space_Nq_"
                      + to_string(glb_N_q) + ".h5";
    write_h5_rvecs(filename,
                   {"propagator_frequencies", "bubble_frequencies",
                    "RealValuesOfPropagator", "ImaginaryValuesOfPropagator",
                    "RealValuesOfSingleScale", "ImaginaryValuesOfSingleScale",
                    "RealValuesOfBubble", "ImaginaryValuesOfBubble",
                    "RealValuesOfDottedBubble", "ImaginaryValuesOfDottedBubble"},
                   {g.selfenergy.frequencies.w, Bubble.fermionic_grid.w,
                    prop.real(), prop.imag(), single_scale.real(), single_scale.imag(),
                    Bubble.FermionicBubble.real(), Bubble.FermionicBubble.imag(),
                    DotBubble.FermionicBubble.real(), DotBubble.FermionicBubble.imag()});
#endif
    };

#endif //KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
