//
// Created by nepomuk on 08.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
#define KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H

#include <iostream>
#include <string>
#include "../parameters/master_parameters.hpp"
#include "../data_structures.hpp"
#include "../correlation_functions/state.hpp"
#include "../perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "../bubble/bubble_function.hpp"
#include "../correlation_functions/two_point/propagator.hpp"
#include "../utilities/write_data2file.hpp"            // write vectors into hdf5 file
#include "../utilities/util.hpp"
#include "../symmetries/Keldysh_symmetries.hpp"
#include "../bubble/bubble.hpp"
#include "../bubble/precalculated_bubble.hpp"


template<typename Q>
class test_PrecalculateBubble{
public:
    Propagator<Q> g;
    Propagator<Q> s;

    test_PrecalculateBubble(): g(Lambda_ini, 'g'), s(Lambda_ini, 's'){
        State<comp> testing_state (Lambda_ini);
        testing_state.initialize();
        sopt_state(testing_state, Lambda_ini);

        g = Propagator<Q>(Lambda_ini, testing_state.selfenergy, 'g');
        s = Propagator<Q>(Lambda_ini, testing_state.selfenergy, 's');
        // g = Propagator(Lambda_ini, 'g');
        // s = Propagator(Lambda_ini, 's');
    }

    void perform_test();
    double iterate_through_channels();
    double find_largest_deviation_from_bubble(bool dot_in, char channel);
    void save_data(std::string& filename,
                   vec<Q>& ValuesOfPreBubble, vec<Q>& ValuesOfUsualBubble,
                   vec<double>& AbsoluteDeviations);
    std::string build_filename(bool dot_in, char channel);
    void test_for_zero_value(Q& value, int& number_of_zero_values);

};

template<typename Q> void test_PrecalculateBubble<Q>::perform_test(){
    double max_error = iterate_through_channels();
    if (max_error != 0.){
        std::string deviation = std::to_string(max_error);
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
    double larger_error = std::max(std::abs(error_undot), std::abs(error_dot));
    if (larger_error > max_error){
        max_error = larger_error;
    }
    return max_error;
}

template<typename Q> double test_PrecalculateBubble<Q>::find_largest_deviation_from_bubble(const bool dot_in, const char channel){
    PrecalculatedBubble<Q> Pre_Bubble(g, s, dot_in, channel);
    Bubble<Q> Usual_Bubble (g, s, dot_in);

    vec<Q> ValuesOfPreBubble (glb_number_of_Keldysh_components_bubble * nBOS * nFER * n_in);
    vec<Q> ValuesOfUsualBubble (glb_number_of_Keldysh_components_bubble * nBOS * nFER * n_in);

    vec<double> AbsoluteDeviations (glb_number_of_Keldysh_components_bubble * nBOS * nFER * n_in);

    double largest_deviation = 0;
    int number_of_zero_pre_results = 0;

    for (int iK = 0; iK < glb_number_of_Keldysh_components_bubble; ++iK) {
        for (int iw = 0; iw < nBOS; ++iw) {
            const double w = g.selfenergy.frequencies.get_ws(iw);
            for (int ivpp = 0; ivpp < nFER; ++ivpp) {
                const double vpp = g.selfenergy.frequencies.get_ws(ivpp);
                for (int i_in = 0; i_in < n_in; ++i_in) {
                    Q Bubble_Value = Usual_Bubble.value(iK, w, vpp, i_in, channel);
                    Q Pre_Bubble_Value = Pre_Bubble.value(iK, w, vpp, i_in);
                    //test_for_zero_value(Pre_Bubble_Value, number_of_zero_pre_results);

                    ValuesOfPreBubble[Pre_Bubble.composite_index(iK, iw, ivpp, i_in)] = Pre_Bubble_Value;
                    ValuesOfUsualBubble[Pre_Bubble.composite_index(iK, iw, ivpp, i_in)] = Bubble_Value;

                    Q deviation = Pre_Bubble_Value - Bubble_Value;
                    AbsoluteDeviations[Pre_Bubble.composite_index(iK, iw, ivpp, i_in)] = std::abs(deviation);
                    if (std::abs(deviation) > largest_deviation) {
                        largest_deviation = std::abs(deviation);
                        // std::string dev_string = std::to_string(largest_deviation);
                        // std::cout << dev_string << "\n";
                    }
                }
            }
        }
    }
    std::cout << "Number of zero results encountered for precalculated bubble = "
    + std::to_string(number_of_zero_pre_results) + "\n";
    std::string filename = build_filename(dot_in, channel);
    save_data(filename, ValuesOfPreBubble, ValuesOfUsualBubble, AbsoluteDeviations);
    return largest_deviation;
}

template<typename Q>
void
test_PrecalculateBubble<Q>::save_data(std::string& filename,
                                      vec<Q>& ValuesOfPreBubble,
                                      vec<Q>& ValuesOfUsualBubble,
                                      vec<double>& AbsoluteDeviations) {
write_h5_rvecs(filename,
               {"Frequencies",
                "RealValuesOfPrecalculatedBubble", "RealValuesOfUsualBubble",
                "ImaginaryValuesOfPrecalculatedBubble", "ImaginaryValuesOfUsualBubble",
                "AbsoluteDeviations"},
               {g.selfenergy.frequencies.ws,
                ValuesOfPreBubble.real(), ValuesOfUsualBubble.real(),
                ValuesOfPreBubble.imag(), ValuesOfUsualBubble.imag(),
                AbsoluteDeviations});
}

template<typename Q>
std::string test_PrecalculateBubble<Q>::build_filename(bool dot_in, char channel) {
    std::string filename = "/scratch-local/Nepomuk.Ritz/testing_data/true_test_PrecalculateBubble_";
    //std::string filename = "../../Data/test_PrecalculateBubble_";
    filename += channel;
    if (dot_in) {filename += "_dot";}
    filename += "_nK=" + std::to_string(glb_number_of_Keldysh_components_bubble)
                + "_nBOS=" + std::to_string(nBOS)
                + "_nFER=" + std::to_string(nFER)
                + "_n_in=" + std::to_string(n_in)
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
                     + std::to_string(value.real()) + "\n";
        std::cout << "Im Pre-Bubble = "
                     + std::to_string(value.imag()) + "\n";
    }
}

template<typename Q>
class Runtime_comparison{
    Propagator<Q> g;
    Propagator<Q> s;
    PrecalculatedBubble<Q> Pre_Bubble; // free versions of the bubbles
    Bubble<Q> Usual_Bubble;
public:
    Runtime_comparison(): g (Lambda_ini, 'g'), s (Lambda_ini, 's'),
                       Pre_Bubble (g, s, 0, 'a'), Usual_Bubble (g, s, false){
     State<comp> testing_state (Lambda_ini);
     testing_state.initialize();
     sopt_state(testing_state, Lambda_ini);

     g = Propagator<Q> (Lambda_ini, testing_state.selfenergy,'g'); // this is done to obtain the frequency grid
     s = Propagator<Q> (Lambda_ini, testing_state.selfenergy,'s');
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
    std::string filename = "/scratch-local/Nepomuk.Ritz/testing_data/runtime_comparisons.h5";
    write_h5_rvecs(filename,
                   {"runtimes_usual", "runtimes_precalculated"},
                   {times_usual, times_precalculated});
}

template<typename Q>
double Runtime_comparison<Q>::run_iterations(int iterations, bool precalculated) {
    double starting_time = get_time();
    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (int iK = 0; iK < glb_number_of_Keldysh_components_bubble; ++iK) {
            for (int iw = 0; iw < nBOS; ++iw) {
                const double w = g.selfenergy.frequencies.get_ws(iw);
                for (int ivpp = 0; ivpp < nFER; ++ivpp) {
                    const double vpp = g.selfenergy.frequencies.get_ws(ivpp);
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
    double Lambda = 0.01;
    Propagator<comp> g (Lambda, 'g');
    Propagator<comp> s (Lambda, 's');

    double starting_time = get_time();
    PrecalculatedBubble<comp> DotBubble (g, s, true);
    double end_time = get_time();
    double diff = (end_time - starting_time); // time given in seconds
    std::cout << "Time for differentiated Bubble = " << diff << " s." << "\n";

    starting_time =  get_time();
    PrecalculatedBubble<comp> Bubble (g, s, false);
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
            double v = g.selfenergy.frequencies.get_ws(iv);
            g_R[iv * n_in + i_in] = g.valsmooth(0, v, i_in);
            g_K[iv * n_in + i_in] = g.valsmooth(1, v, i_in);
            s_R[iv * n_in + i_in] = s.valsmooth(0, v, i_in);
            s_K[iv * n_in + i_in] = s.valsmooth(1, v, i_in);
        }
    }

    std::string filename = "/scratch-local/Nepomuk.Ritz/testing_data/KELDYSH_bubble_in_mom_space_Nq_"
                      + std::to_string(glb_N_q) + "_Lambda_" + std::to_string(Lambda) + ".h5";
    write_h5_rvecs(filename,
                   {"propagator_frequencies", "bubble_frequencies",
                    "RealValuesOfRetardedPropagator", "ImaginaryValuesOfRetardedPropagator",
                    "RealValuesOfKeldyshPropagator", "ImaginaryValuesOfKeldyshPropagator",
                    "RealValuesOfRetardedSingleScale", "ImaginaryValuesOfRetardedSingleScale",
                    "RealValuesOfKeldyshSingleScale", "ImaginaryValuesOfKeldyshSingleScale",
                    "RealValuesOfBubble", "ImaginaryValuesOfBubble",
                    "RealValuesOfDottedBubble", "ImaginaryValuesOfDottedBubble"},
                   {g.selfenergy.frequencies.get_ws_vec(), Bubble.fermionic_grid.get_ws_vec(),
                    g_R.real(), g_R.imag(), g_K.real(), g_K.imag(),
                    s_R.real(), s_R.imag(), s_K.real(), s_K.imag(),
                    Bubble.FermionicBubble.real(), Bubble.FermionicBubble.imag(),
                    DotBubble.FermionicBubble.real(), DotBubble.FermionicBubble.imag()});
#else
    vec<comp> prop (nFER * n_in);
    vec<comp> single_scale (nFER * n_in);
    for (int iv = 0; iv < nFER; ++iv) {
        for (int i_in = 0; i_in < n_in; ++i_in) {
            double v = g.selfenergy.frequencies.get_ws(iv);
            prop[iv * n_in + i_in]         = g.valsmooth(0, v, i_in);
            single_scale[iv * n_in + i_in] = s.valsmooth(0, v, i_in);
        }
    }

    std::string filename = "/scratch-local/Nepomuk.Ritz/testing_data/FFT_parallelized_full_bubble_in_mom_space_Nq_"
                      + std::to_string(glb_N_q) + ".h5";
    write_h5_rvecs(filename,
                   {"propagator_frequencies", "bubble_frequencies",
                    "RealValuesOfPropagator", "ImaginaryValuesOfPropagator",
                    "RealValuesOfSingleScale", "ImaginaryValuesOfSingleScale",
                    "RealValuesOfBubble", "ImaginaryValuesOfBubble",
                    "RealValuesOfDottedBubble", "ImaginaryValuesOfDottedBubble"},
                   {g.selfenergy.frequencies.get_ws_vec(), Bubble.fermionic_grid.get_ws_vec(),
                    prop.real(), prop.imag(), single_scale.real(), single_scale.imag(),
                    Bubble.FermionicBubble.real(), Bubble.FermionicBubble.imag(),
                    DotBubble.FermionicBubble.real(), DotBubble.FermionicBubble.imag()});
#endif
    }

void save_PreBubble_in_freq_space(const PrecalculatedBubble<comp> Pi, const int i_in){
    const std::string directory = "/project/th-scratch/n/Nepomuk.Ritz/PhD_data/SOPT/bare_bubble/";
    const std::string filename  = directory + "bare_bubble_on_fermionic_grid_in_" + std::to_string(i_in) + "n_in_" + std::to_string(n_in) + ".h5";
    vec<comp> pi (glb_number_of_Keldysh_components_bubble * nFER * nFER);
    for (int iK_bubble = 0; iK_bubble < glb_number_of_Keldysh_components_bubble; ++iK_bubble) {
        for (int iv1 = 0; iv1 < nFER; ++iv1) {
            for (int iv2 = 0; iv2 < nFER; ++iv2) {
                const double v1 = Pi.fermionic_grid.get_ws(iv1);
                const double v2 = Pi.fermionic_grid.get_ws(iv2);
                comp val = Pi.value_on_fermionic_grid(iK_bubble, v1, v2, i_in);
                pi[iK_bubble * nFER * nFER + iv1 * nFER + iv2] = val;
            }
        }
    }
    write_h5_rvecs(filename,
                   {"frequencies", "RealBubble", "ImagBubble"},
                   {Pi.fermionic_grid.get_ws_vec(), pi.real(), pi.imag()});
}

#endif //KELDYSH_MFRG_TESTING_TEST_PRECALCULATEBUBBLE_H
