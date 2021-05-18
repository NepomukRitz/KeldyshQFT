//
// Created by nepomuk on 27.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
#define KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H

#define _USE_MATH_DEFINES
#include<cmath>
#include<complex>
#include<tuple>
#include<fftw3.h>
#include <cassert>

#include "data_structures.h"
#include "propagator.h"
#include "parameters.h"

#include "write_data2file.h"            // write vectors into hdf5 file (for testing purposes)


int momentum_index(int n_x, int n_y);
std::tuple<int, int> get_n_x_and_n_y(int n);
std::tuple<int, int> reduce_index_to_one_eighth_of_BZ(int x, int y);

int momentum_index(int n_x, int n_y){
    // TODO>: Cover assertions in unit test.
    assert (n_x >= 0);
    assert (n_y >= 0);
    assert (n_x < N_q);
    assert (n_y <= n_x);
    auto n_xd = (double) n_x;
    auto n_yd = (double) n_y;
    double n;
    n = n_yd + n_xd * (n_xd + 1) / 2;
    return (int) n;
}

std::tuple<int, int> get_n_x_and_n_y(int n) {
    auto nd = (double) n;
    double n_xd;
    n_xd = (sqrt(1 + 8 * nd) - 1) / 2;

    int n_x;
    n_x = (int) n_xd; // Always rounds down, i.e. this is a downstairs Gauss-bracket, as required.

    int n_y;
    double n_yd;
    n_xd = (double) n_x;
    n_yd = nd - n_xd * (n_xd + 1) / 2;
    n_y = (int) std::round(n_yd); // n_yd should already be an very close to an integer.
                                  // Use round to obtain this integer and cast to int type.

    return std::make_tuple(n_x, n_y);
}

std::tuple<double, double> get_k_x_and_k_y(int n){
    int n_x; int n_y;
    double k_x; double k_y;

    std::tie(n_x, n_y) = get_n_x_and_n_y(n);
    k_x = M_PI * n_x / (N_q - 1);
    k_y = M_PI * n_y / (N_q - 1);

    return std::make_tuple(k_x, k_y);
}

/*
 * Given a pair of indices labelling any point in the square of the BZ,
 * what is the corresponding index in the part that is actually stored?
 */
std::tuple<int, int> reduce_index_to_one_eighth_of_BZ(int x, int y){
    // TODO: Cover assertions in unit test.
    assert (x >= 0);
    assert (y >= 0);
    assert (x < 2 * (N_q - 1));
    assert (y < 2 * (N_q - 1));
    int n_x = x;
    int n_y = y;
    if (n_x > (N_q - 1)){n_x = 2 * (N_q - 1) - n_x;} // mirror along x-direction
    if (n_y > (N_q - 1)){n_y = 2 * (N_q - 1) - n_y;} // mirror along y-direction
    if (n_y > n_x){ //mirror along diagonal (i.e. swap)
        int temp_n_x = n_x;
        n_x = n_y;
        n_y = temp_n_x;
    }
    return std::make_tuple(n_x, n_y);
};

class Minimal_2D_FFT_Machine {
    int points_per_dimension = 2 * (N_q - 1); // N_q set globally
    int N_FFT = points_per_dimension * points_per_dimension;

    fftw_complex *input_propagator, *output_propagator_in_real_space, *input_bubble, *output_bubble_in_momentum_space,
                 *first_propagator_in_real_space, *second_propagator_in_real_space;
    fftw_plan propagator_to_real_space_plan;
    fftw_plan bubble_to_momentum_space_plan;

    vec<comp> return_bubble_values = vec<comp> (N);

    void allocate_memory();
    void create_plans();

    void initialize_input_data(const vec<comp>& g_values);
    void transform_propagator_to_real_space();
    void calculate_swave_bubble_in_real_space();
    void transform_bubble_to_momentum_space();
    void prepare_bubble_values_to_return();

public:
    Minimal_2D_FFT_Machine() {
        allocate_memory();
        create_plans();
    }
    ~Minimal_2D_FFT_Machine(){
        fftw_cleanup();
    }

    vec<comp> compute_swave_bubble(vec<comp>& first_g_values, vec<comp>& second_g_values);
};

vec<comp> Minimal_2D_FFT_Machine::compute_swave_bubble(vec<comp> &first_g_values, vec<comp> &second_g_values) {
    initialize_input_data(first_g_values);

    transform_propagator_to_real_space();
    for (int i = 0; i < N_FFT; ++i) {
        first_propagator_in_real_space[i][0] = output_propagator_in_real_space[i][0];
        first_propagator_in_real_space[i][1] = output_propagator_in_real_space[i][1];
    }

    initialize_input_data(second_g_values);
    transform_propagator_to_real_space();
    for (int i = 0; i < N_FFT; ++i) {
        second_propagator_in_real_space[i][0] = output_propagator_in_real_space[i][0];
        second_propagator_in_real_space[i][1] = output_propagator_in_real_space[i][1];
    }

    calculate_swave_bubble_in_real_space();
    transform_bubble_to_momentum_space();
    prepare_bubble_values_to_return();

    return return_bubble_values;
}

void Minimal_2D_FFT_Machine::initialize_input_data(const vec<comp> &g_values) {
    // TODO: Currently only the s-wave bubble out of two full propagators can be computed. Extend to differentiated bubbles!
    for (int x = 0; x < points_per_dimension; ++x) {
        for (int y = 0; y < points_per_dimension; ++y) {
            int n_x, n_y;
            std::tie(n_x, n_y) = reduce_index_to_one_eighth_of_BZ(x, y);
            // This is where the actual values of the propagator have to be accessed! Include normalization factor here!
            input_propagator[x * points_per_dimension + y][0] = g_values[momentum_index(n_x, n_y)].real() / N_FFT; // Real part
            input_propagator[x * points_per_dimension + y][1] = g_values[momentum_index(n_x, n_y)].imag() / N_FFT; // Imaginary part
        }
    }
}

void Minimal_2D_FFT_Machine::transform_propagator_to_real_space() {
    fftw_execute(propagator_to_real_space_plan);
    // Actual computation of FFT. Fourier-transformed propagator now in output_propagator_in_real_space.
}

void Minimal_2D_FFT_Machine::calculate_swave_bubble_in_real_space() {
    // Perform element-wise multiplication of the propagators in real space to obtain the s-wave bubble.
    for (int i = 0; i < N_FFT; ++i) {
        input_bubble[i][0] = first_propagator_in_real_space[i][0] * second_propagator_in_real_space[i][0]
                           - first_propagator_in_real_space[i][1] * second_propagator_in_real_space[i][1]; // Real part
        input_bubble[i][1] = first_propagator_in_real_space[i][0] * second_propagator_in_real_space[i][1]
                           + first_propagator_in_real_space[i][1] * second_propagator_in_real_space[i][0]; // Imaginary part
    }
}

void Minimal_2D_FFT_Machine::transform_bubble_to_momentum_space() {
    fftw_execute(bubble_to_momentum_space_plan); // Bubble in momentum space now in output_bubble_in_momentum_space.
}

void Minimal_2D_FFT_Machine::prepare_bubble_values_to_return() {
    for (int n_x = 0; n_x < N_q; ++n_x) {
        for (int n_y = 0; n_y < n_x+1; ++n_y) {
            return_bubble_values[momentum_index(n_x, n_y)] =
                    output_bubble_in_momentum_space[n_x * points_per_dimension + n_y][0] + glb_i *
                    output_bubble_in_momentum_space[n_x * points_per_dimension + n_y][1];
        }
    }
}

void Minimal_2D_FFT_Machine::allocate_memory() {
    input_propagator = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);
    output_propagator_in_real_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);
    input_bubble = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);
    output_bubble_in_momentum_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);

    first_propagator_in_real_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);
    second_propagator_in_real_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);
}

void Minimal_2D_FFT_Machine::create_plans() {
    propagator_to_real_space_plan = fftw_plan_dft_2d(points_per_dimension, points_per_dimension,
                                                     input_propagator, output_propagator_in_real_space,
                                                     FFTW_BACKWARD, FFTW_MEASURE);

    bubble_to_momentum_space_plan = fftw_plan_dft_2d(points_per_dimension, points_per_dimension,
                                                     input_bubble, output_bubble_in_momentum_space,
                                                     FFTW_FORWARD, FFTW_MEASURE);
}


#endif //KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
