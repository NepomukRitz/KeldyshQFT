//
// Created by nepomuk on 27.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
#define KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H

#include<tuple>
#include<cmath>
#include<complex>
#include<fftw3.h>
#include "data_structures.h"

int totalNumberOfTransferMomentumPointsToStore();
int momentum_index(int n_x, int n_y);
std::tuple<int, int> get_n_x_and_n_y(int n);

int N_q = 129; // Number of transfer momentum points in one dimension. TODO: Define this globally in parameters.
int N = totalNumberOfTransferMomentumPointsToStore();

int momentum_index(int n_x, int n_y){
    // TODO>: assert n_x < N_q; n_y <= n_x
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

int totalNumberOfTransferMomentumPointsToStore(){
    auto N_qd = (double) N_q;
    double N_q_full = N_qd * (N_qd + 1) / 2;
    return (int) N_q_full;
}


class Minimal_2D_FFT_Machine {
    int points_per_dimension = 2 * (N_q - 1); // N_q set globally
    int total_number_of_points = points_per_dimension * points_per_dimension;

    fftw_complex *input_propagator, *output_propagator_in_real_space, *input_bubble, *output_bubble_in_momentum_space;
    fftw_plan propagator_to_real_space_plan;
    fftw_plan bubble_to_momentum_space_plan;

    void allocate_memory();

    void create_plans();

    void initialize_input_data(); // Currently trivial
    void transform_propagator_to_real_space();

    void calculate_swave_bubble_in_real_space();

    void transform_bubble_to_momentum_space();

    void cleanup();

public:
    Minimal_2D_FFT_Machine() {
        allocate_memory();
        create_plans();
    }
    ~Minimal_2D_FFT_Machine(){
        cleanup();
    }

    void compute_swave_bubble();
};

void Minimal_2D_FFT_Machine::allocate_memory() {
    input_propagator = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total_number_of_points);
    output_propagator_in_real_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total_number_of_points);
    input_bubble = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total_number_of_points);
    output_bubble_in_momentum_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total_number_of_points);
}

void Minimal_2D_FFT_Machine::create_plans() {
    propagator_to_real_space_plan = fftw_plan_dft_2d(points_per_dimension, points_per_dimension,
                                                     input_propagator, output_propagator_in_real_space,
                                                     FFTW_BACKWARD, FFTW_MEASURE);

    bubble_to_momentum_space_plan = fftw_plan_dft_2d(points_per_dimension, points_per_dimension,
                                                     input_bubble, output_bubble_in_momentum_space,
                                                     FFTW_FORWARD, FFTW_MEASURE);
}

void Minimal_2D_FFT_Machine::initialize_input_data() {
    for (int x = 0; x < points_per_dimension; ++x) {
        for (int y = 0; y < points_per_dimension; ++y) {
            // Currently trivial. This is where the actual values of the propagator have to be accessed! Include normalization factor here!
            input_propagator[x * points_per_dimension + y][0] = 0; // Real part
            input_propagator[x * points_per_dimension + y][1] = 0; // Imaginary part
        }
    }
}

void Minimal_2D_FFT_Machine::transform_propagator_to_real_space() {
    fftw_execute(propagator_to_real_space_plan);
    // Actual computation of FFT. Fourier-transformed propagator now in output_propagator_in_real_space.
}

void Minimal_2D_FFT_Machine::calculate_swave_bubble_in_real_space() {
    // Perform element-wise multiplication of the propagators in real space to obtain the s-wave bubble.
    for (int i = 0; i < total_number_of_points; ++i) {
        input_bubble[i][0] = output_propagator_in_real_space[i][0] * output_propagator_in_real_space[i][0]; // Real part
        input_bubble[i][1] = output_propagator_in_real_space[i][1] * output_propagator_in_real_space[i][1]; // Imaginary part
    }
}

void Minimal_2D_FFT_Machine::transform_bubble_to_momentum_space() {
    fftw_execute(bubble_to_momentum_space_plan); // Bubble in momentum space now in output_bubble_in_momentum_space.
}

void Minimal_2D_FFT_Machine::compute_swave_bubble() {
    initialize_input_data();
    transform_propagator_to_real_space();
    calculate_swave_bubble_in_real_space();
    transform_bubble_to_momentum_space();
    std::cout << "S-wave bubble calculated!" << "\n";
}

void Minimal_2D_FFT_Machine::cleanup() {
    fftw_destroy_plan(propagator_to_real_space_plan); fftw_destroy_plan(bubble_to_momentum_space_plan);
    fftw_free(input_propagator); fftw_free(output_propagator_in_real_space);
    fftw_free(input_bubble); fftw_free(output_bubble_in_momentum_space);
}



#endif //KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
