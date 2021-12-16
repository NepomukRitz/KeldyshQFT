//
// Created by nepomuk on 27.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
#define KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H

#include<cmath>
#include<complex>
#include<fftw3.h>
#include <cassert>

#include "../data_structures.hpp"
#include "../parameters/master_parameters.hpp"

#include "../utilities/write_data2file.hpp"            // write vectors into hdf5 file (for testing purposes)


int momentum_index(int n_x, int n_y);
void get_n_x_and_n_y(int n, int& n_x, int& n_y);
void get_k_x_and_k_y(int n, double& k_x, double& k_y);
void reduce_index_to_one_eighth_of_BZ(int x, int y, int& n_x, int& n_y);


class Minimal_2D_FFT_Machine {
    fftw_complex *input_propagator, *output_propagator_in_real_space, *input_bubble, *output_bubble_in_momentum_space,
                 *first_propagator_in_real_space, *second_propagator_in_real_space;
    fftw_plan propagator_to_real_space_plan;
    fftw_plan bubble_to_momentum_space_plan;

    vec<comp> return_bubble_values = vec<comp> (glb_N_transfer);

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

    vec<comp> compute_swave_bubble(const vec<comp>& first_g_values, const vec<comp>& second_g_values);
};


#endif //KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
