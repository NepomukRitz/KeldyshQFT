//
// Created by nepomuk on 27.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_TEST_MOMENTUM_GRID_H
#define KELDYSH_MFRG_TESTING_TEST_MOMENTUM_GRID_H

#include "../../grids/momentum_grid.h"
#include "../../data_structures.h"
#include "../../r_vertex.h"
#include <random>

TEST_CASE("momentum grid 1", "Test that the conversions between a flat momentum index "
                             "and the two parametrizing one eighth of the BZ actually work."){
    for (int n_x = 0; n_x < glb_N_q; ++n_x) {
        for (int n_y = 0; n_y < n_x+1; ++n_y) {
            int n = momentum_index(n_x, n_y);
            int n_x_recalc, n_y_recalc;
            get_n_x_and_n_y(n, n_x_recalc, n_y_recalc);
            CHECK (n_x - n_x_recalc == 0);
            CHECK (n_y - n_y_recalc == 0);
        }
    }
}

TEST_CASE("momentum grid 2", "Test that the indices n_x and n_y set by get_n_x_and_n_y "
                             "lie on the correct wedge of the BZ."){
    for (int n = 0; n < glb_N_transfer; ++n) {
        int n_x, n_y;
        get_n_x_and_n_y(n, n_x, n_y);
        CHECK (n_x >= 0);
        CHECK (n_y >= 0);
        CHECK (n_x < glb_N_q);
        CHECK (n_y <= n_x);
    }
}

TEST_CASE("momentum grid 3", "Test that the indices n_x and n_y set by reduce_index_to_one_eighth_of_BZ "
                             "lie on the correct wedge of the BZ."){
    for (int x = 0; x < 2 * (glb_N_q - 1); ++x) {
        for (int y = 0; y < 2 * (glb_N_q - 1); ++y) {
            int n_x, n_y;
            reduce_index_to_one_eighth_of_BZ(x, y, n_x, n_y);
            CHECK (n_x >= 0);
            CHECK (n_y >= 0);
            CHECK (n_x < glb_N_q);
            CHECK (n_y <= n_x);
        }
    }
}

TEST_CASE("FFT Machine", "Test that the backward FFT followed by the forward FFT keeps the data invariant."){
    // To generate random numbers
    constexpr int MIN = -10e7;
    constexpr int MAX = 10e7;
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(MIN, MAX);

    int N_FFT = 4 * (glb_N_q - 1) * (glb_N_q - 1);
    vec<comp> v_in1 (glb_N_transfer);
    vec<comp> v_in2 (glb_N_transfer);
    Minimal_2D_FFT_Machine FFT_calculator;

    for (int n = 0; n < glb_N_transfer; ++n) {v_in1[n] = distr(eng) + glb_i * distr(eng);}
    v_in2[0] = N_FFT; // Effective delta-function to become constantly one in real space.
    vec<comp> v_out = FFT_calculator.compute_swave_bubble(v_in1, v_in2);
    vec<comp> v_diff = v_out - v_in1;
    //std::cout << v_diff.max_norm() << "\n";
    CHECK (v_diff.max_norm() < 10e-8); // We require precision to 10e-15 measured w.r.t. the largest value.

    for (int n = 0; n < glb_N_transfer; ++n) {
        v_in1[n] = 0;
        v_in2[n] = distr(eng) + glb_i * distr(eng);
    }
    v_in1[0] = N_FFT;
    v_out = FFT_calculator.compute_swave_bubble(v_in1, v_in2);
    v_diff = v_out - v_in2;
    //std::cout << v_diff.max_norm() << "\n";
    CHECK (v_diff.max_norm() < 10e-8); // We require precision to 10e-15 measured w.r.t. the largest value.
}

TEST_CASE("K1 BZ Average", "Test that the q-sum over the reduced BZ is normalized properly for K1."){
    rvert<comp> test_vertex_comp('a', Lambda_ini);
    rvert<double> test_vertex_double('a', Lambda_ini);

    for (int i_in = 0; i_in < n_in; ++i_in) {
        test_vertex_comp.K1_setvert(0, 0, i_in, 1.);
        test_vertex_double.K1_setvert(0, 0, i_in, 1.);
    }

    comp K1_BZ_sum_comp = test_vertex_comp.K1_BZ_average(0, 0);
    double K1_BZ_sum_double = test_vertex_double.K1_BZ_average(0, 0);

    CHECK(K1_BZ_sum_comp.real() - 1. < 10e-30);
    CHECK(K1_BZ_sum_comp.imag() < 10e-30);
    CHECK(K1_BZ_sum_double - 1. < 10e-30);
}

#endif //KELDYSH_MFRG_TESTING_TEST_MOMENTUM_GRID_H
