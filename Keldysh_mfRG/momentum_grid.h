//
// Created by nepomuk on 27.04.21.
//

#ifndef KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
#define KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H

#include<tuple>
#include<cmath>

int totalNumberOfTransferMomentumPointsToStore();
int momentum_index(int n_x, int n_y);
std::tuple<int, int> get_n_x_and_n_y(int n);

int N_q = 4000; // Number of transfer momentum points in one dimension. TODO: Define this globally in parameters.
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



#endif //KELDYSH_MFRG_TESTING_MOMENTUM_GRID_H
