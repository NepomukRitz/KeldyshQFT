//
// Created by Marcel on 16.11.2021.
//

#ifndef MAIN_CPP_FPP_GRIDS_HPP
#define MAIN_CPP_FPP_GRIDS_HPP

#include "../data_structures.hpp"

class FPP_Grid {
private:
    std::vector<int> Ns_grid;
    rvec specific_points;
    //double lower;
    //double upper;
    bool negative_values;
    bool zero_included;

    // 0: linear, 1: logarithmic
    std::vector<int> grid_types;

public:
    rvec grid_points;

    double calculate_grid_point (double upper, double lower, int i, int n_points, int grid_type){
        double result;
        if (grid_type == 0) { // linear grid
            result = upper + (i - n_points)*(upper - lower)/(n_points-1);
        }
        else if (grid_type == 1) { // logarithmic grid
            result = upper * exp((i - n_points)*log(upper/lower)/(n_points-1));
        }
        else {
            std::cout << "grid type not defined\n";
        }
        return result;
    }

    double calculate_grid_idx (double upper, double lower, double grid_point, int n_points, int grid_type) const {
        double result;
        if (grid_type == 0) { // linear grid
            result = (grid_point*(1-n_points)+lower*n_points-upper)/(lower-upper);
        }
        else if (grid_type == 1) { // logarithmic grid
            result = (n_points-1)*log(grid_point/upper)/log(upper/lower)+n_points;
        }
        else {
            std::cout << "grid type not defined\n";
        }
        return result;
    }
    /*
    FPP_Grid(){
        Ns_grid = {0};
        specific_points = {0,1};
        negative_values = 0;
        zero_included = 0;
        grid_types = {0};
        grid_points = {0};
    };*/

    FPP_Grid() = delete;
    FPP_Grid(const std::vector<int> Ns_grid_in, rvec specific_points_in, bool negative_values_in, bool zero_included_in, std::vector<int> grid_types_in)
    :Ns_grid(Ns_grid_in), specific_points(specific_points_in), negative_values(negative_values_in), zero_included(zero_included_in), grid_types(grid_types_in){
        assert(Ns_grid.size() == grid_types.size());
        assert(Ns_grid.size()+1 == specific_points.size());
        double grid_point;
        double upper, lower;
        int N_grid;
        int grid_type;
        if (negative_values) {
            int j_back;
            for (int j = 0; j < Ns_grid.size(); ++j) {
                j_back = Ns_grid.size() - j - 1;
                upper = specific_points[j_back + 1];
                lower = specific_points[j_back];
                N_grid = Ns_grid[j_back] + 1;
                grid_type = grid_types[j_back];
                /*
                if (grid_type == 0) {
                    grid_measure = upper-lower;
                }
                else if (grid_type == 1){
                    grid_measure = log(upper / lower);
                }
                else {
                    std::cout << "wrong grid type\n";
                }*/
                int i_back;
                for (int i = 1; i < N_grid; ++i) {
                    i_back = N_grid + 1 - i;
                    grid_point = -calculate_grid_point(upper,lower,i_back,N_grid,grid_type);
                    /* if (grid_type == 0){
                        grid_point = -(upper + (i_back - N_grid) * grid_measure / (N_grid - 1));
                    }
                    else if (grid_type == 1){
                        grid_point = -upper * exp((i_back - N_grid) * grid_measure / (N_grid - 1));
                    }
                    else {
                        std::cout << "wrong grid type\n";
                    } */
                    grid_points.push_back(grid_point);
                }
            }
            grid_points.push_back(-specific_points[0]);
        }
        if (zero_included) {
            grid_points.push_back(0.);
        }
        for (int j = 0; j < Ns_grid.size(); ++j) { // positive values
            upper = specific_points[j + 1];
            lower = specific_points[j];
            N_grid = Ns_grid[j] + 1;
            grid_type = grid_types[j];
            /*
            if (grid_type == 1) {
                grid_measure = log(upper / lower);
            }
            else if (grid_type == 0){
                grid_measure = upper-lower;
            }*/
            for (int i = 1; i < N_grid; ++i) {
                grid_point = calculate_grid_point(upper,lower,i,N_grid,grid_type);
                /* if (grid_type == 0){
                    grid_point = upper + (i - N_grid) * grid_measure / (N_grid - 1);
                }
                else if (grid_type == 1){
                    grid_point = upper * exp((i - N_grid) * grid_measure / (N_grid - 1));
                }
                else {
                    std::cout << "wrong grid type\n";
                } */
                grid_points.push_back(grid_point);
            }
        }
        grid_points.push_back(specific_points[specific_points.size() - 1]);
    };

    auto size() const -> int {
        return grid_points.size();
    }

    auto back() const -> double {
        int i = grid_points.size()-1;
        double result;
        result = grid_points[i];
        return result;
    }

    auto operator[] (int i) const -> double {
        return grid_points[i];
    }

    auto grid_transf_inv(double x) const -> int {
        double grid_idx;
        double upper, lower;
        int N_grid;
        int grid_type;

        int sum_Ns_neg = 0;
        int sum_Ns_pos = 0;

        // specific_points.back()
        if (x >= specific_points.back()) { // may be out of range
            return grid_points.size()-1;
        }
        if (negative_values) {
            if (x < -specific_points[0]) {
                sum_Ns_neg = -1;
                int j_back;
                for (int j = 0; j < Ns_grid.size(); ++j) {
                    j_back = Ns_grid.size() - j - 1;
                    N_grid = Ns_grid[j_back] + 1;
                    grid_type = grid_types[j_back];
                    upper = specific_points[j_back + 1];
                    lower = specific_points[j_back];
                    if (x >= -upper and x < -lower) {
                        grid_idx = sum_Ns_neg + N_grid + 1 - calculate_grid_idx(upper, lower, -x, N_grid, grid_type);
                        if (grid_idx < 0){
                            std::cout << "error!\n";
                        }
                        return (int)round(grid_idx);
                    }
                    sum_Ns_neg += N_grid - 1;
                }
            }
            sum_Ns_neg = (int)round(grid_points.size()/2);
            if (x >= -specific_points[0] and x < 0) {
                if (sum_Ns_neg-1 < 0){
                    std::cout << "error!\n";
                }
                return sum_Ns_neg-1;
            }
        }
        if (!zero_included and negative_values) {
            sum_Ns_neg -= 1;
        }
        if ((x >= 0) and (x < specific_points[0])) {
            if (sum_Ns_neg < 0){
                std::cout << "error!\n";
            }
            return sum_Ns_neg;
        }
        if (!zero_included and !negative_values) {
            sum_Ns_pos = -1;
        }
        for (int j = 0; j < Ns_grid.size(); ++j) { // positive values
            upper = specific_points[j + 1];
            lower = specific_points[j];
            N_grid = Ns_grid[j] + 1;
            grid_type = grid_types[j];
            if ((x >= lower) and (x < upper)) {
                grid_idx = sum_Ns_neg + sum_Ns_pos + calculate_grid_idx(upper,lower,x,N_grid,grid_type);
                if (grid_idx < 0){
                    std::cout << "error!\n";
                }
                return (int)round(grid_idx);
            }
            sum_Ns_pos += N_grid - 1;
        }
        if (x < -specific_points.back()) {
            std::cout << "out of grid!\n";
            return 0;
        }

        /*
        if (negative_values) {
            int j_back;
            for (int j = 0; j < Ns_grid.size(); ++j){
                j_back = Ns_grid.size() - j - 1;
                sum_Ns_neg += Ns_grid[j_back];
                if (x > specific_points[j_back] and x < specific_points[j_back+1]) {
                    upper = specific_points[j_back + 1];
                    lower = specific_points[j_back];
                    N_grid = Ns_grid[j_back] + 1;
                    grid_idx = -calculate_grid_idx(upper,lower,x,N_grid,grid_type) + sum_Ns_neg + sum_Ns_pos + int(zero_included);
                }
                sum_Ns_pos += N_grid;
            }
            if (x < 0) {
                // FILL NEGATIVE VALUES
            }
        }
        if (x < 0 and !negative_values) {
            std::cout << "x negative, but no negative grid points\n";
            return 0;
        }
        if (x == 0 and !negative_values) {
            return 0;
        }
        if (x == 0 and negative_values) {
            return sum_Ns_neg + int(zero_included);
        }
        if (x < 0 and negative_values) {

        }
        if (x > 0) {
            for (int j = 0; j < Ns_grid.size(); ++j){
                if (x > specific_points[j] and x < specific_points[j+1]) {
                    upper = specific_points[j+1];
                    lower = specific_points[j];
                    N_grid = Ns_grid[j] + 1;
                    grid_idx = calculate_grid_idx(upper,lower,x,N_grid,grid_type) + sum_Ns_neg + sum_Ns_pos + int(zero_included);
                }
                sum_Ns_pos += N_grid - 1;
            }
            return int(grid_idx);
        }
         */
    }

    auto fconv(double x) const -> int {
        //assert((x >= grid_points[0]) and (x <= grid_points.back());
        /*for (int i = 0; i < grid_points.size(); ++i){
            if (x > grid_points[i] and x < grid_points[i+1]){
                return i;
            }
        }*/
        return grid_transf_inv(x);
    }

    auto get_ws(const int i) const -> double {return grid_points[i];}
};

/*template<int Dim, typename Q, typename Grid_Class, typename F>
class Function_on_Grid {
private:
    std::array<Grid_Class,Dim> grids;
    F func;

public:
    Function_on_Grid(std::array<Grid_Class,Dim> grids_in, F func_in):grids(grids_in), func(func_in){};

    auto operator() (int i) const -> double {
        return func(grids[0][i]);
    }

    auto operator() (int i, int j) const -> double {
        return func(grids[0][i],grids[1][j]);
    }

    auto operator() (int i, int j, int k) const -> double {
        return func(grids[0][i],grids[1][j],grids[2][k]);
    }
};*/

class Function_test {
public:
    Function_test() = default;

    auto operator() (double x) const -> double {
        double result = 1./(sqrt(2*M_PI))*2*exp(-x*x/2);
        return result;
    }

    auto operator() (double x, double y) const -> double {
        double result = 1./(sqrt(2*M_PI))*2*exp(-(x*x+y*y)/2);
        return result;
    }

    auto operator() (double x, double y, double z) const -> double {
        double result = 1./(sqrt(2*M_PI))*2*exp(-(x*x+y*y+z*z)/2);
        return result;
    }
};

class Function2D_test {
public:
    Function2D_test() = default;

    auto operator() (double x, double y) const -> double {
        double result = 1./(sqrt(2*M_PI))*2*exp(-x*x/2);
        return result;
    }
};

#endif //MAIN_CPP_FPP_GRIDS_HPP
