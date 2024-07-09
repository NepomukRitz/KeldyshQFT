#ifndef KELDYSH_MFRG_READ_NRG_DATA_HPP
#define KELDYSH_MFRG_READ_NRG_DATA_HPP

#include "data_structures.hpp"
#include "multidimensional/multiarray.hpp"
#include "H5Cpp.h"
#include <cassert>

multidimensional::multiarray<double,7> read_raw_NRG_vertex_component(const std::string& FILENAME,
                                                                     const std::string& DATASET_NAME);

multidimensional::multiarray<double,4> normalize_NRG_vertex_component(const multidimensional::multiarray<double,7>& vertex_component,
                                                                      double U);

multidimensional::multiarray<double,3> read_raw_NRG_selfenergy(const std::string& FILENAME,
                                                               const std::string& DATASET_NAME);

multidimensional::multiarray<double,2> normalize_NRG_selfenergy(const multidimensional::multiarray<double,3>& selfenergy_component,
                                                                double Hartree_shift_in_units_of_U = 0.0);

std::vector<double> read_NRG_frequency(const std::string& FILENAME,
                                       const std::string& DATASET_NAME);

void check_NRG_input(const std::string& NRG_FILENAME, double U_over_Delta, double T_in);


#endif //KELDYSH_MFRG_READ_NRG_DATA_HPP
