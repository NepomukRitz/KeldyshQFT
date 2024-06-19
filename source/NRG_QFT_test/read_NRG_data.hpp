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

std::vector<double> read_raw_NRG_frequency(const std::string& FILENAME,
                                           const std::string& DATASET_NAME);

std::vector<double> normalize_NRG_frequency(const std::vector<double>& frequency, double Delta);



#endif //KELDYSH_MFRG_READ_NRG_DATA_HPP
