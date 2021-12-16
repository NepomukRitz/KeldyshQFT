#ifndef KELDYSH_MFRG_WRITE_DATA2FILE_HPP
#define KELDYSH_MFRG_WRITE_DATA2FILE_HPP

#include "../data_structures.hpp" // real/complex vector classes
#include <fstream>           // standard file input/output
#include "H5Cpp.h"           // HDF5 package

// write vectors of real numbers into text (.dat) file
void write_dat_rvecs(std::string path, std::initializer_list<rvec> rvec_list);
void test_write_dat_rvecs();

/**
 * Write vectors of real numbers into hdf5 (.h5) file
 * @param path          : Path to file in which the data is going to be stored
 * @param key_list      : List of the keys, given in the format {"...", "..",...}
 * @param rvec_list     : List of the vectors to be printed. Must be in the order of key_list
 * Notice the vectors to be stored must be real!
 */
void write_h5_rvecs(std::string path, std::initializer_list<std::string> key_list, std::initializer_list<rvec> rvec_list);
void test_write_h5_rvecs();

#endif // KELDYSH_MFRG_WRITE_DATA2FILE_HPP