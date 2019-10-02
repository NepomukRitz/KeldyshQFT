//
// Created by Elias Walter on 2019-10-01.
//

#ifndef TRY_HDF5_SAVE_HDF5_H
#define TRY_HDF5_SAVE_HDF5_H

#include <iostream>
#include <string>
#include "H5Cpp.h"

using namespace std;
using namespace H5;

/**
 * Create a new HDF5 file
 * @param FILE_NAME : file name
 */
void create_hdf5_file(const H5std_string FILE_NAME) {
    // Try block to detect exceptions raised by any of the calls inside it
    // Exception handling due to HDF5 C++ examples
    try {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        Exception::dontPrint();

        // create a HDF5 file
        // arguments given to constructor:
        // - file name
        // - file access mode: H5F_ACC_TRUNC (if the file already exists, contents will be overwritten)
        //                     H5F_ACC_EXCL (if the file already exists, creation will fail)
        // - file creation property list
        // - file access property list
        H5File file (FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    }
        // catch failure caused by the H5File operations
    catch(FileIException error) {
        error.printErrorStack();
    }

        // catch failure caused by the DataSet operations
    catch(DataSetIException error) {
        error.printErrorStack();
    }

        // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error) {
        error.printErrorStack();
    }
}



void write_hdf5(const H5std_string FILE_NAME) {
    // Try block to detect exceptions raised by any of the calls inside it
    // Exception handling due to HDF5 C++ examples
    try {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        Exception::dontPrint();

        // open an existing HDF5 file
        // arguments given to constructor:
        // - file name
        // - file access mode: H5F_ACC_RDWR (read-write)
        //                     H5F_ACC_RDONLY (read-only)
        // - file creation property list
        // - file access property list
        H5File file (FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT);

        Group vertex (file.createGroup("/Vertex"));
        Group selfenergy (file.createGroup("/SelfEnergy"));



    }
        // catch failure caused by the H5File operations
    catch(FileIException error) {
        error.printErrorStack();
    }

        // catch failure caused by the DataSet operations
    catch(DataSetIException error) {
        error.printErrorStack();
    }

        // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error) {
        error.printErrorStack();
    }
}

#endif //TRY_HDF5_SAVE_HDF5_H
