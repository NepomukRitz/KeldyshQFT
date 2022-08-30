#include "hdf5_routines.hpp"

// Create the memory data type for storing complex numbers in file


H5::H5File create_hdf_file(const std::string & filename) {
    if (mpi_world_rank() == 0) {
        H5::H5File outfile(filename, H5F_ACC_TRUNC);
        return outfile;
    }
    else {
        H5::H5File file;
        return file;
    }
}
H5::H5File open_hdf_file_readOnly(const std::string & filename) {
    if (mpi_world_rank() == 0) {
        H5::H5File outfile(filename, H5F_ACC_RDONLY);
        return outfile;
    }
    else {
        H5::H5File file;
        return file;
    }
}
H5::H5File open_hdf_file_readWrite(const std::string & filename) {
    if (mpi_world_rank() == 0) {
        H5::H5File outfile(filename, H5F_ACC_RDWR);
        return outfile;
    }
    else {
        H5::H5File file;
        return file;
    }
}

void close_hdf_file(H5::H5File & file) {
    if (mpi_world_rank() == 0) {
        file.close();
    }
}


H5::CompType def_mtype_comp() {
    H5::CompType mtype(sizeof(h5_comp));
    mtype.insertMember(RE, HOFFSET(h5_comp, re), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember(IM, HOFFSET(h5_comp, im), H5::PredType::NATIVE_DOUBLE);
    return mtype;
}
H5::DSetCreatPropList def_proplist_comp() {
    h5_comp fillvalue_vert;
    fillvalue_vert.re = 0;
    fillvalue_vert.im = 0;
    H5::DSetCreatPropList plist_vert;
    H5::CompType mtype_comp = def_mtype_comp();
    plist_vert.setFillValue(mtype_comp, &fillvalue_vert);
    return plist_vert;
}

hsize_t h5_cast(int dim) {
    return static_cast<hsize_t>(dim);
}






rvec read_Lambdas_from_hdf(const H5std_string FILE_NAME){



    // Open the file. Access rights: read-only
    H5::H5File file = open_hdf_file_readOnly(FILE_NAME);

    H5::DataSet lambda_dataset = file.openDataSet("lambdas");


    // Create the data spaces for the data sets in file and for buffer objects
    H5::DataSpace dataSpace_Lambda = lambda_dataset.getSpace();
    //
    //Get the dimension size of each dimension in the dataspace and
    //display them.
    //
    hsize_t dims_out[2];
    dataSpace_Lambda.getSimpleExtentDims( dims_out, NULL);
    //std::cout << "rank " << rank << ", dimensions " <<
    //     (unsigned long)(dims_out[0]) << " x " <<
    //     (unsigned long)(dims_out[1]) << std::endl;


    H5::DataSpace dataSpace_Lambda_buffer(1, dims_out);

    /// load data into buffer
    double Lambdas_arr[dims_out[0]];
    lambda_dataset.read(Lambdas_arr, H5::PredType::NATIVE_DOUBLE,
                        dataSpace_Lambda_buffer, dataSpace_Lambda);

    rvec Lambdas (Lambdas_arr, Lambdas_arr + sizeof(Lambdas_arr) / sizeof(double));

    // Terminate
    lambda_dataset.close();
    close_hdf_file(file);

    return Lambdas;


}


bool test_read_write_data_hdf(bool verbose) {
    //if (verbose) utils::print("Testing HDF input and output:", true);

    vec<comp> vector_original(10, glb_i);
    std::array<std::size_t,3> dims_mularr = {1,2,3};
    multidimensional::multiarray<comp,3> mularr_original(dims_mularr, glb_i);

    std::string testfile_name = "test.h5";
    //if (verbose) utils::print("File name: " + testfile_name, true);

    // write to hdf file
    H5::H5File file_out = create_hdf_file(testfile_name.c_str());
    //if (verbose) utils::print("Created file.", true);
    //H5::H5std_string FILE_NAME(data_dir + "test.h5");
    H5::Group group_out(file_out.createGroup("/data"));
    //if (verbose) utils::print("Created group.", true);

    write_to_hdf(group_out, "vector", vector_original, false);
    write_to_hdf(group_out, "mularr", mularr_original, false);
    //if (verbose) utils::print("Written to hdf.", true);

    // read from hdf file
    hid_t file_id = H5Fopen(testfile_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5::Group group = H5Gopen(file_id, "data", H5P_DEFAULT);
    //if (verbose) utils::print("Loaded group.", true);
    vec<comp> vector_input;
    read_from_hdf<comp>(group, "vector", vector_input);
    //if (verbose) utils::print("Loaded vector.", true);
    multidimensional::multiarray<comp,3> mularr_input;
    read_from_hdf<comp,3>(group, "mularr", mularr_input);
    //if (verbose) utils::print("Loaded multiarray.", true);

    //if (verbose) utils::print("Length of original vec: ", vector_original.size(), true);
    //if (verbose) utils::print("Length of loaded vec: ", vector_input.size(), true);

    comp deviation_vec = (vector_input - vector_original).max_norm();
    if (std::abs(deviation_vec) < 1e-10) utils::print("Read / write of vector to HDF file successful.", true);
    else utils::print("PROBLEM during read / write of vector to HDF file.", true);

    comp deviation_mularr=(mularr_input -mularr_original).maxabs();
    if (std::abs(deviation_mularr) < 1e-10) utils::print("Read / write of multiarray to HDF file successful.", true);
    else utils::print("PROBLEM during read / write of multiarray to HDF file.", true);



    //H5::H5std_string FILE_NAME(data_dir + "test.h5");
    H5::Group group_out_L(file_out.createGroup("/data_with_LambdaLayers"));
    //if (verbose) putils::print("Created group.", true);

    write_to_hdf_LambdaLayer(group_out_L, "vector", vector_original, 1, 10,  false);
    write_to_hdf_LambdaLayer(group_out_L, "mularr", mularr_original, 1, 10,  false);
    //if (verbose) putils::print("Written to hdf LambdaLayer.", true);

    // read from hdf file
    //hid_t file_id = H5Fopen(testfile_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5::Group group_L = H5Gopen(file_id, "data_with_LambdaLayers", H5P_DEFAULT);
    vec<comp> vector_input_L;
    read_from_hdf_LambdaLayer<comp>(group_L, "vector", vector_input_L, 1);
    multidimensional::multiarray<comp,3> mularr_input_L;
    read_from_hdf_LambdaLayer<comp>(group_L, "mularr", mularr_input_L, 1);


    comp deviation_vec_L = (vector_input_L - vector_original).max_norm();
    comp deviation_mularr_L=(mularr_input_L -mularr_original).maxabs();


    bool passed = true;
    if (verbose) {
        if (std::abs(deviation_vec_L) < 1e-10)
            utils::print("Read / write of vector to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            utils::print("PROBLEM during read / write of vector to LambdaLayer of HDF file.", true);
        }
        if (std::abs(deviation_mularr_L) < 1e-10)
            utils::print("Read / write of multiarray to LambdaLayer ofHDF file successful.", true);
        else {
            passed = false;
            utils::print("PROBLEM during read / write to LambdaLayer of multiarray to HDF file.", true);
        }
    }

    return passed;
}


