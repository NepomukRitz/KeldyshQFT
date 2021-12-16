#include "write_data2file.hpp"

void write_dat_rvecs(std::string path, std::initializer_list<rvec> rvec_list) {
    std::ofstream file; file.open(path); // create file at path
    for (auto myvec : rvec_list) { // iterate through vectors
        for (int i = 0; i<myvec.size(); ++i) file << myvec[i] << " "; // write data into file
        file << std::endl; // newline
    }
    file.close(); // close file
}
void test_write_dat_rvecs() {
    rvec x{0., 1., 2.};
    rvec y{-0.5, 0.5, 1.5, 4};
    write_dat_rvecs("test_write_dat_rvecs.dat", {x, y});
}

void write_h5_rvecs(std::string path, std::initializer_list<std::string> key_list, std::initializer_list<rvec> rvec_list) {
    H5::H5File myfile(path, H5F_ACC_TRUNC); // create file at path
    std::vector<std::string> keys; // store keys in a vector
    for (auto mykey : key_list) keys.push_back(mykey); // fill keys vector
    int i = 0; // iterator for keys
    hsize_t dim_vec[1]; // dimension of vector, to be updated
    for (auto myvec : rvec_list) { // iterate through vectors
        dim_vec[0] = myvec.size(); // updated vector length
        H5::DataSpace mydataspace(1, dim_vec); // create dataspace to store vector
        H5::DataSet mydataset = myfile.createDataSet(keys[i], H5::PredType::NATIVE_DOUBLE, mydataspace); // put vector in dataset
        mydataset.write(&myvec[0], H5::PredType::NATIVE_DOUBLE); // write dataset into file
        mydataset.close(); // close dataset
        mydataspace.close(); // close dataspace
        ++i; // iterator for keys
    }
    myfile.close(); // close file
}
void test_write_h5_rvecs() {
    rvec x{0., 1., 2.};
    rvec y{-0.5, 0.5, 1.5, 4};
    write_h5_rvecs("test_write_h5_rvecs.h5", {"xvector", "yvector"}, {x, y});
}
