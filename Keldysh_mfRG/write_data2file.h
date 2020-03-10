#ifndef WRITE_DATA2FILE_H
#define WRITE_DATA2FILE_H

#include "H5Cpp.h"

// write vectors of real numbers into text (.dat) file
void write_dat_rvecs(string path, initializer_list<rvec> rvec_list) {
    ofstream file; file.open(path); // create file at path
    for (auto myvec : rvec_list) { // iterate through vectors
    	for (int i = 0; i<myvec.size(); ++i) file << myvec[i] << " "; // write data into file
    	file << endl; // newline
    }
    file.close(); // close file
}
void test_write_dat_rvecs() {
    rvec x{0., 1., 2.};
    rvec y{-0.5, 0.5, 1.5, 4};
    write_dat_rvecs("test_write_dat_rvecs.dat", {x, y});
}

// write vectors of real numbers into hdf5 (.h5) file
void write_h5_rvecs(string path, initializer_list<string> key_list, initializer_list<rvec> rvec_list) {
	H5::H5File myfile(path, H5F_ACC_TRUNC); // create file at path
	vector<string> keys; // store keys in a vector
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

#endif // WRITE_DATA2FILE_H