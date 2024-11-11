#include "read_NRG_data.hpp"

const std::vector<std::vector<int>> Keldysh_indices = {
        {0, 0, 0, 0},   // iK = 0
        {0, 0, 0, 1},   // iK = 1
        {0, 0, 1, 0},   // iK = 2
        {0, 0, 1, 1},   // iK = 3
        {0, 1, 0, 0},   // iK = 4
        {0, 1, 0, 1},   // iK = 5
        {0, 1, 1, 0},   // iK = 6
        {0, 1, 1, 1},   // iK = 7
        {1, 0, 0, 0},   // iK = 8
        {1, 0, 0, 1},   // iK = 9
        {1, 0, 1, 0},   // iK = 10
        {1, 0, 1, 1},   // iK = 11
        {1, 1, 0, 0},   // iK = 12
        {1, 1, 0, 1},   // iK = 13
        {1, 1, 1, 0},   // iK = 14
        {1, 1, 1, 1}    // iK = 15
};

multidimensional::multiarray<double,7> read_raw_NRG_vertex_component(const std::string& FILENAME,
                                                                     const std::string& DATASET_NAME){
    H5::H5File NRG_file(FILENAME, H5F_ACC_RDONLY);
    H5::DataSet NRG_dataset = NRG_file.openDataSet(DATASET_NAME);
    H5::DataSpace NRG_dataspace = NRG_dataset.getSpace();

    int rank = NRG_dataspace.getSimpleExtentNdims();
    assert(rank==7);

    hsize_t dims[7];
    NRG_dataspace.getSimpleExtentDims(dims, nullptr);

    std::array<std::size_t, 7> length = {dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6]};

    multidimensional::multiarray<double,7> data(length);
    NRG_dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);

    //normalize data by U:
    H5::DataSet NRG_U_set = NRG_file.openDataSet("meta_physical/U");
    double NRG_U;
    NRG_U_set.read(&NRG_U, H5::PredType::NATIVE_DOUBLE);
    for (double & NRG_vertex : data) {
        NRG_vertex = NRG_vertex / NRG_U;
    }
    return data;
}

multidimensional::multiarray<double,4> normalize_NRG_vertex_component(const multidimensional::multiarray<double,7>& vertex_component){

    const size_t Nw  = vertex_component.length()[4];
    const size_t Nv  = vertex_component.length()[5];
    const size_t Nvp = vertex_component.length()[6];

    //initialize array to be returned
    std::array<std::size_t, 4> length = {16, Nw, Nv, Nvp};
    multidimensional::multiarray<double, 4> NRG_component(length);

#pragma omp parallel for schedule(static)
    for (int iK = 0; iK < 16; ++iK) {
        const std::vector<int>& K = Keldysh_indices[iK];
        for (int iw = 0; iw < Nw; ++iw) {
            for (int iv = 0; iv < Nv; ++iv) {
                for (int ivp = 0; ivp < Nvp; ++ivp) {
                    NRG_component.at(iK, iw, iv, ivp) = - vertex_component.at(K[3], K[1], K[2], K[0], iw, iv, ivp);
                    // global minus sign, switch middle Keldysh components.
                }
            }
        }
    }
    return NRG_component;
}

multidimensional::multiarray<double,4> read_NRG_vertex_component(const std::string& FILENAME,
                                                                 const std::string& DATASET_NAME){
    return normalize_NRG_vertex_component(read_raw_NRG_vertex_component(FILENAME, DATASET_NAME));
}

multidimensional::multiarray<double,3> read_raw_NRG_selfenergy(const std::string& FILENAME,
                                                               const std::string& DATASET_NAME){
    H5::H5File NRG_file(FILENAME, H5F_ACC_RDONLY);
    H5::DataSet NRG_dataset = NRG_file.openDataSet(DATASET_NAME);
    H5::DataSpace NRG_dataspace = NRG_dataset.getSpace();

    int rank = NRG_dataspace.getSimpleExtentNdims();
    assert(rank==3);

    hsize_t dims[3];
    NRG_dataspace.getSimpleExtentDims(dims, nullptr);

    std::array<std::size_t, 3> length = {dims[0], dims[1], dims[2]};

    multidimensional::multiarray<double,3> data(length);
    NRG_dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);

    //normalize data by U:
    H5::DataSet NRG_U_set = NRG_file.openDataSet("meta_physical/U");

    double NRG_U;
    NRG_U_set.read(&NRG_U, H5::PredType::NATIVE_DOUBLE);

    for (double & NRG_selfenergy : data) {
        NRG_selfenergy = NRG_selfenergy / NRG_U;
    }

    return data;
}

multidimensional::multiarray<double,2> normalize_NRG_selfenergy(const multidimensional::multiarray<double,3>& raw_selfenergy,
                                                                const double Hartree_shift_in_units_of_U){
    const size_t Nv  = raw_selfenergy.length()[2];

    //initialize array to be returned
    std::array<std::size_t, 2> length = {2, Nv};    // Keldysh component, frequencies
    multidimensional::multiarray<double, 2> selfenergy(length);

    for (int iv = 0; iv < Nv; ++iv) {
        selfenergy.at(0, iv) = raw_selfenergy.at(1, 0, iv) - Hartree_shift_in_units_of_U;   // retarded component
        selfenergy.at(1, iv) = raw_selfenergy.at(0, 0, iv);                                 // Keldysh component
    }

    return selfenergy;
}

std::vector<double> read_NRG_frequency(const std::string& FILENAME,
                                       const std::string& DATASET_NAME,
                                       const bool data_from_MuNRG){
    H5::H5File NRG_file(FILENAME, H5F_ACC_RDONLY);
    H5::DataSet NRG_dataset = NRG_file.openDataSet(DATASET_NAME);
    H5::DataSpace NRG_dataspace = NRG_dataset.getSpace();

    int rank = NRG_dataspace.getSimpleExtentNdims();
    assert(rank==1);

    hsize_t dims[1];
    NRG_dataspace.getSimpleExtentDims(dims, nullptr);
    std::vector<double> NRG_frequencies(dims[0]);

    NRG_dataset.read(NRG_frequencies.data(), H5::PredType::NATIVE_DOUBLE);

    /// normalize frequencies by U:
    const std::string U_name = (data_from_MuNRG ? "meta_physical/U" : "Hamilton_parameters/U");
    H5::DataSet NRG_U_set = NRG_file.openDataSet(U_name);
    double NRG_U;
    NRG_U_set.read(&NRG_U, H5::PredType::NATIVE_DOUBLE);
    for (double & NRG_frequency : NRG_frequencies) {
        NRG_frequency = NRG_frequency / NRG_U;
    }

    return NRG_frequencies;
}

void check_NRG_input(const std::string& NRG_FILENAME, const double U_over_Delta, const double T_in){
    H5::H5File NRG_file(NRG_FILENAME, H5F_ACC_RDONLY);
    H5::DataSet NRG_Delta_set   = NRG_file.openDataSet("meta_physical/Delta");
    H5::DataSet NRG_T_set       = NRG_file.openDataSet("meta_physical/T");
    H5::DataSet NRG_U_set       = NRG_file.openDataSet("meta_physical/U");

    double Delta;
    double T;
    double U;
    NRG_Delta_set.read(&Delta, H5::PredType::NATIVE_DOUBLE);
    NRG_T_set.read(&T, H5::PredType::NATIVE_DOUBLE);
    NRG_U_set.read(&U, H5::PredType::NATIVE_DOUBLE);

    assert(std::abs(U/Delta - U_over_Delta) < 1e-9);
    assert(std::abs(T/U - T_in) < 1e-9);

}



