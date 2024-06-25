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

    return data;
}

multidimensional::multiarray<double,4> normalize_NRG_vertex_component(const multidimensional::multiarray<double,7>& vertex_component,
                                                                      const double U){

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
                    NRG_component.at(iK, iw, iv, ivp) = - vertex_component.at(K[3], K[1], K[2], K[0], iw, iv, ivp) / U;
                    // global minus sign, switch middle Keldysh components, normalize by U. TODO: Check that this is now correct.
                }
            }
        }
    }
    return NRG_component;
}

std::vector<double> read_raw_NRG_frequency(const std::string& FILENAME,
                                           const std::string& DATASET_NAME){
    H5::H5File NRG_file(FILENAME, H5F_ACC_RDONLY);
    H5::DataSet NRG_dataset = NRG_file.openDataSet(DATASET_NAME);
    H5::DataSpace NRG_dataspace = NRG_dataset.getSpace();

    int rank = NRG_dataspace.getSimpleExtentNdims();
    assert(rank==1);

    hsize_t dims[1];
    NRG_dataspace.getSimpleExtentDims(dims, nullptr);
    std::vector<double> NRG_frequencies(dims[0]);

    NRG_dataset.read(NRG_frequencies.data(), H5::PredType::NATIVE_DOUBLE);

    //normalize frequencies by Δ:
    H5::DataSet NRG_Delta_set = NRG_file.openDataSet("meta_physical/Delta");
    H5::DataSpace NRG_Delta_space = NRG_Delta_set.getSpace();

    double NRG_Delta;
    NRG_Delta_set.read(&NRG_Delta, H5::PredType::NATIVE_DOUBLE);

    for (double & NRG_frequency : NRG_frequencies) {
        NRG_frequency = NRG_frequency / NRG_Delta;
    }
    return NRG_frequencies;
}

std::vector<double>
normalize_NRG_frequencies(const std::vector<double> &frequencies, const double U_over_Delta) {
    // in our code, we need the frequencies normalized w.r.t. U, as this is our energy unit.
    std::vector<double> normalized_frequencies(frequencies.size());
    for (int iv = 0; iv < normalized_frequencies.size(); ++iv) {
        normalized_frequencies.at(iv) = frequencies.at(iv) / U_over_Delta;
    }
    return normalized_frequencies;
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



