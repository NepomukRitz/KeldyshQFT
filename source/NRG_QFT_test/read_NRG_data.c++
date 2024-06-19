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
    multidimensional::multiarray<double, 4> data(length);

#pragma omp parallel for schedule(static)
    for (int iK = 0; iK < 16; ++iK) {
        const std::vector<int>& K = Keldysh_indices[iK];
        for (int iw = 0; iw < Nw; ++iw) {
            for (int iv = 0; iv < Nv; ++iv) {
                for (int ivp = 0; ivp < Nvp; ++ivp) {
                    data.at(iK, iw, iv, ivp) = - vertex_component.at(K[3], K[1], K[2], K[0], iw, iv, ivp) / U;
                    // global minus sign, switch middle Keldysh components, normalize by U. TODO: Check that this is now correct.
                }
            }
        }
    }
    return data;
}

std::vector<double> read_raw_NRG_frequency(const std::string& FILENAME,
                                           const std::string& DATASET_NAME){
    // TODO: Implement!
}

std::vector<double> normalize_NRG_frequency(const std::vector<double>& frequency, const double Delta){
    // TODO: Implement!
}


