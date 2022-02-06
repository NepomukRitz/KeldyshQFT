/**
 * Functions to write/read a State object to/from an HDF5 file.
 */

#ifndef KELDYSH_MFRG_HDF5_ROUTINES_HPP
#define KELDYSH_MFRG_HDF5_ROUTINES_HPP

#include <stdexcept>
#include <cmath>
#include <vector>
#include "util.hpp"               // printing text
#include "../parameters/master_parameters.hpp"         // system parameters (necessary for vector lengths etc.)
#include "../data_structures.hpp"    // comp data type, std::real/complex vector class
#include "../grids/frequency_grid.hpp"     // store frequency grid parameters
#include "H5Cpp.h"              // HDF5 functions
#include "../correlation_functions/state.hpp"
#include "../multidimensional/multiarray.hpp"

//template<typename Q> class State;
/// TODO: Save frequency grids and frequency parameters for each channel individually

#ifdef USE_MPI
#include "mpi_setup.hpp"          // mpi routines: when using mpi, only the process with ID 0 writes into file
#endif

// TODO(medium): Currently, global parameters are used to set the size of the buffer arrays.
//  Thus, in order to properly read data from a file, global parameters need to be the same as in the file
//  --> fix this: read buffer sizes (and dims) from file (-> Marc; Write this once and use for several projects!)
//  Also, always save parameters.h and FrequencyGrid.h (and whereever elso global parameters are stored)
//  for each computation and load it as well.

/// --- Constants concerning HDF5 data format --- ///

// Dataset dimensions
const int RANK_K1 = 2;
const int RANK_K2 = 2;
const int RANK_K3 = 2;
const int RANK_irreducible = 2;
const int RANK_self = 2;
const int RANK_freqs = 2;

const int N_freq_params = 64; // number of parameters for frequency grid

// Names of the individual datasets within the hdf5 file
const H5std_string	DATASET_irred("irred");

const H5std_string	DATASET_K1_a("K1_a");
const H5std_string	DATASET_K1_p("K1_p");
const H5std_string	DATASET_K1_t("K1_t");

const H5std_string	DATASET_K2_a("K2_a");
const H5std_string	DATASET_K2_p("K2_p");
const H5std_string	DATASET_K2_t("K2_t");

const H5std_string	DATASET_K3_a("K3_a");
const H5std_string	DATASET_K3_p("K3_p");
const H5std_string	DATASET_K3_t("K3_t");

const H5std_string	SELF_LIST("selflist");
const H5std_string	HARTREE("hartree");
const H5std_string	LAMBDA_LIST("lambdas");
const H5std_string  BFREQS_LISTa ("bfreqs_a");
const H5std_string  BFREQS_LISTp ("bfreqs_p");
const H5std_string  BFREQS_LISTt ("bfreqs_t");
const H5std_string  BFREQS2_LISTa("bfreqs2_a");
const H5std_string  BFREQS2_LISTp("bfreqs2_p");
const H5std_string  BFREQS2_LISTt("bfreqs2_t");
const H5std_string  BFREQS3_LISTa("bfreqs3_a");
const H5std_string  BFREQS3_LISTp("bfreqs3_p");
const H5std_string  BFREQS3_LISTt("bfreqs3_t");
const H5std_string  FFREQS_LIST ("ffreqs");
const H5std_string  FFREQS2_LISTa("ffreqs2_a");
const H5std_string  FFREQS2_LISTp("ffreqs2_p");
const H5std_string  FFREQS2_LISTt("ffreqs2_t");
const H5std_string  FFREQS3_LISTa("ffreqs3_a");
const H5std_string  FFREQS3_LISTp("ffreqs3_p");
const H5std_string  FFREQS3_LISTt("ffreqs3_t");
const H5std_string  FREQ_PARAMS("freq_params");
const H5std_string  PARAM_LIST("parameters");
const H5std_string  RE( "re" );
const H5std_string  IM( "im" );

/// --- Definitions of necessary data types --- ///

// Define struct to save complex numbers in hdf5 file
typedef struct h5_comp {
    double re; // std::real part
    double im; // imaginary part
} h5_comp;

// Create the memory data type for storing complex numbers in file
H5::CompType def_mtype_comp();

const H5::CompType mtype_comp = def_mtype_comp();
H5::DSetCreatPropList def_proplist_comp();
const H5::DSetCreatPropList plist_vert_comp = def_proplist_comp();

namespace hdf5_impl {
    /// Create new dataset of suitable datatype
    template <typename Q, typename H5Object,
            std::enable_if_t<
                    std::is_same_v<Q, double> ||
                    std::is_same_v<Q, comp>||
                    std::is_same_v<Q, int>||
                    std::is_same_v<Q, char>,
                            bool> = true>
    H5::DataSet create_Dataset(H5Object& group, const H5std_string& dataset_name, H5::DataSpace& file_space) {
    if constexpr(std::is_same_v<Q, double>) {
            H5::DataSet mydataset = group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, file_space);
            return mydataset;
    }
    else if constexpr(std::is_same_v<Q, comp>) {
        H5::DataSet mydataset = group.createDataSet(dataset_name, mtype_comp, file_space, plist_vert_comp);
        return mydataset;
    }
    else if constexpr(std::is_same_v<Q, int>) {
        H5::DataSet mydataset = group.createDataSet(dataset_name, H5::PredType::NATIVE_INT, file_space);
        return mydataset;
    }
    else if constexpr(std::is_same_v<Q, char>) {
        H5::DataSet mydataset = group.createDataSet(dataset_name, H5::PredType::NATIVE_CHAR, file_space);
        return mydataset;
    }
    }

    /// Open existing dataset of suitable datatype
    template <typename H5Object>
    H5::DataSet open_Dataset(H5Object& group, const H5std_string& dataset_name) {
            H5::DataSet mydataset = group.openDataSet(dataset_name);
            return mydataset;
    }


    template <typename container,
            std::enable_if_t<
                    std::is_same_v<typename container::value_type, double> ||
                    std::is_same_v<typename container::value_type, comp> ||
                    std::is_same_v<typename container::value_type, int> ||
                    std::is_same_v<typename container::value_type, char>,
    bool> = true>
    void write_data_to_Dataset(const container data, H5::DataSet& dataset) {
        if constexpr(std::is_same_v<typename container::value_type, double>) dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
        else if constexpr(std::is_same_v<typename container::value_type, comp>) dataset.write(data.data(), mtype_comp);
        else if constexpr(std::is_same_v<typename container::value_type, int>) dataset.write(data.data(), H5::PredType::NATIVE_INT);
        else if constexpr(std::is_same_v<typename container::value_type, char>) dataset.write(data.data(), H5::PredType::NATIVE_CHAR);
    }

    template <typename container,
            std::enable_if_t<
                    std::is_same_v<typename container::value_type, double> ||
                    std::is_same_v<typename container::value_type, comp> ||
                    std::is_same_v<typename container::value_type, int> ||
                    std::is_same_v<typename container::value_type, char>,
    bool> = true>
    void write_data_to_Dataset(const container data, H5::DataSet& dataset, H5::DataSpace& mem_space, H5::DataSpace& file_space) {
        if constexpr(std::is_same_v<typename container::value_type, double>) dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, mem_space, file_space);
        else if constexpr(std::is_same_v<typename container::value_type, comp>) dataset.write(data.data(), mtype_comp, mem_space, file_space);
        else if constexpr(std::is_same_v<typename container::value_type, int>) dataset.write(data.data(), H5::PredType::NATIVE_INT, mem_space, file_space);
        else if constexpr(std::is_same_v<typename container::value_type, char>) dataset.write(data.data(), H5::PredType::NATIVE_CHAR, mem_space, file_space);
    }


    template<typename Q, std::size_t depth, typename H5object, typename container, typename dimensions_type>
    void write_to_hdf_LambdaLayer_impl(H5object& group, const H5std_string& dataset_name, const container& data, const dimensions_type& length, const hsize_t Lambda_it, const hsize_t numberLambda_layers, const bool data_set_exists) {
        assert(Lambda_it < numberLambda_layers);

        /// create arrays with information on dimensions
        const hsize_t RANK = depth;
        hsize_t start[RANK+1];      // starting location for the hyperslab
        hsize_t stride[RANK+1];     // number of elements to separate each element or block
        hsize_t count[RANK+1];      // number of elements or blocks to select along each dimension
        hsize_t block[RANK+1];      // size of the block selected from the dataspace
        hsize_t dims_file[RANK+1];       // size of the full dataspace (including all lambda layers)
        hsize_t dims_mem[RANK];       // size of the full dataspace (including all lambda layers)
        for (int i = 0; i < RANK+1; i++) {
            stride[i] = 1;
            block[i] = 1;
        }
        count[0] = 1;           // dimension of Lambda layer
        start[0] = Lambda_it;
        dims_file[0] = numberLambda_layers;

        for (int i = 1; i < RANK+1; i++) {
            count[i] = length[i-1];
            start[i] = 0;
            dims_file[i] = length[i-1];
            dims_mem[i-1] = length[i-1];
        }

        // create dataspaces and select hyperslab
        H5::DataSpace file_space(RANK+1, dims_file);
        H5::DataSpace mem_space(RANK, dims_mem);
        file_space.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        // create or open dataset in HDF5group/file
        H5::DataSet mydataset;
        if (!data_set_exists) mydataset = hdf5_impl::create_Dataset<Q,H5object>(group, dataset_name, file_space);
        else mydataset = group.openDataSet(dataset_name);//= hdf5_impl::open_Dataset<H5object>(group, dataset_name);

        //write data to dataset
        hdf5_impl::write_data_to_Dataset(data, mydataset, mem_space, file_space);

        mydataset.close();
        mem_space.close();
        file_space.close();
    };


    template<typename Q, std::size_t depth, typename H5object>
    std::vector<Q> read_from_hdf_impl(const H5object& group, const H5std_string& dataset_name, std::array<std::size_t,depth>& dims_result) {
        H5::DataSet dataset = hdf5_impl::open_Dataset(group, dataset_name);
        H5::DataSpace file_space = dataset.getSpace();

        // get the size of the dataset
        hsize_t dims_file[depth];
        hsize_t rank = file_space.getSimpleExtentDims(dims_file, NULL);
        //std::cout<<"rank: "<<rank<<std::endl; // this is the correct number of values
        assert(rank == depth);  // Important assertion! Makes sure that the desired multiarray depth matches the data in the file


        hsize_t dims_mem[depth];
        hsize_t dims_mem_flat = 1;
        for (hsize_t i = 0; i < depth; i++) {
            dims_mem[i] = dims_file[i];
            dims_mem_flat *= dims_file[i];
        }
        H5::DataSpace dataSpace_buffer(depth, dims_mem);
        //std::cout<<"Datasize of loaded (flat): "<<dims_mem_flat<<std::endl; // this is the correct number of values


        // create a vector the same size as the dataset
        std::vector<Q> result;
        result.resize(dims_mem_flat);
        //std::cout<<"Vectsize: "<<result.size()<<std::endl;


        if constexpr(std::is_same_v<Q,double>) dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE, dataSpace_buffer, file_space);
        else if constexpr(std::is_same_v<Q,comp>) dataset.read(result.data(), mtype_comp, dataSpace_buffer, file_space);   /// Funktioniert das auch für comp?
        else if constexpr(std::is_same_v<Q,int>) dataset.read(result.data(), H5::PredType::NATIVE_INT, dataSpace_buffer, file_space);
        else if constexpr(std::is_same_v<Q,char>) dataset.read(result.data(), H5::PredType::NATIVE_CHAR, dataSpace_buffer, file_space);


        for (int i = 0; i < depth; i++) dims_result[i] = dims_file[i];
        return result;
    }


    template<typename Q, std::size_t depth, typename H5object>
    std::vector<Q> read_from_hdf_LambdaLayer_impl(const H5object& group, const H5std_string& dataset_name, std::array<std::size_t,depth>& dims_result, const int Lambda_it) {

        // open file_space to read from
        H5::DataSet dataset = hdf5_impl::open_Dataset(group, dataset_name);
        H5::DataSpace file_space = dataset.getSpace();

        // get the size of the file_space (including all Lambda layers)
        hsize_t dims_file[depth+1];       // size of the full dataspace (including all lambda layers)
        hsize_t rank = file_space.getSimpleExtentDims(dims_file, NULL);
        assert(Lambda_it < dims_file[0]); // if this fails, Lambda_it exceeded the available Lambda layers

        // create arrays with information on dimensions
        const hsize_t RANK = depth;
        hsize_t start[RANK+1];      // starting location for the hyperslab
        hsize_t stride[RANK+1];     // number of elements to separate each element or block
        hsize_t count[RANK+1];      // number of elements or blocks to select along each dimension
        hsize_t block[RANK+1];      // size of the block selected from the dataspace
        hsize_t dims_mem[RANK];       // size of the full dataspace (including all lambda layers)
        for (int i = 0; i < RANK+1; i++) {
            stride[i] = 1;
            block[i] = 1;
        }
        count[0] = 1;           // dimension of Lambda layer
        start[0] = Lambda_it;
        for (int i = 1; i < RANK+1; i++) {
            count[i] = dims_file[i];
            start[i] = 0;
            dims_result[i-1] = dims_file[i];
            dims_mem[i-1] = dims_file[i];
        }
        hsize_t dims_mem_flat = 1;
        for (hsize_t i = 0; i < depth; i++) {
            dims_mem_flat *= dims_file[i+1];
        }

        assert(rank == depth+1);  // Important assertion! Makes sure that the desired multiarray depth matches the data in the file

        // select hyperslab to read from
        file_space.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        H5::DataSet mydataset;

        H5::DataSpace dataSpace_buffer(depth, dims_mem);

        // create a vector the same size as the dataset
        std::vector<Q> result(dims_mem_flat);


        if constexpr(std::is_same_v<Q,double>) dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE, dataSpace_buffer, file_space);
        else if constexpr(std::is_same_v<Q,comp>) dataset.read(result.data(), mtype_comp, dataSpace_buffer, file_space);   /// Funktioniert das auch für comp?
        else if constexpr(std::is_same_v<Q,int>) dataset.read(result.data(), H5::PredType::NATIVE_INT, dataSpace_buffer, file_space);
        else if constexpr(std::is_same_v<Q,char>) dataset.read(result.data(), H5::PredType::NATIVE_CHAR, dataSpace_buffer, file_space);

        return result;
    }


}


/// Write scalar to HDF group/file as Attribute
template<typename Q, typename H5object,
        std::enable_if_t<
                                std::is_same_v<Q, double> ||
                                std::is_same_v<Q, comp> ||
                                std::is_same_v<Q, int> ||
                                std::is_same_v<Q, char>,
                bool> = true>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const Q& data, const bool data_set_exists) {
    H5::DataSpace file_space(H5S_SCALAR);

    if constexpr(std::is_same_v<Q, double>) {
        H5::Attribute attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_DOUBLE, file_space );
        attr.write( H5::PredType::NATIVE_DOUBLE, &data );
    }
    else if constexpr(std::is_same_v<Q, comp>) {
        H5::Attribute attr = group.createAttribute( dataset_name, mtype_comp, file_space );
        attr.write( mtype_comp, &data );
    }
    else if constexpr(std::is_same_v<Q, int>) {
        H5::Attribute attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_INT, file_space );
        attr.write( H5::PredType::NATIVE_INT, &data );
    }
    else if constexpr(std::is_same_v<Q, char>) {
        H5::Attribute attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_CHAR, file_space );
        attr.write( H5::PredType::NATIVE_CHAR, &data );
    }

    file_space.close();
};
/// Write vector to HDF group/file
template <typename Q, typename H5object>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const std::vector<Q>& data, const bool data_set_exists) {
    hsize_t dims[1] = {data.size()};
    H5::DataSpace file_space(1, dims);
    H5::DataSet mydataset;
    if (!data_set_exists) mydataset = hdf5_impl::create_Dataset<Q,H5object>(group, dataset_name, file_space);
    else mydataset = group.openDataSet(dataset_name); // hdf5_impl::open_Dataset<H5object>(group, dataset_name);

    hdf5_impl::write_data_to_Dataset(data, mydataset);

    mydataset.close();
    file_space.close();
}
/// Write multiarray to HDF group/file
template<typename Q, std::size_t depth, typename H5object>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const multidimensional::multiarray<Q, depth>& data, const bool data_set_exists) {
    hsize_t dims[depth];
    std::copy(std::begin(data.length()), std::end(data.length()), std::begin(dims));
    H5::DataSpace file_space(depth, dims);

    H5::DataSet mydataset;
    if (!data_set_exists) mydataset = hdf5_impl::create_Dataset<Q,H5object>(group, dataset_name, file_space);
    else mydataset = group.openDataSet(dataset_name); // hdf5_impl::open_Dataset<H5object>(group, dataset_name);
    hdf5_impl::write_data_to_Dataset(data, mydataset);

    mydataset.close();
    file_space.close();
};

/// Write vector to Lambda layer of HDF group/file
template <typename Q, typename H5object>
void write_to_hdf_LambdaLayer(H5object& group, const H5std_string& dataset_name, const std::vector<Q>& data, const hsize_t Lambda_it, const hsize_t numberLambda_layers, const bool data_set_exists) {
    assert(Lambda_it < numberLambda_layers);
    std::array<std::size_t,1> length = {data.size()};
    hdf5_impl::write_to_hdf_LambdaLayer_impl<Q,1,H5object>(group, dataset_name, data, length, Lambda_it, numberLambda_layers, data_set_exists);
}
/// Write multiarray to Lambda layer of HDF group/file
template<typename Q, std::size_t depth, typename H5object>
void write_to_hdf_LambdaLayer(H5object& group, const H5std_string& dataset_name, const multidimensional::multiarray<Q, depth>& data, const hsize_t Lambda_it, const hsize_t numberLambda_layers, const bool data_set_exists) {
    assert(Lambda_it < numberLambda_layers);
    typename multidimensional::multiarray<Q, depth>::dimensions_type length = data.length();
    hdf5_impl::write_to_hdf_LambdaLayer_impl<Q,depth,H5object>(group, dataset_name, data, length, Lambda_it, numberLambda_layers, data_set_exists);
};





/// Read multiarray from HDF group/file
template<typename Q, std::size_t depth, typename H5object>
void read_from_hdf(const H5object& group, const H5std_string& dataset_name, multidimensional::multiarray<Q,depth>& result) {

    typename multidimensional::multiarray<Q,depth>::dimensions_type dims_result;
    std::vector<Q> result_vec = hdf5_impl::read_from_hdf_impl<Q,depth,H5object>(group, dataset_name, dims_result);
    result = multidimensional::multiarray<Q,depth>(dims_result, result_vec);
}
/// Read vector from HDF group/file
template<typename Q, typename H5object>
void read_from_hdf(const H5object& group, const H5std_string& dataset_name, std::vector<Q>& result) {
    std::array<std::size_t, 1> dims_result;
    result = hdf5_impl::read_from_hdf_impl<Q,1,H5object>(group, dataset_name, dims_result);
}


/// Read multiarray from Lambda layer of HDF group/file
template<typename Q, std::size_t depth, typename H5object>
void read_from_hdf_LambdaLayer(const H5object& group, const H5std_string& dataset_name, multidimensional::multiarray<Q,depth>& result, const int Lambda_it) {

    typename multidimensional::multiarray<Q,depth>::dimensions_type dims_result;
    std::vector<Q> result_vec = hdf5_impl::read_from_hdf_LambdaLayer_impl<Q,depth,H5object>(group, dataset_name, dims_result, Lambda_it);
    result = multidimensional::multiarray<Q,depth>(dims_result, result_vec);
}
/// Read vector from Lambda layer of HDF group/file
template<typename Q, typename H5object>
void read_from_hdf_LambdaLayer(const H5object& group, const H5std_string& dataset_name, std::vector<Q>& result, const int Lambda_it) {
    std::array<std::size_t, 1> dims_result;
    result = hdf5_impl::read_from_hdf_LambdaLayer_impl<Q,1,H5object>(group, dataset_name, dims_result, Lambda_it);
}


namespace hdf5_impl {

    void write_freqparams_to_hdf_LambdaLayer(H5::Group& group, const FrequencyGrid& freqgrid, int Lambda_it, int numberLambdaLayers, bool file_exists, bool verbose=true);
    void init_freqgrid_from_hdf_LambdaLayer(H5::Group& group, FrequencyGrid& freqgrid, int Lambda_it, double Lambda);

    template<typename Q>
    void write_state_to_hdf_LambdaLayer(const H5std_string& filename, const State<Q>& state, const int Lambda_it, const int numberLambdaLayers, const bool file_exists, const bool verbose=true) {
        H5::H5File file_out;
        if (!file_exists) {
            // If file doesn't exist, create a new file using the default property lists.
            file_out = H5::H5File(filename, H5F_ACC_TRUNC);
        }
        else {
            // If file exists, open existing file. Access rights: read/write
            file_out = H5::H5File(filename, H5F_ACC_RDWR);
        }

        /// Save state data
        write_to_hdf_LambdaLayer<double>(file_out, LAMBDA_LIST, std::vector<double>({state.Lambda}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, SELF_LIST, state.selfenergy.Sigma, Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, HARTREE, std::vector<Q>({state.selfenergy.asymp_val_R}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_irred, state.vertex.irred().get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K1_a, state.vertex.avertex().K1.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K1_p, state.vertex.pvertex().K1.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K1_t, state.vertex.tvertex().K1.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS_LISTa, state.vertex.avertex().K1.K1_get_freqGrid().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS_LISTp, state.vertex.pvertex().K1.K1_get_freqGrid().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS_LISTt, state.vertex.tvertex().K1.K1_get_freqGrid().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
#if MAX_DIAG_CLASS>1
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2_a, state.vertex.avertex().K2.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2_p, state.vertex.pvertex().K2.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2_t, state.vertex.tvertex().K2.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2_LISTa, state.vertex.avertex().K2.K2_get_freqGrid_b().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2_LISTp, state.vertex.pvertex().K2.K2_get_freqGrid_b().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2_LISTt, state.vertex.tvertex().K2.K2_get_freqGrid_b().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2_LISTa, state.vertex.avertex().K2.K2_get_freqGrid_f().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2_LISTp, state.vertex.pvertex().K2.K2_get_freqGrid_f().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2_LISTt, state.vertex.tvertex().K2.K2_get_freqGrid_f().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
#endif
#if MAX_DIAG_CLASS>2
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K3_a, state.vertex.avertex().K3.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K3_p, state.vertex.pvertex().K3.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K3_t, state.vertex.tvertex().K3.get_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS3_LISTa, state.vertex.avertex().K3.K3_get_freqGrid_b().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS3_LISTp, state.vertex.pvertex().K3.K3_get_freqGrid_b().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS3_LISTt, state.vertex.tvertex().K3.K3_get_freqGrid_b().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTa, state.vertex.avertex().K3.K3_get_freqGrid_f().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTp, state.vertex.pvertex().K3.K3_get_freqGrid_f().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTt, state.vertex.tvertex().K3.K3_get_freqGrid_f().get_ws_vec(), Lambda_it, numberLambdaLayers, file_exists);
#endif
        if(!file_exists) {
            /// Write used parameters for documentation purpose
            H5::Group group_params(file_out.createGroup(PARAM_LIST));
            write_to_hdf(group_params, "REG", REG, file_exists);
            write_to_hdf(group_params, "Gamma", glb_Gamma, file_exists);
            write_to_hdf(group_params, "MAX_DIAG_CLASS", MAX_DIAG_CLASS, file_exists);
            write_to_hdf(group_params, "N_LOOPS", N_LOOPS, file_exists);
            write_to_hdf(group_params, "T", glb_T, file_exists);
            write_to_hdf(group_params, "mu", glb_mu, file_exists);
            write_to_hdf(group_params, "U", glb_U, file_exists);
            write_to_hdf(group_params, "epsilon", glb_epsilon, file_exists);
            write_to_hdf(group_params, "V", glb_V, file_exists);
            write_to_hdf(group_params, "ODEsolver", ODEsolver, file_exists);
        }

        H5::Group group_freqparams;
        H5::Group group_freqparams_bfreqsa;
        H5::Group group_freqparams_bfreqsp;
        H5::Group group_freqparams_bfreqst;
        H5::Group group_freqparams_bfreqs2a;
        H5::Group group_freqparams_bfreqs2p;
        H5::Group group_freqparams_bfreqs2t;
        H5::Group group_freqparams_bfreqs3a;
        H5::Group group_freqparams_bfreqs3p;
        H5::Group group_freqparams_bfreqs3t;
        H5::Group group_freqparams_ffreqs2a;
        H5::Group group_freqparams_ffreqs2p;
        H5::Group group_freqparams_ffreqs2t;
        H5::Group group_freqparams_ffreqs3a;
        H5::Group group_freqparams_ffreqs3p;
        H5::Group group_freqparams_ffreqs3t;
    if (!file_exists) {
        group_freqparams = file_out.createGroup(FREQ_PARAMS);
        group_freqparams_bfreqsa = group_freqparams.createGroup(BFREQS_LISTa );
        group_freqparams_bfreqsp = group_freqparams.createGroup(BFREQS_LISTp );
        group_freqparams_bfreqst = group_freqparams.createGroup(BFREQS_LISTt );
        group_freqparams_bfreqs2a = group_freqparams.createGroup(BFREQS2_LISTa);
        group_freqparams_bfreqs2p = group_freqparams.createGroup(BFREQS2_LISTp);
        group_freqparams_bfreqs2t = group_freqparams.createGroup(BFREQS2_LISTt);
        group_freqparams_bfreqs3a = group_freqparams.createGroup(BFREQS3_LISTa);
        group_freqparams_bfreqs3p = group_freqparams.createGroup(BFREQS3_LISTp);
        group_freqparams_bfreqs3t = group_freqparams.createGroup(BFREQS3_LISTt);
        group_freqparams_ffreqs2a = group_freqparams.createGroup(FFREQS2_LISTa);
        group_freqparams_ffreqs2p = group_freqparams.createGroup(FFREQS2_LISTp);
        group_freqparams_ffreqs2t = group_freqparams.createGroup(FFREQS2_LISTt);
        group_freqparams_ffreqs3a = group_freqparams.createGroup(FFREQS3_LISTa);
        group_freqparams_ffreqs3p = group_freqparams.createGroup(FFREQS3_LISTp);
        group_freqparams_ffreqs3t = group_freqparams.createGroup(FFREQS3_LISTt);
    }
    else {
        group_freqparams = file_out.openGroup(FREQ_PARAMS);
        group_freqparams_bfreqsa = group_freqparams.openGroup(BFREQS_LISTa );
        group_freqparams_bfreqsp = group_freqparams.openGroup(BFREQS_LISTp );
        group_freqparams_bfreqst = group_freqparams.openGroup(BFREQS_LISTt );
        group_freqparams_bfreqs2a = group_freqparams.openGroup(BFREQS2_LISTa);
        group_freqparams_bfreqs2p = group_freqparams.openGroup(BFREQS2_LISTp);
        group_freqparams_bfreqs2t = group_freqparams.openGroup(BFREQS2_LISTt);
        group_freqparams_bfreqs3a = group_freqparams.openGroup(BFREQS3_LISTa);
        group_freqparams_bfreqs3p = group_freqparams.openGroup(BFREQS3_LISTp);
        group_freqparams_bfreqs3t = group_freqparams.openGroup(BFREQS3_LISTt);
        group_freqparams_ffreqs2a = group_freqparams.openGroup(FFREQS2_LISTa);
        group_freqparams_ffreqs2p = group_freqparams.openGroup(FFREQS2_LISTp);
        group_freqparams_ffreqs2t = group_freqparams.openGroup(FFREQS2_LISTt);
        group_freqparams_ffreqs3a = group_freqparams.openGroup(FFREQS3_LISTa);
        group_freqparams_ffreqs3p = group_freqparams.openGroup(FFREQS3_LISTp);
        group_freqparams_ffreqs3t = group_freqparams.openGroup(FFREQS3_LISTt);

    }
        /// Write frequency parameters
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqsa, state.vertex.avertex().K1.K1_get_freqGrid()  , Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqsp, state.vertex.pvertex().K1.K1_get_freqGrid()  , Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqst, state.vertex.tvertex().K1.K1_get_freqGrid()  , Lambda_it, numberLambdaLayers, file_exists, verbose);
#if MAX_DIAG_CLASS>1
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2a,state.vertex.avertex().K2.K2_get_freqGrid_b(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2p,state.vertex.pvertex().K2.K2_get_freqGrid_b(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2t,state.vertex.tvertex().K2.K2_get_freqGrid_b(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2a,state.vertex.avertex().K2.K2_get_freqGrid_f(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2p,state.vertex.pvertex().K2.K2_get_freqGrid_f(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2t,state.vertex.tvertex().K2.K2_get_freqGrid_f(), Lambda_it, numberLambdaLayers, file_exists, verbose);
#endif
#if MAX_DIAG_CLASS>2
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs3a,state.vertex.avertex().K3.K3_get_freqGrid_b(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs3p,state.vertex.pvertex().K3.K3_get_freqGrid_b(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs3t,state.vertex.tvertex().K3.K3_get_freqGrid_b(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3a,state.vertex.avertex().K3.K3_get_freqGrid_f(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3p,state.vertex.pvertex().K3.K3_get_freqGrid_f(), Lambda_it, numberLambdaLayers, file_exists, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3t,state.vertex.tvertex().K3.K3_get_freqGrid_f(), Lambda_it, numberLambdaLayers, file_exists, verbose);
#endif

        file_out.close();
    }

}

/// Create file with fixed number of Lambda layers and save state to first Lambda layer
template <typename Q>
void write_state_to_hdf(const H5std_string FILE_NAME, const State<Q>& state_in, const int Lambda_size, const bool verbose=true) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        hdf5_impl::write_state_to_hdf_LambdaLayer(FILE_NAME, state_in, 0, Lambda_size, false);
        if (verbose) {
            print("Successfully saved in hdf5 file: ", FILE_NAME);
            print_add(" in Lambda-layer ", 0, false);
            print_add("", true);
        }
    }
}

/// Open file and save state to a specified Lambda layer
template <typename Q>
void add_state_to_hdf(const H5std_string FILE_NAME, const Q& state_in, int Lambda_it, const bool verbose=true) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        multidimensional::multiarray<double,2> Lambdas;
        H5::H5File file_out(FILE_NAME, H5F_ACC_RDONLY);
        read_from_hdf<double>(file_out, LAMBDA_LIST, Lambdas);
        file_out.close();
        const std::size_t Lambda_size = Lambdas.size();

        if (Lambda_it < Lambda_size) {
            hdf5_impl::write_state_to_hdf_LambdaLayer(FILE_NAME, state_in, Lambda_it, Lambda_size, true);
            if (verbose) {
                print("Successfully saved in hdf5 file: ", FILE_NAME);
                print_add(" in Lambda-layer ", Lambda_it, false);
                print_add("", true);
            }
        } else {
            print("\t\t  ERROR: Cannot write to file ", FILE_NAME, " since Lambda layer", Lambda_it,
                  " is out of range.", true);
        }
    }
}

/// Read state from specified Lambda layer of hdf file
State<state_datatype> read_state_from_hdf_LambdaLayer(const H5std_string& filename, int Lambda_it) ;



/**
 * Class containing buffer lengths and arrays to buffer data for selfenergy and irreducible vertex
 * as well as K1, K2, K3 classes.
 * Constructor creates new arrays, destructor deletes them.
 */
class Buffer {
public:
    double * lambda;
    double * freq_params;
    double * bfreqsa_buffer;
    double * bfreqsp_buffer;
    double * bfreqst_buffer;
    double * ffreqs_buffer;
    const size_t self_dim = dimsSE_flat;    // length of self-energy buffer
#ifdef KELDYSH_FORMALISM
    const int irred_dim = 16 * n_in;                                  // length of irreducible vertex buffer
#else
    const int irred_dim = n_in;                                       // length of irreducible vertex buffer
#endif
    h5_comp * selfenergy;
    h5_comp * irreducible_class;

#if MAX_DIAG_CLASS >= 1
    const size_t K1_dim = dimsK1_flat;                         // length of K1 buffer
    h5_comp * K1_class_a;
    h5_comp * K1_class_p;
    h5_comp * K1_class_t;
#endif
#if MAX_DIAG_CLASS >= 2
    double * bfreqs2a_buffer;
    double * ffreqs2a_buffer;
    double * bfreqs2p_buffer;
    double * ffreqs2p_buffer;
    double * bfreqs2t_buffer;
    double * ffreqs2t_buffer;
    const size_t K2_dim = dimsK2_flat;               // length of K2 buffer
    h5_comp * K2_class_a;
    h5_comp * K2_class_p;
    h5_comp * K2_class_t;
#endif
#if MAX_DIAG_CLASS >= 3
    double * bfreqs3a_buffer;
    double * ffreqs3a_buffer;
    double * bfreqs3p_buffer;
    double * ffreqs3p_buffer;
    double * bfreqs3t_buffer;
    double * ffreqs3t_buffer;
    const size_t K3_dim = dimsK3_flat;    // length of K3 buffer
    h5_comp * K3_class_a;
    h5_comp * K3_class_p;
    h5_comp * K3_class_t;
#endif

    Buffer() {
        lambda = new double [1];
        freq_params = new double[N_freq_params];
        bfreqsa_buffer = new double[nBOS];                        // create buffer for bosonic frequencies
        bfreqsp_buffer = new double[nBOS];                        // create buffer for bosonic frequencies
        bfreqst_buffer = new double[nBOS];                        // create buffer for bosonic frequencies
        ffreqs_buffer  = new double[nFER];                        // create buffer for fermionic frequencies
        selfenergy = new h5_comp[self_dim];                           // create buffer for self-energy
        irreducible_class = new h5_comp[irred_dim];              // create buffer for irreducible vertex
#if MAX_DIAG_CLASS >= 1
        K1_class_a = new h5_comp[K1_dim];                        // create buffer for K1_a
        K1_class_p = new h5_comp[K1_dim];                        // create buffer for K1_p
        K1_class_t = new h5_comp[K1_dim];                        // create buffer for K1_t
#endif
#if MAX_DIAG_CLASS >= 2
        bfreqs2a_buffer = new double[nBOS2];                      // create buffer for bosonic frequencies for K2
        ffreqs2a_buffer = new double[nFER2];                      // create buffer for fermionic frequencies for K2
        bfreqs2p_buffer = new double[nBOS2];                      // create buffer for bosonic frequencies for K2
        ffreqs2p_buffer = new double[nFER2];                      // create buffer for fermionic frequencies for K2
        bfreqs2t_buffer = new double[nBOS2];                      // create buffer for bosonic frequencies for K2
        ffreqs2t_buffer = new double[nFER2];                      // create buffer for fermionic frequencies for K2
        K2_class_a = new h5_comp[K2_dim];                        // create buffer for K2_a
        K2_class_p = new h5_comp[K2_dim];                        // create buffer for K2_p
        K2_class_t = new h5_comp[K2_dim];                        // create buffer for K2_t
#endif
#if MAX_DIAG_CLASS >= 3
        bfreqs3a_buffer = new double[nBOS3];                      // create buffer for bosonic frequencies for K3
        ffreqs3a_buffer = new double[nFER3];                      // create buffer for fermionic frequencies for K3
        bfreqs3p_buffer = new double[nBOS3];                      // create buffer for bosonic frequencies for K3
        ffreqs3p_buffer = new double[nFER3];                      // create buffer for fermionic frequencies for K3
        bfreqs3t_buffer = new double[nBOS3];                      // create buffer for bosonic frequencies for K3
        ffreqs3t_buffer = new double[nFER3];                      // create buffer for fermionic frequencies for K3
        K3_class_a = new h5_comp[K3_dim];                        // create buffer for K3_a
        K3_class_p = new h5_comp[K3_dim];                        // create buffer for K3_p
        K3_class_t = new h5_comp[K3_dim];                        // create buffer for K3_t
#endif
    }

    ~Buffer() {
        delete[] freq_params;
        delete[] bfreqsa_buffer;
        delete[] bfreqsp_buffer;
        delete[] bfreqst_buffer;
        delete[] ffreqs_buffer;
        delete[] selfenergy;
        delete[] irreducible_class;
#if MAX_DIAG_CLASS >= 1
        delete[] K1_class_a;
        delete[] K1_class_p;
        delete[] K1_class_t;
#endif
#if MAX_DIAG_CLASS >= 2
        delete[] bfreqs2a_buffer;
        delete[] ffreqs2a_buffer;
        delete[] bfreqs2p_buffer;
        delete[] ffreqs2p_buffer;
        delete[] bfreqs2t_buffer;
        delete[] ffreqs2t_buffer;
        delete[] K2_class_a;
        delete[] K2_class_p;
        delete[] K2_class_t;
#endif
#if MAX_DIAG_CLASS >= 3
        delete[] bfreqs3a_buffer;
        delete[] ffreqs3a_buffer;
        delete[] bfreqs3p_buffer;
        delete[] ffreqs3p_buffer;
        delete[] bfreqs3t_buffer;
        delete[] ffreqs3t_buffer;
        delete[] K3_class_a;
        delete[] K3_class_p;
        delete[] K3_class_t;
#endif
    }

template <typename Q>
    void initialize(const State<Q>& state_in) {
        //print("Starting to copy to buffer...", true);
        FrequencyGrid bfreqsa = state_in.vertex.avertex().K1.K1_get_freqGrid();
        FrequencyGrid bfreqsp = state_in.vertex.pvertex().K1.K1_get_freqGrid();
        FrequencyGrid bfreqst = state_in.vertex.tvertex().K1.K1_get_freqGrid();
        FrequencyGrid ffreqs = state_in.selfenergy.frequencies;
        freq_params[0] = (double) bfreqsa.N_w;
        freq_params[1] = bfreqsa.w_upper;
        freq_params[2] = bfreqsa.w_lower;
        freq_params[3] = bfreqsa.W_scale;
        freq_params[4] = (double) ffreqs.N_w;
        freq_params[5] = ffreqs.w_upper;
        freq_params[6] = ffreqs.w_lower;
        freq_params[7] = ffreqs.W_scale;
        freq_params[24] = (double) bfreqsp.N_w;
        freq_params[25] = bfreqsp.w_upper;
        freq_params[26] = bfreqsp.w_lower;
        freq_params[27] = bfreqsp.W_scale;
        freq_params[28] = (double) bfreqst.N_w;
        freq_params[29] = bfreqst.w_upper;
        freq_params[30] = bfreqst.w_lower;
        freq_params[31] = bfreqst.W_scale;

    convert_vec_to_type(bfreqsa.get_ws_vec(), bfreqsa_buffer);
    convert_vec_to_type(bfreqsp.get_ws_vec(), bfreqsp_buffer);
    convert_vec_to_type(bfreqst.get_ws_vec(), bfreqst_buffer);
    convert_vec_to_type(ffreqs.get_ws_vec(), ffreqs_buffer);

        for (int i=0; i<self_dim; ++i) {                        // write self-energy into buffer
#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM) and not defined(HUBBARD)
            // in the particle-hole symmetric case in Matsubara we only save the imaginary part of the selfenergy
            selfenergy[i].re = glb_U/2.;
            selfenergy[i].im = state_in.selfenergy.acc(i);
#else
            selfenergy[i].re = std::real(state_in.selfenergy.acc(i));
            selfenergy[i].im = std::imag(state_in.selfenergy.acc(i));
#endif
        }
        for (int i=0; i<irred_dim; ++i) {                       // write irreducible vertex into buffer
            irreducible_class[i].re = std::real(state_in.vertex.irred().acc(i));
            irreducible_class[i].im = std::imag(state_in.vertex.irred().acc(i));
        }
#if MAX_DIAG_CLASS >= 1
        for(int i=0; i<K1_dim; ++i) {                                // write K1 into buffer
            K1_class_a[i].re = std::real(state_in.vertex.avertex().K1.acc(i));
            K1_class_a[i].im = std::imag(state_in.vertex.avertex().K1.acc(i));

            K1_class_p[i].re = std::real(state_in.vertex.pvertex().K1.acc(i));
            K1_class_p[i].im = std::imag(state_in.vertex.pvertex().K1.acc(i));

            K1_class_t[i].re = std::real(state_in.vertex.tvertex().K1.acc(i));
            K1_class_t[i].im = std::imag(state_in.vertex.tvertex().K1.acc(i));
        }
#endif
#if MAX_DIAG_CLASS >= 2
    FrequencyGrid bfreqs2a = state_in.vertex.avertex().K2.K2_get_freqGrid_b();
    FrequencyGrid ffreqs2a = state_in.vertex.avertex().K2.K2_get_freqGrid_f();
    FrequencyGrid bfreqs2p = state_in.vertex.pvertex().K2.K2_get_freqGrid_b();
    FrequencyGrid ffreqs2p = state_in.vertex.pvertex().K2.K2_get_freqGrid_f();
    FrequencyGrid bfreqs2t = state_in.vertex.tvertex().K2.K2_get_freqGrid_b();
    FrequencyGrid ffreqs2t = state_in.vertex.tvertex().K2.K2_get_freqGrid_f();
    freq_params[8]  = (double) bfreqs2a.N_w;
    freq_params[9]  = bfreqs2a.w_upper;
    freq_params[10] = bfreqs2a.w_lower;
    freq_params[11] = bfreqs2a.W_scale;
    freq_params[12] = (double) ffreqs2a.N_w;
    freq_params[13] = ffreqs2a.w_upper;
    freq_params[14] = ffreqs2a.w_lower;
    freq_params[15] = ffreqs2a.W_scale;
    freq_params[32]  = (double) bfreqs2p.N_w;
    freq_params[33]  = bfreqs2p.w_upper;
    freq_params[34] = bfreqs2p.w_lower;
    freq_params[35] = bfreqs2p.W_scale;
    freq_params[36] = (double) ffreqs2p.N_w;
    freq_params[37] = ffreqs2p.w_upper;
    freq_params[38] = ffreqs2p.w_lower;
    freq_params[39] = ffreqs2p.W_scale;
    freq_params[40]  = (double) bfreqs2t.N_w;
    freq_params[41]  = bfreqs2t.w_upper;
    freq_params[42] = bfreqs2t.w_lower;
    freq_params[43] = bfreqs2t.W_scale;
    freq_params[44] = (double) ffreqs2t.N_w;
    freq_params[45] = ffreqs2t.w_upper;
    freq_params[46] = ffreqs2t.w_lower;
    freq_params[47] = ffreqs2t.W_scale;

    convert_vec_to_type(bfreqs2a.get_ws_vec(), bfreqs2a_buffer);
    convert_vec_to_type(ffreqs2a.get_ws_vec(), ffreqs2a_buffer);
    convert_vec_to_type(bfreqs2p.get_ws_vec(), bfreqs2p_buffer);
    convert_vec_to_type(ffreqs2p.get_ws_vec(), ffreqs2p_buffer);
    convert_vec_to_type(bfreqs2t.get_ws_vec(), bfreqs2t_buffer);
    convert_vec_to_type(ffreqs2t.get_ws_vec(), ffreqs2t_buffer);

        for(int i=0; i<K2_dim; ++i) {                                // write K2 into buffer
            K2_class_a[i].re = std::real(state_in.vertex.avertex().K2.acc(i));
            K2_class_a[i].im = std::imag(state_in.vertex.avertex().K2.acc(i));

            K2_class_p[i].re = std::real(state_in.vertex.pvertex().K2.acc(i));
            K2_class_p[i].im = std::imag(state_in.vertex.pvertex().K2.acc(i));

            K2_class_t[i].re = std::real(state_in.vertex.tvertex().K2.acc(i));
            K2_class_t[i].im = std::imag(state_in.vertex.tvertex().K2.acc(i));
        }
#endif
#if MAX_DIAG_CLASS >= 3
        FrequencyGrid bfreqs3a = state_in.vertex.avertex().K3.K3_get_freqGrid_b();
        FrequencyGrid ffreqs3a = state_in.vertex.avertex().K3.K3_get_freqGrid_f();
        FrequencyGrid bfreqs3p = state_in.vertex.pvertex().K3.K3_get_freqGrid_b();
        FrequencyGrid ffreqs3p = state_in.vertex.pvertex().K3.K3_get_freqGrid_f();
        FrequencyGrid bfreqs3t = state_in.vertex.tvertex().K3.K3_get_freqGrid_b();
        FrequencyGrid ffreqs3t = state_in.vertex.tvertex().K3.K3_get_freqGrid_f();
        freq_params[16] = (double) bfreqs3a.N_w;
        freq_params[17] = bfreqs3a.w_upper;
        freq_params[18] = bfreqs3a.w_lower;
        freq_params[19] = bfreqs3a.W_scale;
        freq_params[20] = (double) ffreqs3a.N_w;
        freq_params[21] = ffreqs3a.w_upper;
        freq_params[22] = ffreqs3a.w_lower;
        freq_params[23] = ffreqs3a.W_scale;
        freq_params[48] = (double) bfreqs3p.N_w;
        freq_params[49] = bfreqs3p.w_upper;
        freq_params[50] = bfreqs3p.w_lower;
        freq_params[51] = bfreqs3p.W_scale;
        freq_params[52] = (double) ffreqs3p.N_w;
        freq_params[53] = ffreqs3p.w_upper;
        freq_params[54] = ffreqs3p.w_lower;
        freq_params[55] = ffreqs3p.W_scale;
        freq_params[56] = (double) bfreqs3t.N_w;
        freq_params[57] = bfreqs3t.w_upper;
        freq_params[58] = bfreqs3t.w_lower;
        freq_params[59] = bfreqs3t.W_scale;
        freq_params[60] = (double) ffreqs3t.N_w;
        freq_params[61] = ffreqs3t.w_upper;
        freq_params[62] = ffreqs3t.w_lower;
        freq_params[63] = ffreqs3t.W_scale;

    convert_vec_to_type(bfreqs3a.get_ws_vec(), bfreqs3a_buffer);
    convert_vec_to_type(ffreqs3a.get_ws_vec(), ffreqs3a_buffer);
    convert_vec_to_type(bfreqs3p.get_ws_vec(), bfreqs3p_buffer);
    convert_vec_to_type(ffreqs3p.get_ws_vec(), ffreqs3p_buffer);
    convert_vec_to_type(bfreqs3t.get_ws_vec(), bfreqs3t_buffer);
    convert_vec_to_type(ffreqs3t.get_ws_vec(), ffreqs3t_buffer);

        for(int i=0; i<K3_dim; ++i) {                                // write K3 into buffer
            K3_class_a[i].re = std::real(state_in.vertex.avertex().K3.acc(i));
            K3_class_a[i].im = std::imag(state_in.vertex.avertex().K3.acc(i));

            K3_class_p[i].re = std::real(state_in.vertex.pvertex().K3.acc(i));
            K3_class_p[i].im = std::imag(state_in.vertex.pvertex().K3.acc(i));

            K3_class_t[i].re = std::real(state_in.vertex.tvertex().K3.acc(i));
            K3_class_t[i].im = std::imag(state_in.vertex.tvertex().K3.acc(i));
        }
#endif
        //print("Buffer ready. Preparing for saving into Hdf5 file...", true);
    }
};

// Wrapper for static cast from int to hsize_t, used in class Dims
hsize_t h5_cast(int dim);

/**
 * Class containing dimension arrays for data sets in file and for the buffers.
 * Constructor initializes them to fixed values specified via buffer lengths.
 */
class Dims {
public:
    hsize_t Lambda[1];
    hsize_t freq_params_dims[2];
    hsize_t freq_params_buffer_dims[2];
    hsize_t bfreqs_dims[2];
    hsize_t ffreqs_dims[2];
    hsize_t bfreqs_buffer_dims[1];
    hsize_t ffreqs_buffer_dims[1];
    hsize_t params[1];
    hsize_t selfenergy[2];
    hsize_t selfenergy_buffer[1];
    hsize_t irreducible[2];
    hsize_t irreducible_buffer[1];
#if MAX_DIAG_CLASS >= 1
    hsize_t K1[2];
    hsize_t K1_buffer[1];
#endif
#if MAX_DIAG_CLASS >= 2
    hsize_t bfreqs2_dims[2];
    hsize_t ffreqs2_dims[2];
    hsize_t bfreqs2_buffer_dims[1];
    hsize_t ffreqs2_buffer_dims[1];
    hsize_t K2[2];
    hsize_t K2_buffer[1];
#endif
#if MAX_DIAG_CLASS >= 3
    hsize_t bfreqs3_dims[2];
    hsize_t ffreqs3_dims[2];
    hsize_t bfreqs3_buffer_dims[1];
    hsize_t ffreqs3_buffer_dims[1];
    hsize_t K3[2];
    hsize_t K3_buffer[1];
#endif

    Dims (Buffer& buffer, long Lambda_size) :
#if MAX_DIAG_CLASS >= 1
        K1 {h5_cast(Lambda_size), h5_cast(buffer.K1_dim)},
        K1_buffer {h5_cast(buffer.K1_dim)},
#endif
#if MAX_DIAG_CLASS >= 2
        bfreqs2_dims {h5_cast(Lambda_size), h5_cast(nBOS2)},
        ffreqs2_dims {h5_cast(Lambda_size), h5_cast(nFER2)},
        bfreqs2_buffer_dims {h5_cast(nBOS2)},
        ffreqs2_buffer_dims {h5_cast(nFER2)},
        K2 {h5_cast(Lambda_size), h5_cast(buffer.K2_dim)},
        K2_buffer {h5_cast(buffer.K2_dim)},
#endif
#if MAX_DIAG_CLASS >= 3
        bfreqs3_dims {h5_cast(Lambda_size), h5_cast(nBOS3)},
        ffreqs3_dims {h5_cast(Lambda_size), h5_cast(nFER3)},
        bfreqs3_buffer_dims {h5_cast(nBOS3)},
        ffreqs3_buffer_dims {h5_cast(nFER3)},
        K3 {h5_cast(Lambda_size), h5_cast(buffer.K3_dim)},
        K3_buffer {h5_cast(buffer.K3_dim)},
#endif
        Lambda {h5_cast(Lambda_size)},
        freq_params_dims {h5_cast(Lambda_size), h5_cast(N_freq_params)},
        freq_params_buffer_dims {h5_cast(N_freq_params)},
        bfreqs_dims {h5_cast(Lambda_size), h5_cast(nBOS)},
        ffreqs_dims {h5_cast(Lambda_size), h5_cast(nFER)},
        bfreqs_buffer_dims {h5_cast(nBOS)},
        ffreqs_buffer_dims {h5_cast(nFER)},
        params {h5_cast(param_size)},
        selfenergy {h5_cast(Lambda_size), h5_cast(buffer.self_dim)},
        selfenergy_buffer {h5_cast(buffer.self_dim)},
        irreducible {h5_cast(Lambda_size), h5_cast(buffer.irred_dim)},
        irreducible_buffer {h5_cast(buffer.irred_dim)}
    {}
};

/**
 * Class containing HDF5 data sets for all objects that should be stored in the file.
 * There are two sets of members:
 *  - Pointers: when writing to a new file, new data sets are created, which require access via pointers.
 *  - Objects: when writing to an existing file, DataSet objects are read from the file, which are stored in this class.
 * There are also two constructors:
 *  - When writing file, new data sets need to be created if the file does not yet exist.
 *  - When reading file, existing data sets are opened from the file.
 */
class DataSets {
public:
    // Data sets for new file: Pointers
    H5::DataSet *lambda_p, *self_p, *irred_p;
    H5::DataSet *freq_params_p;
    H5::DataSet *bfreqsa_p, *bfreqsp_p, *bfreqst_p, *ffreqs_p, *params_p;
#if MAX_DIAG_CLASS >= 1
    H5::DataSet *K1_a_p, *K1_p_p, *K1_t_p;
#endif
#if MAX_DIAG_CLASS >= 2
    H5::DataSet *bfreqs2a_p, *bfreqs2p_p, *bfreqs2t_p, *ffreqs2a_p, *ffreqs2p_p, *ffreqs2t_p;
    H5::DataSet *K2_a_p, *K2_p_p, *K2_t_p;
#endif
#if MAX_DIAG_CLASS >= 3
    H5::DataSet *bfreqs3a_p, *bfreqs3p_p, *bfreqs3t_p, *ffreqs3a_p, *ffreqs3p_p, *ffreqs3t_p;
    H5::DataSet *K3_a_p, *K3_p_p, *K3_t_p;
#endif

    // Data sets from existing file: Objects
    H5::DataSet lambda, self, irred;
    H5::DataSet freq_params;
    H5::DataSet bfreqsa_dataset, bfreqsp_dataset, bfreqst_dataset, ffreqs_dataset;
#if MAX_DIAG_CLASS >= 1
    H5::DataSet K1_a, K1_p, K1_t;
#endif
#if MAX_DIAG_CLASS >= 2
    H5::DataSet bfreqs2a_dataset, bfreqs2p_dataset, bfreqs2t_dataset, ffreqs2a_dataset, ffreqs2p_dataset, ffreqs2t_dataset;
    H5::DataSet K2_a, K2_p, K2_t;
#endif
#if MAX_DIAG_CLASS >= 3
    H5::DataSet bfreqs3a_dataset, bfreqs3p_dataset, bfreqs3t_dataset, ffreqs3a_dataset, ffreqs3p_dataset, ffreqs3t_dataset;
    H5::DataSet K3_a, K3_p, K3_t;
#endif

    /**
     * Initialize the data sets for writing to file:
     *  - When writing to a new file, create new data sets.
     *  - When writing to an existing file, open existing data sets.
     * @param file            : File to which data sets are added / from which they are read.
     * @param file_exists     : To check if the file already existed.
     * @param dataSpaces_...  : Data spaces, needed when creating new data sets.
     * @param mtype_comp      : Data type of selfenergy and vertex (complex number).
     * @param plist_vert      : Initial value for vertex data sets.
     */
    DataSets(H5::H5File* file, bool file_exists,
             H5::DataSpace& dataSpaces_Lambda,
             H5::DataSpace& dataSpaces_selfenergy,
             H5::DataSpace& dataSpaces_irreducible,
             H5::DataSpace& dataSpaces_freq_params,
             H5::DataSpace& dataSpaces_bfreqsa,
             H5::DataSpace& dataSpaces_bfreqsp,
             H5::DataSpace& dataSpaces_bfreqst,
             H5::DataSpace& dataSpaces_ffreqs,
             H5::DataSpace& dataSpaces_params,
#if MAX_DIAG_CLASS >= 1
             H5::DataSpace& dataSpaces_K1_a,
             H5::DataSpace& dataSpaces_K1_p,
             H5::DataSpace& dataSpaces_K1_t,
#endif
#if MAX_DIAG_CLASS >= 2
             H5::DataSpace& dataSpaces_bfreqs2a,
             H5::DataSpace& dataSpaces_ffreqs2a,
             H5::DataSpace& dataSpaces_bfreqs2p,
             H5::DataSpace& dataSpaces_ffreqs2p,
             H5::DataSpace& dataSpaces_bfreqs2t,
             H5::DataSpace& dataSpaces_ffreqs2t,
             H5::DataSpace& dataSpaces_K2_a,
             H5::DataSpace& dataSpaces_K2_p,
             H5::DataSpace& dataSpaces_K2_t,
#endif
#if MAX_DIAG_CLASS >= 3
             H5::DataSpace& dataSpaces_bfreqs3a,
             H5::DataSpace& dataSpaces_ffreqs3a,
             H5::DataSpace& dataSpaces_bfreqs3p,
             H5::DataSpace& dataSpaces_ffreqs3p,
             H5::DataSpace& dataSpaces_bfreqs3t,
             H5::DataSpace& dataSpaces_ffreqs3t,
             H5::DataSpace& dataSpaces_K3_a,
             H5::DataSpace& dataSpaces_K3_p,
             H5::DataSpace& dataSpaces_K3_t,
#endif
             H5::CompType mtype_comp, H5::DSetCreatPropList plist_vert) {
        if (!file_exists) {
            lambda_p = new H5::DataSet(
                    file->createDataSet(LAMBDA_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_Lambda));
            self_p = new H5::DataSet(file->createDataSet(SELF_LIST, mtype_comp, dataSpaces_selfenergy));
            irred_p = new H5::DataSet(
                    file->createDataSet(DATASET_irred, mtype_comp, dataSpaces_irreducible, plist_vert));
            freq_params_p = new H5::DataSet(
                    file->createDataSet(FREQ_PARAMS, H5::PredType::NATIVE_DOUBLE, dataSpaces_freq_params));
            bfreqsa_p = new H5::DataSet(
                    file->createDataSet(BFREQS_LISTa, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqsa));
            bfreqsp_p = new H5::DataSet(
                    file->createDataSet(BFREQS_LISTp, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqsp));
            bfreqst_p = new H5::DataSet(
                    file->createDataSet(BFREQS_LISTt, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqst));
            ffreqs_p = new H5::DataSet(
                    file->createDataSet(FFREQS_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs));
            params_p = new H5::DataSet(
                    file->createDataSet(PARAM_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_params));
#if MAX_DIAG_CLASS >= 1
            // Create the datasets in file:
            K1_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K1_a, mtype_comp, dataSpaces_K1_a, plist_vert));
            K1_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K1_p, mtype_comp, dataSpaces_K1_p, plist_vert));
            K1_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K1_t, mtype_comp, dataSpaces_K1_t, plist_vert));
#endif
#if MAX_DIAG_CLASS >= 2
            // Create the datasets in file:
            bfreqs2a_p = new H5::DataSet(
                    file->createDataSet(BFREQS2_LISTa, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs2a));
            bfreqs2p_p = new H5::DataSet(
                    file->createDataSet(BFREQS2_LISTp, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs2p));
            bfreqs2t_p = new H5::DataSet(
                    file->createDataSet(BFREQS2_LISTt, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs2t));
            ffreqs2a_p = new H5::DataSet(
                    file->createDataSet(FFREQS2_LISTa, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs2a));
            ffreqs2p_p = new H5::DataSet(
                    file->createDataSet(FFREQS2_LISTp, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs2p));
            ffreqs2t_p = new H5::DataSet(
                    file->createDataSet(FFREQS2_LISTt, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs2t));
            K2_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_a, mtype_comp, dataSpaces_K2_a, plist_vert));
            K2_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_p, mtype_comp, dataSpaces_K2_p, plist_vert));
            K2_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_t, mtype_comp, dataSpaces_K2_t, plist_vert));
#endif
#if MAX_DIAG_CLASS >= 3
            // Create the datasets in file:
            bfreqs3a_p = new H5::DataSet(
                    file->createDataSet(BFREQS3_LISTa, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs3a));
            bfreqs3p_p = new H5::DataSet(
                    file->createDataSet(BFREQS3_LISTp, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs3p));
            bfreqs3t_p = new H5::DataSet(
                    file->createDataSet(BFREQS3_LISTt, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs3t));
            ffreqs3a_p = new H5::DataSet(
                    file->createDataSet(FFREQS3_LISTa, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs3a));
            ffreqs3p_p = new H5::DataSet(
                    file->createDataSet(FFREQS3_LISTp, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs3p));
            ffreqs3t_p = new H5::DataSet(
                    file->createDataSet(FFREQS3_LISTt, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs3t));
            K3_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K3_a, mtype_comp, dataSpaces_K3_a, plist_vert));
            K3_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K3_p, mtype_comp, dataSpaces_K3_p, plist_vert));
            K3_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K3_t, mtype_comp, dataSpaces_K3_t, plist_vert));
#endif
        }
        else {
            lambda = file->openDataSet("lambdas");
            self = file->openDataSet("selflist");
            irred = file->openDataSet("irred");
            freq_params = file->openDataSet("freq_params");
            bfreqsa_dataset = file->openDataSet("bfreqs_a");
            bfreqsp_dataset = file->openDataSet("bfreqs_p");
            bfreqst_dataset = file->openDataSet("bfreqs_t");
            ffreqs_dataset = file->openDataSet("ffreqs");
#if MAX_DIAG_CLASS >=1
            K1_a = file->openDataSet("K1_a");
            K1_p = file->openDataSet("K1_p");
            K1_t = file->openDataSet("K1_t");
#endif
#if MAX_DIAG_CLASS >=2
            bfreqs2a_dataset = file->openDataSet("bfreqs2_a");
            ffreqs2a_dataset = file->openDataSet("ffreqs2_a");
            bfreqs2p_dataset = file->openDataSet("bfreqs2_p");
            ffreqs2p_dataset = file->openDataSet("ffreqs2_p");
            bfreqs2t_dataset = file->openDataSet("bfreqs2_t");
            ffreqs2t_dataset = file->openDataSet("ffreqs2_t");
            K2_a = file->openDataSet("K2_a");
            K2_p = file->openDataSet("K2_p");
            K2_t = file->openDataSet("K2_t");
#endif
#if MAX_DIAG_CLASS >=3
            bfreqs3a_dataset = file->openDataSet("bfreqs3_a");
            ffreqs3a_dataset = file->openDataSet("ffreqs3_a");
            bfreqs3p_dataset = file->openDataSet("bfreqs3_p");
            ffreqs3p_dataset = file->openDataSet("ffreqs3_p");
            bfreqs3t_dataset = file->openDataSet("bfreqs3_t");
            ffreqs3t_dataset = file->openDataSet("ffreqs3_t");
            K3_a = file->openDataSet("K3_a");
            K3_p = file->openDataSet("K3_p");
            K3_t = file->openDataSet("K3_t");
#endif
        }
    }
    /**
     * Initialize the data sets when reading from file.
     * @param file : File from which data sets are loaded.
     */
    DataSets(H5::H5File* file) {
        lambda = file->openDataSet("lambdas");
        self = file->openDataSet("selflist");
        irred = file->openDataSet("irred");
        try {   // storing frequency gri parameters was implemented later --> old files do not have it
            freq_params = file->openDataSet("freq_params");
        }
        catch (H5::FileIException error) {
            error.printErrorStack();
        }
        bfreqsa_dataset = file->openDataSet("bfreqs_a");
        bfreqsp_dataset = file->openDataSet("bfreqs_p");
        bfreqst_dataset = file->openDataSet("bfreqs_t");
        ffreqs_dataset = file->openDataSet("ffreqs");
#if MAX_DIAG_CLASS >=1
        K1_a = file->openDataSet("K1_a");
        K1_p = file->openDataSet("K1_p");
        K1_t = file->openDataSet("K1_t");
#endif
#if MAX_DIAG_CLASS >=2
        bfreqs2a_dataset = file->openDataSet("bfreqs2_a");
        ffreqs2a_dataset = file->openDataSet("ffreqs2_a");
        bfreqs2p_dataset = file->openDataSet("bfreqs2_p");
        ffreqs2p_dataset = file->openDataSet("ffreqs2_p");
        bfreqs2t_dataset = file->openDataSet("bfreqs2_t");
        ffreqs2t_dataset = file->openDataSet("ffreqs2_t");
        K2_a = file->openDataSet("K2_a");
        K2_p = file->openDataSet("K2_p");
        K2_t = file->openDataSet("K2_t");
#endif
#if MAX_DIAG_CLASS >=3
        bfreqs3a_dataset = file->openDataSet("bfreqs3_a");
        ffreqs3a_dataset = file->openDataSet("ffreqs3_a");
        bfreqs3p_dataset = file->openDataSet("bfreqs3_p");
        ffreqs3p_dataset = file->openDataSet("ffreqs3_p");
        bfreqs3t_dataset = file->openDataSet("bfreqs3_t");
        ffreqs3t_dataset = file->openDataSet("ffreqs3_t");
        K3_a = file->openDataSet("K3_a");
        K3_p = file->openDataSet("K3_p");
        K3_t = file->openDataSet("K3_t");
#endif
    }
    /**
     * Close all data sets when writing/reading is finished.
     */
    void close(bool file_exists) {
        if (!file_exists) {
            lambda_p -> close();
            freq_params_p -> close();
            bfreqsa_p -> close();
            bfreqsp_p -> close();
            bfreqst_p -> close();
            ffreqs_p -> close();
            params_p -> close();
            irred_p  -> close();
            self_p   -> close();

#if MAX_DIAG_CLASS >= 1
            K1_a_p -> close();
            K1_p_p -> close();
            K1_t_p -> close();
#endif
#if MAX_DIAG_CLASS >= 2
            bfreqs2a_p -> close();
            ffreqs2a_p -> close();
            bfreqs2p_p -> close();
            ffreqs2p_p -> close();
            bfreqs2t_p -> close();
            ffreqs2t_p -> close();
            K2_a_p -> close();
            K2_p_p -> close();
            K2_t_p -> close();
#endif
#if MAX_DIAG_CLASS >= 3
            bfreqs3a_p -> close();
            ffreqs3a_p -> close();
            bfreqs3p_p -> close();
            ffreqs3p_p -> close();
            bfreqs3t_p -> close();
            ffreqs3t_p -> close();
            K3_a_p -> close();
            K3_p_p -> close();
            K3_t_p -> close();
#endif
        }
        else {
            self.close();
            irred.close();
            freq_params.close();
            bfreqsa_dataset.close();
            bfreqsp_dataset.close();
            bfreqst_dataset.close();
            ffreqs_dataset.close();

#if MAX_DIAG_CLASS >=1
            K1_a.close();
            K1_p.close();
            K1_t.close();
#endif
#if MAX_DIAG_CLASS >=2
            bfreqs2a_dataset.close();
            ffreqs2a_dataset.close();
            bfreqs2p_dataset.close();
            ffreqs2p_dataset.close();
            bfreqs2t_dataset.close();
            ffreqs2t_dataset.close();
            K2_a.close();
            K2_p.close();
            K2_t.close();
#endif
#if MAX_DIAG_CLASS >=3
            bfreqs3a_dataset.close();
            ffreqs3a_dataset.close();
            bfreqs3p_dataset.close();
            ffreqs3p_dataset.close();
            bfreqs3t_dataset.close();
            ffreqs3t_dataset.close();
            K3_a.close();
            K3_p.close();
            K3_t.close();
#endif
        }
    }
};


/// --- Functions for reading data from file --- ///

/**
 * Read Lambdas from an existing file
 * @param FILE_NAME   : File name
 * @return            : Lambdas in a vector of doubles
 */
rvec read_Lambdas_from_hdf(H5std_string FILE_NAME);

/**
 * Initialize the frequency grids of the result using the parameters stored in buffer.
 * @param result : Empty State object into which result is copied.
 * @param buffer : Buffer from which result is read. Should contain data read from a file.
 */
template <typename Q>
void result_set_frequency_grids(State<Q>& result, Buffer& buffer) {
    // create new frequency grids
    FrequencyGrid bfreqsa ('b', 1, Lambda_ini);
    FrequencyGrid bfreqsp ('b', 1, Lambda_ini);
    FrequencyGrid bfreqst ('b', 1, Lambda_ini);
    FrequencyGrid ffreqs ('f', 1, Lambda_ini);
    // read grid parameters from buffer
    bfreqsa.N_w = (int)buffer.freq_params[0];
    bfreqsa.w_upper = buffer.freq_params[1];
    bfreqsa.w_lower = buffer.freq_params[2];
    bfreqsa.W_scale = buffer.freq_params[3];
    ffreqs.N_w = (int)buffer.freq_params[4];
    ffreqs.w_upper = buffer.freq_params[5];
    ffreqs.w_lower = buffer.freq_params[6];
    ffreqs.W_scale = buffer.freq_params[7];
    bfreqsp.N_w = (int)buffer.freq_params[24];
    bfreqsp.w_upper = buffer.freq_params[25];
    bfreqsp.w_lower = buffer.freq_params[26];
    bfreqsp.W_scale = buffer.freq_params[27];
    bfreqst.N_w = (int)buffer.freq_params[28];
    bfreqst.w_upper = buffer.freq_params[29];
    bfreqst.w_lower = buffer.freq_params[30];
    bfreqst.W_scale = buffer.freq_params[31];
    // initialize grids
    bfreqsa.initialize_grid();
    bfreqsp.initialize_grid();
    bfreqst.initialize_grid();
    ffreqs.initialize_grid();
    // copy grids to result
    result.selfenergy.frequencies = ffreqs;
    result.vertex.avertex().K1.frequencies_K1.b = bfreqsa;
    result.vertex.pvertex().K1.frequencies_K1.b = bfreqsp;
    result.vertex.tvertex().K1.frequencies_K1.b = bfreqst;
#if MAX_DIAG_CLASS >= 2
    FrequencyGrid bfreqs2a ('b', 2, Lambda_ini);
    FrequencyGrid ffreqs2a ('f', 2, Lambda_ini);
    FrequencyGrid bfreqs2p ('b', 2, Lambda_ini);
    FrequencyGrid ffreqs2p ('f', 2, Lambda_ini);
    FrequencyGrid bfreqs2t ('b', 2, Lambda_ini);
    FrequencyGrid ffreqs2t ('f', 2, Lambda_ini);
    bfreqs2a.N_w = (int)buffer.freq_params[8];
    bfreqs2a.w_upper = buffer.freq_params[9];
    bfreqs2a.w_lower = buffer.freq_params[10];
    bfreqs2a.W_scale = buffer.freq_params[11];
    ffreqs2a.N_w = (int)buffer.freq_params[12];
    ffreqs2a.w_upper = buffer.freq_params[13];
    ffreqs2a.w_lower = buffer.freq_params[14];
    ffreqs2a.W_scale = buffer.freq_params[15];
    bfreqs2p.N_w = (int)buffer.freq_params[32];
    bfreqs2p.w_upper = buffer.freq_params[33];
    bfreqs2p.w_lower = buffer.freq_params[34];
    bfreqs2p.W_scale = buffer.freq_params[35];
    ffreqs2p.N_w = (int)buffer.freq_params[36];
    ffreqs2p.w_upper = buffer.freq_params[37];
    ffreqs2p.w_lower = buffer.freq_params[38];
    ffreqs2p.W_scale = buffer.freq_params[39];
    bfreqs2t.N_w = (int)buffer.freq_params[40];
    bfreqs2t.w_upper = buffer.freq_params[41];
    bfreqs2t.w_lower = buffer.freq_params[42];
    bfreqs2t.W_scale = buffer.freq_params[43];
    ffreqs2t.N_w = (int)buffer.freq_params[44];
    ffreqs2t.w_upper = buffer.freq_params[45];
    ffreqs2t.w_lower = buffer.freq_params[46];
    ffreqs2t.W_scale = buffer.freq_params[47];
    bfreqs2a.initialize_grid();
    ffreqs2a.initialize_grid();
    bfreqs2p.initialize_grid();
    ffreqs2p.initialize_grid();
    bfreqs2t.initialize_grid();
    ffreqs2t.initialize_grid();
    result.vertex.avertex().K2.frequencies_K2.b = bfreqs2a;
    result.vertex.pvertex().K2.frequencies_K2.b = bfreqs2a;
    result.vertex.tvertex().K2.frequencies_K2.b = bfreqs2p;
    result.vertex.avertex().K2.frequencies_K2.f = ffreqs2p;
    result.vertex.pvertex().K2.frequencies_K2.f = ffreqs2t;
    result.vertex.tvertex().K2.frequencies_K2.f = ffreqs2t;
#endif
#if MAX_DIAG_CLASS >= 3
    FrequencyGrid bfreqs3a ('b', 3, Lambda_ini);
    FrequencyGrid ffreqs3a ('f', 3, Lambda_ini);
    FrequencyGrid bfreqs3p ('b', 3, Lambda_ini);
    FrequencyGrid ffreqs3p ('f', 3, Lambda_ini);
    FrequencyGrid bfreqs3t ('b', 3, Lambda_ini);
    FrequencyGrid ffreqs3t ('f', 3, Lambda_ini);
    bfreqs3a.N_w = (int)buffer.freq_params[16];
    bfreqs3a.w_upper = buffer.freq_params[17];
    bfreqs3a.w_lower = buffer.freq_params[18];
    bfreqs3a.W_scale = buffer.freq_params[19];
    ffreqs3a.N_w = (int)buffer.freq_params[20];
    ffreqs3a.w_upper = buffer.freq_params[21];
    ffreqs3a.w_lower = buffer.freq_params[22];
    ffreqs3a.W_scale = buffer.freq_params[23];
    bfreqs3p.N_w = (int)buffer.freq_params[48];
    bfreqs3p.w_upper = buffer.freq_params[49];
    bfreqs3p.w_lower = buffer.freq_params[50];
    bfreqs3p.W_scale = buffer.freq_params[51];
    ffreqs3p.N_w = (int)buffer.freq_params[52];
    ffreqs3p.w_upper = buffer.freq_params[53];
    ffreqs3p.w_lower = buffer.freq_params[54];
    ffreqs3p.W_scale = buffer.freq_params[55];
    bfreqs3t.N_w = (int)buffer.freq_params[56];
    bfreqs3t.w_upper = buffer.freq_params[57];
    bfreqs3t.w_lower = buffer.freq_params[58];
    bfreqs3t.W_scale = buffer.freq_params[59];
    ffreqs3t.N_w = (int)buffer.freq_params[60];
    ffreqs3t.w_upper = buffer.freq_params[61];
    ffreqs3t.w_lower = buffer.freq_params[62];
    ffreqs3t.W_scale = buffer.freq_params[63];
    bfreqs3a.initialize_grid();
    ffreqs3a.initialize_grid();
    bfreqs3p.initialize_grid();
    ffreqs3p.initialize_grid();
    bfreqs3t.initialize_grid();
    ffreqs3t.initialize_grid();
    result.vertex.avertex().K3.frequencies_K3.b = bfreqs3a;
    result.vertex.pvertex().K3.frequencies_K3.b = bfreqs3a;
    result.vertex.tvertex().K3.frequencies_K3.b = bfreqs3p;
    result.vertex.avertex().K3.frequencies_K3.f = ffreqs3p;
    result.vertex.pvertex().K3.frequencies_K3.f = ffreqs3t;
    result.vertex.tvertex().K3.frequencies_K3.f = ffreqs3t;
#endif
}

/**
 * Copy results that are read from a file to a buffer into a State object.
 * @param result : Empty State object into which result is copied.
 * @param buffer : Buffer from which result is read. Should contain data read from a file.
 */
template <typename Q>
void copy_buffer_to_result(State<Q>& result, Buffer& buffer) {
    Q val; // buffer value

    result.Lambda = *buffer.lambda;

    for (int i=0; i<buffer.self_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.selfenergy[i].re, buffer.selfenergy[i].im};
#else
        val = buffer.selfenergy[i].im;
#endif
        result.selfenergy.direct_set(i, val);
    }
    for (int i=0; i<buffer.irred_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.irreducible_class[i].re, buffer.irreducible_class[i].im};
#else
        val = buffer.irreducible_class[i].re;
#endif
        result.vertex.irred().direct_set(i, val);
    }
#if MAX_DIAG_CLASS >= 1
    for (int i=0; i<buffer.K1_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K1_class_a[i].re, buffer.K1_class_a[i].im};
#else
        val = buffer.K1_class_a[i].re;
#endif
        result.vertex.avertex().K1.direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K1_class_p[i].re, buffer.K1_class_p[i].im};
#else
        val = buffer.K1_class_p[i].re;
#endif
        result.vertex.pvertex().K1.direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K1_class_t[i].re, buffer.K1_class_t[i].im};
#else
        val = buffer.K1_class_t[i].re;
#endif
        result.vertex.tvertex().K1.direct_set(i, val);
    }
#endif
#if MAX_DIAG_CLASS >= 2
    for (int i=0; i<buffer.K2_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K2_class_a[i].re, buffer.K2_class_a[i].im};
#else
        val = buffer.K2_class_a[i].re;
#endif
        result.vertex.avertex().K2.direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K2_class_p[i].re, buffer.K2_class_p[i].im};
#else
        val = buffer.K2_class_p[i].re;
#endif
        result.vertex.pvertex().K2.direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K2_class_t[i].re, buffer.K2_class_t[i].im};
#else
        val = buffer.K2_class_t[i].re;
#endif
        result.vertex.tvertex().K2.direct_set(i, val);
    }
#endif
#if MAX_DIAG_CLASS >= 3
    for (int i=0; i<buffer.K3_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K3_class_a[i].re, buffer.K3_class_a[i].im};
#else
        val = buffer.K3_class_a[i].re;
#endif
        result.vertex.avertex().K3.direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K3_class_p[i].re, buffer.K3_class_p[i].im};
#else
        val = buffer.K3_class_p[i].re;
#endif
        result.vertex.pvertex().K3.direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K3_class_t[i].re, buffer.K3_class_t[i].im};
#else
        val = buffer.K3_class_t[i].re;
#endif
        result.vertex.tvertex().K3.direct_set(i, val);
    }
#endif
}

/**
 * Read results from an existing file to a State object. Useful when resuming a computation (checkpointing).
 * @param FILE_NAME   : File name.
 * @param Lambda_it   : Lambda iteration from which to load result.
 * @param Lambda_size : Total number of Lambda iterations saved in the file.
 * @param Lambdas     : Vector containing all Lambda values for which results can be saved in file.
 * @return            : State object containing the result.
 */
State<state_datatype> read_hdf(const H5std_string FILE_NAME, size_t Lambda_it);


/// --- Functions for writing data to file --- ///

/**
 * Save data in State object to file.
 * @param FILE_NAME    : File name.
 * @param Lambda_it    : Lambda iteration at which to save the data.
 * @param Lambda_size  : Total number of Lambda iterations saved in the file.
 * @param state_in     : State to be written to the file.
 * @param Lambdas      : Vector containing all Lambda values for which results can be saved in file.
 * @param file_exists  : To check if the file already exists:
 *                        - If not, create new file.
 *                        - If yes, write data to existing file at iteration Lambda_it.
 */
template <typename  Q>
void save_to_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size,
                 const State<Q>& state_in, rvec& Lambdas, bool file_exists, const bool verbose=true) {
    //Try block to detect exceptions raised by any of the calls inside it
    try {
        // Prepare a buffer for writing data into the file
        Buffer buffer;
        buffer.initialize(state_in);    // copy data from state_in into the buffer

        // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
        H5::Exception::dontPrint();

        H5::H5File* file = 0;
        if (!file_exists) {
            // If file doesn't exist, create a new file using the default property lists.
            file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
        }
        else {
            // If file exists, open existing file. Access rights: read/write
            file = new H5::H5File(FILE_NAME, H5F_ACC_RDWR);
        }

        // Create the memory data type for storing complex numbers in file
        H5::CompType mtype_comp = def_mtype_comp();

        // Create the dimension arrays for objects in file and in buffer
        Dims dims(buffer, Lambda_size);

        // Create the data spaces for the data sets in file and for buffer objects
        H5::DataSpace dataSpaces_Lambda(1, dims.Lambda);
        H5::DataSpace dataSpaces_freq_params(RANK_freqs, dims.freq_params_dims);
        H5::DataSpace dataSpaces_bfreqsa(RANK_freqs, dims.bfreqs_dims);
        H5::DataSpace dataSpaces_bfreqsp(RANK_freqs, dims.bfreqs_dims);
        H5::DataSpace dataSpaces_bfreqst(RANK_freqs, dims.bfreqs_dims);
        H5::DataSpace dataSpaces_ffreqs(RANK_freqs, dims.ffreqs_dims);
        H5::DataSpace dataSpaces_params(1, dims.params);
        H5::DataSpace dataSpaces_selfenergy(RANK_self, dims.selfenergy);
        H5::DataSpace dataSpaces_irreducible(RANK_irreducible, dims.irreducible);

        H5::DataSpace dataSpaces_freq_params_buffer(RANK_freqs-1, dims.freq_params_buffer_dims);
        H5::DataSpace dataSpaces_bfreqsa_buffer(RANK_freqs-1, dims.bfreqs_buffer_dims);
        H5::DataSpace dataSpaces_bfreqsp_buffer(RANK_freqs-1, dims.bfreqs_buffer_dims);
        H5::DataSpace dataSpaces_bfreqst_buffer(RANK_freqs-1, dims.bfreqs_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs_buffer(RANK_freqs-1, dims.ffreqs_buffer_dims);
        H5::DataSpace dataSpaces_selfenergy_buffer(RANK_self-1, dims.selfenergy_buffer);
        H5::DataSpace dataSpaces_irreducible_buffer(RANK_irreducible-1, dims.irreducible_buffer);

#if MAX_DIAG_CLASS >= 1
        H5::DataSpace dataSpaces_K1_a(RANK_K1, dims.K1);
        H5::DataSpace dataSpaces_K1_p(RANK_K1, dims.K1);
        H5::DataSpace dataSpaces_K1_t(RANK_K1, dims.K1);

        H5::DataSpace dataSpaces_K1_a_buffer(RANK_K1-1, dims.K1_buffer);
        H5::DataSpace dataSpaces_K1_p_buffer(RANK_K1-1, dims.K1_buffer);
        H5::DataSpace dataSpaces_K1_t_buffer(RANK_K1-1, dims.K1_buffer);
#endif
#if MAX_DIAG_CLASS >= 2
        H5::DataSpace dataSpaces_bfreqs2a(RANK_freqs, dims.bfreqs2_dims);
        H5::DataSpace dataSpaces_ffreqs2a(RANK_freqs, dims.ffreqs2_dims);
        H5::DataSpace dataSpaces_bfreqs2p(RANK_freqs, dims.bfreqs2_dims);
        H5::DataSpace dataSpaces_ffreqs2p(RANK_freqs, dims.ffreqs2_dims);
        H5::DataSpace dataSpaces_bfreqs2t(RANK_freqs, dims.bfreqs2_dims);
        H5::DataSpace dataSpaces_ffreqs2t(RANK_freqs, dims.ffreqs2_dims);

        H5::DataSpace dataSpaces_bfreqs2a_buffer(RANK_freqs-1, dims.bfreqs2_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs2a_buffer(RANK_freqs-1, dims.ffreqs2_buffer_dims);
        H5::DataSpace dataSpaces_bfreqs2p_buffer(RANK_freqs-1, dims.bfreqs2_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs2p_buffer(RANK_freqs-1, dims.ffreqs2_buffer_dims);
        H5::DataSpace dataSpaces_bfreqs2t_buffer(RANK_freqs-1, dims.bfreqs2_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs2t_buffer(RANK_freqs-1, dims.ffreqs2_buffer_dims);

        H5::DataSpace dataSpaces_K2_a(RANK_K2, dims.K2);
        H5::DataSpace dataSpaces_K2_p(RANK_K2, dims.K2);
        H5::DataSpace dataSpaces_K2_t(RANK_K2, dims.K2);

        H5::DataSpace dataSpaces_K2_a_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_p_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_t_buffer(RANK_K2-1, dims.K2_buffer);
#endif
#if MAX_DIAG_CLASS >= 3
        H5::DataSpace dataSpaces_bfreqs3a(RANK_freqs, dims.bfreqs3_dims);
        H5::DataSpace dataSpaces_ffreqs3a(RANK_freqs, dims.ffreqs3_dims);
        H5::DataSpace dataSpaces_bfreqs3p(RANK_freqs, dims.bfreqs3_dims);
        H5::DataSpace dataSpaces_ffreqs3p(RANK_freqs, dims.ffreqs3_dims);
        H5::DataSpace dataSpaces_bfreqs3t(RANK_freqs, dims.bfreqs3_dims);
        H5::DataSpace dataSpaces_ffreqs3t(RANK_freqs, dims.ffreqs3_dims);

        H5::DataSpace dataSpaces_bfreqs3a_buffer(RANK_freqs-1, dims.bfreqs3_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs3a_buffer(RANK_freqs-1, dims.ffreqs3_buffer_dims);
        H5::DataSpace dataSpaces_bfreqs3p_buffer(RANK_freqs-1, dims.bfreqs3_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs3p_buffer(RANK_freqs-1, dims.ffreqs3_buffer_dims);
        H5::DataSpace dataSpaces_bfreqs3t_buffer(RANK_freqs-1, dims.bfreqs3_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs3t_buffer(RANK_freqs-1, dims.ffreqs3_buffer_dims);

        H5::DataSpace dataSpaces_K3_a(RANK_K3, dims.K3);
        H5::DataSpace dataSpaces_K3_p(RANK_K3, dims.K3);
        H5::DataSpace dataSpaces_K3_t(RANK_K3, dims.K3);

        H5::DataSpace dataSpaces_K3_a_buffer(RANK_K3-1, dims.K3_buffer);
        H5::DataSpace dataSpaces_K3_p_buffer(RANK_K3-1, dims.K3_buffer);
        H5::DataSpace dataSpaces_K3_t_buffer(RANK_K3-1, dims.K3_buffer);
#endif

        // Initial value for vertex data sets // TODO(low): remove?
        h5_comp fillvalue_vert;
        fillvalue_vert.re = 0;
        fillvalue_vert.im = 0;
        H5::DSetCreatPropList plist_vert;
        plist_vert.setFillValue(mtype_comp, &fillvalue_vert);

        // Create the data sets for all data to be saved
        DataSets dataSets(file, file_exists,
                          dataSpaces_Lambda, dataSpaces_selfenergy, dataSpaces_irreducible,
                          dataSpaces_freq_params,
                          dataSpaces_bfreqsp, dataSpaces_bfreqsp, dataSpaces_bfreqst, dataSpaces_ffreqs, dataSpaces_params,
#if MAX_DIAG_CLASS >= 1
                          dataSpaces_K1_a, dataSpaces_K1_p, dataSpaces_K1_t,
#endif
#if MAX_DIAG_CLASS >= 2
                          dataSpaces_bfreqs2a, dataSpaces_ffreqs2a,
                          dataSpaces_bfreqs2p, dataSpaces_ffreqs2p,
                          dataSpaces_bfreqs2t, dataSpaces_ffreqs2t,
                          dataSpaces_K2_a, dataSpaces_K2_p, dataSpaces_K2_t,
#endif
#if MAX_DIAG_CLASS >= 3
                          dataSpaces_bfreqs3a, dataSpaces_ffreqs3a,
                          dataSpaces_bfreqs3p, dataSpaces_ffreqs3p,
                          dataSpaces_bfreqs3t, dataSpaces_ffreqs3t,
                          dataSpaces_K3_a, dataSpaces_K3_p, dataSpaces_K3_t,
#endif
                          mtype_comp, plist_vert);

        //Select hyperslab in the file where the data should be located and after that write buffered data into file.
        hsize_t start[2];
        hsize_t stride[2];
        hsize_t count[2];
        hsize_t block[2];

        start[0] = Lambda_it;
        start[1] = 0;
        for (int i = 0; i < 2; i++) {
            stride[i] = 1;
            block[i] = 1;
        }
        count[0] = 1;

        if (!file_exists) {
            dataSets.lambda_p -> write(Lambdas.data(), H5::PredType::NATIVE_DOUBLE);
            dataSets.params_p -> write(parameter_list, H5::PredType::NATIVE_DOUBLE);
        }
        else
            // overwrite vector containing all values for lambda
            dataSets.lambda.write(Lambdas.data(), H5::PredType::NATIVE_DOUBLE);

        count[1] = N_freq_params;
        dataSpaces_freq_params.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.freq_params_p -> write(buffer.freq_params, H5::PredType::NATIVE_DOUBLE,
                                            dataSpaces_freq_params_buffer, dataSpaces_freq_params);
        else
            dataSets.freq_params.write(buffer.freq_params, H5::PredType::NATIVE_DOUBLE,
                                       dataSpaces_freq_params_buffer, dataSpaces_freq_params);


        count[1] = nBOS;
        dataSpaces_bfreqsa.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqsa_p -> write(buffer.bfreqsa_buffer, H5::PredType::NATIVE_DOUBLE,
                                       dataSpaces_bfreqsa_buffer, dataSpaces_bfreqsa);
        else
            dataSets.bfreqsa_dataset.write(buffer.bfreqsa_buffer, H5::PredType::NATIVE_DOUBLE,
                                          dataSpaces_bfreqsa_buffer, dataSpaces_bfreqsa);

        dataSpaces_bfreqsp.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqsp_p -> write(buffer.bfreqsp_buffer, H5::PredType::NATIVE_DOUBLE,
                                       dataSpaces_bfreqsp_buffer, dataSpaces_bfreqsp);
        else
            dataSets.bfreqsp_dataset.write(buffer.bfreqsp_buffer, H5::PredType::NATIVE_DOUBLE,
                                          dataSpaces_bfreqsp_buffer, dataSpaces_bfreqsp);

        dataSpaces_bfreqst.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqst_p -> write(buffer.bfreqst_buffer, H5::PredType::NATIVE_DOUBLE,
                                       dataSpaces_bfreqst_buffer, dataSpaces_bfreqst);
        else
            dataSets.bfreqst_dataset.write(buffer.bfreqst_buffer, H5::PredType::NATIVE_DOUBLE,
                                          dataSpaces_bfreqst_buffer, dataSpaces_bfreqst);


        count[1] = nFER;
        dataSpaces_ffreqs.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs_p -> write(buffer.ffreqs_buffer, H5::PredType::NATIVE_DOUBLE,
                                       dataSpaces_ffreqs_buffer, dataSpaces_ffreqs);
        else
            dataSets.ffreqs_dataset.write(buffer.ffreqs_buffer, H5::PredType::NATIVE_DOUBLE,
                                          dataSpaces_ffreqs_buffer, dataSpaces_ffreqs);

        count[1] = buffer.self_dim;
        dataSpaces_selfenergy.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.self_p -> write(buffer.selfenergy, mtype_comp,
                                     dataSpaces_selfenergy_buffer, dataSpaces_selfenergy);
        else
            dataSets.self.write(buffer.selfenergy, mtype_comp,
                                dataSpaces_selfenergy_buffer, dataSpaces_selfenergy);

        count[1] = buffer.irred_dim;
        dataSpaces_irreducible.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.irred_p -> write(buffer.irreducible_class, mtype_comp,
                                      dataSpaces_irreducible_buffer, dataSpaces_irreducible);
        else
            dataSets.irred.write(buffer.irreducible_class, mtype_comp,
                                 dataSpaces_irreducible_buffer, dataSpaces_irreducible);

#if MAX_DIAG_CLASS >= 1
        count[1]= buffer.K1_dim;
        dataSpaces_K1_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K1_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K1_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        if (!file_exists) {
            dataSets.K1_a_p -> write(buffer.K1_class_a, mtype_comp, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
            dataSets.K1_p_p -> write(buffer.K1_class_p, mtype_comp, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
            dataSets.K1_t_p -> write(buffer.K1_class_t, mtype_comp, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
        }
        else {
            dataSets.K1_a.write(buffer.K1_class_a, mtype_comp, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
            dataSets.K1_p.write(buffer.K1_class_p, mtype_comp, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
            dataSets.K1_t.write(buffer.K1_class_t, mtype_comp, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
        }
#endif

#if MAX_DIAG_CLASS >= 2
        count[1] = nBOS2;
        dataSpaces_bfreqs2a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs2a_p -> write(buffer.bfreqs2a_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs2a_buffer, dataSpaces_bfreqs2a);
        else
            dataSets.bfreqs2a_dataset.write(buffer.bfreqs2a_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs2a_buffer, dataSpaces_bfreqs2a);

        count[1] = nBOS2;
        dataSpaces_bfreqs2p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs2p_p -> write(buffer.bfreqs2p_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs2p_buffer, dataSpaces_bfreqs2p);
        else
            dataSets.bfreqs2p_dataset.write(buffer.bfreqs2p_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs2p_buffer, dataSpaces_bfreqs2p);

        count[1] = nBOS2;
        dataSpaces_bfreqs2t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs2t_p -> write(buffer.bfreqs2t_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs2t_buffer, dataSpaces_bfreqs2t);
        else
            dataSets.bfreqs2t_dataset.write(buffer.bfreqs2t_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs2t_buffer, dataSpaces_bfreqs2t);


        count[1] = nFER2;
        dataSpaces_ffreqs2a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs2a_p -> write(buffer.ffreqs2a_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs2a_buffer, dataSpaces_ffreqs2a);
        else
            dataSets.ffreqs2a_dataset.write(buffer.ffreqs2a_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs2a_buffer, dataSpaces_ffreqs2a);

        count[1] = nFER2;
        dataSpaces_ffreqs2p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs2p_p -> write(buffer.ffreqs2p_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs2p_buffer, dataSpaces_ffreqs2p);
        else
            dataSets.ffreqs2p_dataset.write(buffer.ffreqs2p_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs2p_buffer, dataSpaces_ffreqs2p);

        count[1] = nFER2;
        dataSpaces_ffreqs2t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs2t_p -> write(buffer.ffreqs2t_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs2t_buffer, dataSpaces_ffreqs2t);
        else
            dataSets.ffreqs2t_dataset.write(buffer.ffreqs2t_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs2t_buffer, dataSpaces_ffreqs2t);

        count[1]= buffer.K2_dim;
        dataSpaces_K2_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K2_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K2_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        if (!file_exists) {
            dataSets.K2_a_p -> write(buffer.K2_class_a, mtype_comp, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
            dataSets.K2_p_p -> write(buffer.K2_class_p, mtype_comp, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
            dataSets.K2_t_p -> write(buffer.K2_class_t, mtype_comp, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
        }
        else {
            dataSets.K2_a.write(buffer.K2_class_a, mtype_comp, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
            dataSets.K2_p.write(buffer.K2_class_p, mtype_comp, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
            dataSets.K2_t.write(buffer.K2_class_t, mtype_comp, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
        }
#endif

#if MAX_DIAG_CLASS >= 3
        count[1] = nBOS3;
        dataSpaces_bfreqs3a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs3a_p -> write(buffer.bfreqs3a_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs3a_buffer, dataSpaces_bfreqs3a);
        else
            dataSets.bfreqs3a_dataset.write(buffer.bfreqs3a_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs3a_buffer, dataSpaces_bfreqs3a);

        count[1] = nBOS3;
        dataSpaces_bfreqs3p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs3p_p -> write(buffer.bfreqs3p_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs3p_buffer, dataSpaces_bfreqs3p);
        else
            dataSets.bfreqs3p_dataset.write(buffer.bfreqs3p_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs3p_buffer, dataSpaces_bfreqs3p);

        count[1] = nBOS3;
        dataSpaces_bfreqs3t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs3t_p -> write(buffer.bfreqs3t_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs3t_buffer, dataSpaces_bfreqs3t);
        else
            dataSets.bfreqs3t_dataset.write(buffer.bfreqs3t_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs3t_buffer, dataSpaces_bfreqs3t);


        count[1] = nFER3;
        dataSpaces_ffreqs3a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs3a_p -> write(buffer.ffreqs3a_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs3a_buffer, dataSpaces_ffreqs3a);
        else
            dataSets.ffreqs3a_dataset.write(buffer.ffreqs3a_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs3a_buffer, dataSpaces_ffreqs3a);

        count[1] = nFER3;
        dataSpaces_ffreqs3p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs3p_p -> write(buffer.ffreqs3p_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs3p_buffer, dataSpaces_ffreqs3p);
        else
            dataSets.ffreqs3p_dataset.write(buffer.ffreqs3p_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs3p_buffer, dataSpaces_ffreqs3p);

        count[1] = nFER3;
        dataSpaces_ffreqs3t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs3t_p -> write(buffer.ffreqs3t_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs3t_buffer, dataSpaces_ffreqs3t);
        else
            dataSets.ffreqs3t_dataset.write(buffer.ffreqs3t_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs3t_buffer, dataSpaces_ffreqs3t);

        count[1]= buffer.K3_dim;
        dataSpaces_K3_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K3_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K3_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        if (!file_exists) {
            dataSets.K3_a_p -> write(buffer.K3_class_a, mtype_comp, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
            dataSets.K3_p_p -> write(buffer.K3_class_p, mtype_comp, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
            dataSets.K3_t_p -> write(buffer.K3_class_t, mtype_comp, dataSpaces_K3_t_buffer, dataSpaces_K3_t);
        }
        else {
            dataSets.K3_a.write(buffer.K3_class_a, mtype_comp, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
            dataSets.K3_p.write(buffer.K3_class_p, mtype_comp, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
            dataSets.K3_t.write(buffer.K3_class_t, mtype_comp, dataSpaces_K3_t_buffer, dataSpaces_K3_t);
        }
#endif
        if (verbose) {
            print("Successfully saved in hdf5 file: ", FILE_NAME);
            if (file_exists) print_add(" in Lambda-layer ", Lambda_it, false);
            print_add("", true);
        }

        // Terminate
        dataSets.close(file_exists);
        file -> close();
        delete file;

    }  // end of try block

    // catch failure caused by the H5File operations
    catch (H5::FileIException error) {
        error.printErrorStack();
        return;
    }

    // catch failure caused by the DataSet operations
    catch (H5::DataSetIException error) {
        error.printErrorStack();
        return;
    }

    // catch failure caused by the DataSpace operations
    catch (H5::DataSpaceIException error) {
        error.printErrorStack();
        return;
    }
}

/**
 * Write the inital state to an HDF5 file.
 * @param FILE_NAME   : File name. Creates a new file if it does not exist.
 *                      If a file with this name already exists, overwrite the existing file.
 * @param Lambda_i    : Inital Lambda value.
 * @param Lambda_size : Total number of Lambda iterations to be saved in the file.
 * @param state_in    : State to be written to the file.
 */
template <typename Q>
void write_hdf(const H5std_string FILE_NAME, double Lambda_i, long Lambda_size, const State<Q>& state_in) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
    int Lambda_it = 0;  // store data as 0th Lambda iteration

    // List with Lambda values where only the first one is non-zero
    rvec Lambdas (Lambda_size);
    Lambdas[0] = Lambda_i;
    for (int i = 1; i < Lambda_size; i++) {
        Lambdas[i] = 0;
    }

    // write data to file
    save_to_hdf(FILE_NAME, Lambda_it, Lambda_size, state_in, Lambdas, false);
    }
}

/**
 * Add the state of a new iteration to an existing HDF5 file.
 * @param FILE_NAME   : File name.
 * @param Lambda_it   : Lambda iteration at which to save the data.
 * @param Lambda_size : Total number of Lambda iterations saved in the file.
 * @param state_in    : State to be written to the file.
 * @param Lambdas     : Vector containing all Lambda values for which results can be saved in file.
 */
template <typename Q>
void add_hdf(const H5std_string FILE_NAME, int Lambda_it, const State<Q>& state_in, rvec& Lambdas, const bool verbose=true) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        long Lambda_size = read_Lambdas_from_hdf(FILE_NAME).size();
        // write data to file if Lambda iteration number is in allowed range, otherwise print error message
        if (Lambda_it < Lambda_size) {
            save_to_hdf(FILE_NAME, Lambda_it, Lambda_size, state_in, Lambdas, true, verbose);
        } else {
            print("Cannot write to file ", FILE_NAME, " since Lambda layer", Lambda_it, " is out of range.", true);
        }
    }
}
/// Overload of above function that only updates the Lambda at iteration Lambda_it
template <typename Q>
void add_hdf(const H5std_string FILE_NAME, const double Lambda_now, const int Lambda_it, const State<Q>& state_in, const bool verbose=true) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        rvec Lambdas = read_Lambdas_from_hdf(FILE_NAME);
        Lambdas[Lambda_it] = Lambda_now; // update Lambda
        add_hdf<Q>(FILE_NAME, Lambda_it, state_in, Lambdas, verbose);
    }
}


/** overload of add_hdf for non-States, does not do anything */
template <typename Q>
void add_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size,
             Q& state_in, rvec& Lambdas) {}

/// --- Test function --- ///

void test_hdf5(H5std_string FILE_NAME, int i, State<state_datatype>& state);

bool test_read_write_data_hdf();
bool test_read_write_state_hdf();

#endif //KELDYSH_MFRG_HDF5_ROUTINES_HPP