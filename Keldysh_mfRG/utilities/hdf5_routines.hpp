/**
 * Functions to write/read a State object to/from an HDF5 file.
 */

#ifndef KELDYSH_MFRG_HDF5_ROUTINES_HPP
#define KELDYSH_MFRG_HDF5_ROUTINES_HPP

#include <stdexcept>
#include <cmath>
#include <vector>
#include "../parameters/master_parameters.hpp"         // system parameters (necessary for vector lengths etc.)
#include "util.hpp"               // printing text
#include "../data_structures.hpp"    // comp data type, std::real/complex vector class
#include "../grids/frequency_grid.hpp"     // store frequency grid parameters
#include "H5Cpp.h"              // HDF5 functions
#include "../correlation_functions/state.hpp"
#include "../multidimensional/multiarray.hpp"
#include "../symmetries/Keldysh_symmetries.hpp"


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

const H5std_string	DATASET_K2b_a("K2b_a");
const H5std_string	DATASET_K2b_p("K2b_p");
const H5std_string	DATASET_K2b_t("K2b_t");

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
const H5std_string  BFREQS2b_LISTa("bfreqs2b_a");
const H5std_string  BFREQS2b_LISTp("bfreqs2b_p");
const H5std_string  BFREQS2b_LISTt("bfreqs2b_t");
const H5std_string  BFREQS3_LISTa("bfreqs3_a");
const H5std_string  BFREQS3_LISTp("bfreqs3_p");
const H5std_string  BFREQS3_LISTt("bfreqs3_t");
const H5std_string  FFREQS_LIST ("ffreqs");
const H5std_string  FFREQS2_LISTa("ffreqs2_a");
const H5std_string  FFREQS2_LISTp("ffreqs2_p");
const H5std_string  FFREQS2_LISTt("ffreqs2_t");
const H5std_string  FFREQS2b_LISTa("ffreqs2b_a");
const H5std_string  FFREQS2b_LISTp("ffreqs2b_p");
const H5std_string  FFREQS2b_LISTt("ffreqs2b_t");
const H5std_string  FFREQS3_LISTa("ffreqs3_a");
const H5std_string  FFREQS3_LISTp("ffreqs3_p");
const H5std_string  FFREQS3_LISTt("ffreqs3_t");
const H5std_string  FFREQS3_LISTa2("ffreqs3_a2");
const H5std_string  FFREQS3_LISTp2("ffreqs3_p2");
const H5std_string  FFREQS3_LISTt2("ffreqs3_t2");
const H5std_string  FREQ_PARAMS("freq_params");
const H5std_string  PARAM_LIST("parameters");
const H5std_string  IS_CONVERGED("is_converged");
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

//const H5::CompType mtype_comp = def_mtype_comp();
H5::DSetCreatPropList def_proplist_comp();
//const H5::DSetCreatPropList plist_vert_comp = def_proplist_comp();

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
        H5::DSetCreatPropList plist_vert = def_proplist_comp();
        H5::CompType mtype_comp = def_mtype_comp();
        H5::DataSet mydataset = group.createDataSet(dataset_name, mtype_comp, file_space, plist_vert);
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
        else if constexpr(std::is_same_v<typename container::value_type, comp>) {
            H5::CompType mtype_comp = def_mtype_comp();
            dataset.write(data.data(), mtype_comp);
        }
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
        else if constexpr(std::is_same_v<typename container::value_type, comp>) {
            H5::CompType mtype_comp = def_mtype_comp();
            dataset.write(data.data(), mtype_comp, mem_space, file_space);
        }
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
        for (hsize_t i = 0; i < RANK+1; i++) {
            stride[i] = 1;
            block[i] = 1;
        }
        count[0] = 1;           // dimension of Lambda layer
        start[0] = Lambda_it;
        dims_file[0] = numberLambda_layers;

        for (hsize_t i = 1; i < RANK+1; i++) {
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
        // choose the length of this array a bit bigger, so we can also load vectors of shape (Lambda_size, 1)
        hsize_t dims_file[depth+1];

#ifndef NDEBUG
        hsize_t rank =
#endif
        file_space.getSimpleExtentDims(dims_file, NULL);
        assert(rank == depth or (dims_file[0] and depth==1));  // Important assertion! Makes sure that the desired multiarray depth matches the data in the file


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
        else if constexpr(std::is_same_v<Q,comp>) {
            H5::CompType mtype_comp = def_mtype_comp();
            dataset.read(result.data(), mtype_comp, dataSpace_buffer, file_space);   /// Funktioniert das auch für comp?
        }
        else if constexpr(std::is_same_v<Q,int>) dataset.read(result.data(), H5::PredType::NATIVE_INT, dataSpace_buffer, file_space);
        else if constexpr(std::is_same_v<Q,char>) dataset.read(result.data(), H5::PredType::NATIVE_CHAR, dataSpace_buffer, file_space);


        for (size_t i = 0; i < depth; i++) dims_result[i] = dims_file[i];
        return result;
    }


    template<typename Q, std::size_t depth, typename H5object>
    std::vector<Q> read_from_hdf_LambdaLayer_impl(const H5object& group, const H5std_string& dataset_name, std::array<std::size_t,depth>& dims_result, const int Lambda_it) {

        // open file_space to read from
        H5::DataSet dataset = hdf5_impl::open_Dataset(group, dataset_name);
        H5::DataSpace file_space = dataset.getSpace();

        // get the size of the file_space (including all Lambda layers)
        hsize_t dims_file[depth+1];       // size of the full dataspace (including all lambda layers)
#ifndef NDEBUG
        hsize_t rank =
#endif
        file_space.getSimpleExtentDims(dims_file, NULL);
        assert(rank == depth+1);  // Important assertion! Makes sure that the desired multiarray depth matches the data in the file
        assert(Lambda_it < dims_file[0]); // if this fails, Lambda_it exceeded the available Lambda layers

        // create arrays with information on dimensions
        const hsize_t RANK = depth;
        hsize_t start[RANK+1];      // starting location for the hyperslab
        hsize_t stride[RANK+1];     // number of elements to separate each element or block
        hsize_t count[RANK+1];      // number of elements or blocks to select along each dimension
        hsize_t block[RANK+1];      // size of the block selected from the dataspace
        hsize_t dims_mem[RANK];       // size of the full dataspace (including all lambda layers)
        for (size_t i = 0; i < RANK+1; i++) {
            stride[i] = 1;
            block[i] = 1;
        }
        count[0] = 1;           // dimension of Lambda layer
        start[0] = Lambda_it;
        for (size_t i = 1; i < RANK+1; i++) {
            count[i] = dims_file[i];
            start[i] = 0;
            dims_result[i-1] = dims_file[i];
            dims_mem[i-1] = dims_file[i];
        }
        hsize_t dims_mem_flat = 1;
        for (hsize_t i = 0; i < depth; i++) {
            dims_mem_flat *= dims_file[i+1];
        }


        // select hyperslab to read from
        file_space.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        H5::DataSet mydataset;

        H5::DataSpace dataSpace_buffer(depth, dims_mem);

        // create a vector the same size as the dataset
        std::vector<Q> result(dims_mem_flat);


        if constexpr(std::is_same_v<Q,double>) dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE, dataSpace_buffer, file_space);
        else if constexpr(std::is_same_v<Q,comp>) {
            H5::CompType mtype_comp = def_mtype_comp();
            dataset.read(result.data(), mtype_comp, dataSpace_buffer, file_space);   /// Funktioniert das auch für comp?
        }
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
                                std::is_same_v<Q, bool> ||
                                std::is_same_v<Q, char>,
                bool> = true>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const Q& data, const bool data_set_exists) {
    H5::DataSpace file_space(H5S_SCALAR);

    if constexpr(std::is_same_v<Q, double>) {
        H5::Attribute attr;
        if (!data_set_exists) attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_DOUBLE, file_space );
        else attr = group.openAttribute( dataset_name );
        attr.write( H5::PredType::NATIVE_DOUBLE, &data );
    }
    else if constexpr(std::is_same_v<Q, comp>) {
        H5::CompType mtype_comp = def_mtype_comp();
        H5::Attribute attr;
        if (!data_set_exists) attr = group.createAttribute( dataset_name, mtype_comp, file_space );
        else attr = group.openAttribute( dataset_name );
        attr.write( mtype_comp, &data );
    }
    else if constexpr(std::is_same_v<Q, int>) {
        H5::Attribute attr;
        if (!data_set_exists) attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_INT, file_space );
        else attr = group.openAttribute( dataset_name );
        attr.write( H5::PredType::NATIVE_INT, &data );
    }
    else if constexpr(std::is_same_v<Q, char>) {
        H5::Attribute attr;
        if (!data_set_exists) attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_CHAR, file_space );
        else attr = group.openAttribute( dataset_name );
        attr.write( H5::PredType::NATIVE_CHAR, &data );
    }
    else if constexpr(std::is_same_v<Q, bool>) {
        H5::Attribute attr;
        if (!data_set_exists) attr = group.createAttribute( dataset_name, H5::PredType::NATIVE_HBOOL, file_space );
        else attr = group.openAttribute( dataset_name );
        attr.write( H5::PredType::NATIVE_HBOOL, &data );
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

/// Write Eigen::Matrix to HDF group/file
template<typename Q, typename H5object, int nrows, int ncols>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const Eigen::Matrix<Q,nrows, ncols>& data, const bool data_set_exists) {
    hsize_t dims[2] = {(hsize_t)data.cols(), (hsize_t)data.rows()}; // standard Eigen::Matrix is stored in column-major format
    H5::DataSpace file_space(2, dims);

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





/// Read scalar from HDF group/file as Attribute
template<typename Q, typename H5object,
        std::enable_if_t<
                std::is_same_v<Q, double> ||
                std::is_same_v<Q, comp> ||
                std::is_same_v<Q, int> ||
                std::is_same_v<Q, bool> ||
                std::is_same_v<Q, char>,
                bool> = true>
void read_from_hdf(H5object& group, const H5std_string& dataset_name, Q& result) {
    H5::DataSpace file_space(H5S_SCALAR);

    H5::Attribute attr = group.openAttribute( dataset_name );
    H5::DataType type = attr.getDataType();
    attr.read(type,&result);

    file_space.close();
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
    template<typename gridType>
    void write_freqparams_to_hdf_LambdaLayer(H5::Group& group, const gridType& freqgrid, int Lambda_it, int numberLambdaLayers, bool file_exists, bool verbose);
    template<typename gridType>
    void init_freqgrid_from_hdf_LambdaLayer(H5::Group& group, gridType& freqgrid, int Lambda_it, double Lambda);

    template<typename Q, bool diff>
    void write_state_to_hdf_LambdaLayer(const H5std_string& filename, const State<Q, diff>& state, const int Lambda_it, const int numberLambdaLayers, const std::string write_mode, const bool is_converged=false, const bool verbose=true) {
        H5::H5File file_out;
        //H5::Exception::dontPrint();
        const bool keep_existing_file = (write_mode == "rw");
        bool is_dataset_existent = (write_mode == "rw");
        assert(write_mode == "w" or write_mode == "rw");

        if (!keep_existing_file) {
            // Create a new file using the default property lists. (No matter of file exists)
            file_out = H5::H5File(filename, H5F_ACC_TRUNC);
        }
        else {
            // If file exists, open existing file. Access rights: read/write
            try {
                file_out = H5::H5File(filename, H5F_ACC_RDWR);
            } catch(const H5::FileIException&) {
                // If file does not exist, open new file
                file_out = H5::H5File(filename, H5F_ACC_TRUNC);
                is_dataset_existent = false;
            }
        }

        /// Save state data
        write_to_hdf_LambdaLayer<double>(file_out, LAMBDA_LIST, std::vector<double>({state.Lambda}), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, SELF_LIST, state.selfenergy.Sigma.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, HARTREE, std::vector<Q>({state.selfenergy.asymp_val_R}), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_irred, state.vertex.irred().get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K1_a, state.vertex.avertex().K1.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K1_p, state.vertex.pvertex().K1.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K1_t, state.vertex.tvertex().K1.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS_LIST, state.selfenergy.Sigma.frequencies.  primary_grid.get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS_LISTa, state.vertex.avertex().K1.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS_LISTp, state.vertex.pvertex().K1.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS_LISTt, state.vertex.tvertex().K1.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
#if MAX_DIAG_CLASS>1
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2_a, state.vertex.avertex().K2.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2_p, state.vertex.pvertex().K2.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2_t, state.vertex.tvertex().K2.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2_LISTa, state.vertex.avertex().K2.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2_LISTp, state.vertex.pvertex().K2.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2_LISTt, state.vertex.tvertex().K2.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2_LISTa, state.vertex.avertex().K2.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2_LISTp, state.vertex.pvertex().K2.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2_LISTt, state.vertex.tvertex().K2.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
#if DEBUG_SYMMETRIES
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2b_a, state.vertex.avertex().K2b.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2b_p, state.vertex.pvertex().K2b.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K2b_t, state.vertex.tvertex().K2b.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2b_LISTa, state.vertex.avertex().K2b.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2b_LISTp, state.vertex.pvertex().K2b.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS2b_LISTt, state.vertex.tvertex().K2b.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2b_LISTa, state.vertex.avertex().K2b.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2b_LISTp, state.vertex.pvertex().K2b.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS2b_LISTt, state.vertex.tvertex().K2b.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
#endif
#endif
#if MAX_DIAG_CLASS>2
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K3_a, state.vertex.avertex().K3.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K3_p, state.vertex.pvertex().K3.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<Q>(file_out, DATASET_K3_t, state.vertex.tvertex().K3.get_vec(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS3_LISTa, state.vertex.avertex().K3.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS3_LISTp, state.vertex.pvertex().K3.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, BFREQS3_LISTt, state.vertex.tvertex().K3.frequencies.get_freqGrid_b().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTa, state.vertex.avertex().K3.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTp, state.vertex.pvertex().K3.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTt, state.vertex.tvertex().K3.frequencies.get_freqGrid_f().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTa2,state.vertex.avertex().K3.frequencies.get_freqGrid_3().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTp2,state.vertex.pvertex().K3.frequencies.get_freqGrid_3().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
        write_to_hdf_LambdaLayer<double>(file_out, FFREQS3_LISTt2,state.vertex.tvertex().K3.frequencies.get_freqGrid_3().get_all_frequencies(), Lambda_it, numberLambdaLayers, is_dataset_existent);
#endif
        if(!is_dataset_existent) {
            /// Write used parameters for documentation purpose
            H5::Group group_params(file_out.createGroup(PARAM_LIST));
            write_to_hdf(group_params, "REG", REG, is_dataset_existent);
            write_to_hdf(group_params, "Gamma", glb_Gamma, is_dataset_existent);
            write_to_hdf(group_params, "MAX_DIAG_CLASS", MAX_DIAG_CLASS, is_dataset_existent);
            write_to_hdf(group_params, "N_LOOPS", N_LOOPS, is_dataset_existent);
            write_to_hdf(group_params, "T", glb_T, is_dataset_existent);
            write_to_hdf(group_params, "mu", glb_mu, is_dataset_existent);
            write_to_hdf(group_params, "U", glb_U, is_dataset_existent);
            write_to_hdf(group_params, "epsilon", glb_epsilon, is_dataset_existent);
            write_to_hdf(group_params, "V", glb_V, is_dataset_existent);
            write_to_hdf(group_params, "ODEsolver", ODEsolver, is_dataset_existent);
            write_to_hdf(group_params, "GRID", GRID, is_dataset_existent);
            write_to_hdf(group_params, IS_CONVERGED, is_converged, is_dataset_existent);
        }
        H5::Group group_params(file_out.openGroup(PARAM_LIST));
        write_to_hdf(group_params, "last_Lambda_it", Lambda_it, is_dataset_existent);
        if (is_converged) {
            write_to_hdf(group_params, IS_CONVERGED, is_converged, is_dataset_existent);
        }

        H5::Group group_freqparams;
        H5::Group group_freqparams_ffreqs;
        H5::Group group_freqparams_bfreqsa;
        H5::Group group_freqparams_bfreqsp;
        H5::Group group_freqparams_bfreqst;
        H5::Group group_freqparams_bfreqs2a;
        H5::Group group_freqparams_bfreqs2p;
        H5::Group group_freqparams_bfreqs2t;
        H5::Group group_freqparams_bfreqs2ba;
        H5::Group group_freqparams_bfreqs2bp;
        H5::Group group_freqparams_bfreqs2bt;
        H5::Group group_freqparams_bfreqs3a;
        H5::Group group_freqparams_bfreqs3p;
        H5::Group group_freqparams_bfreqs3t;
        H5::Group group_freqparams_ffreqs2a;
        H5::Group group_freqparams_ffreqs2p;
        H5::Group group_freqparams_ffreqs2t;
        H5::Group group_freqparams_ffreqs2ba;
        H5::Group group_freqparams_ffreqs2bp;
        H5::Group group_freqparams_ffreqs2bt;
        H5::Group group_freqparams_ffreqs3a;
        H5::Group group_freqparams_ffreqs3p;
        H5::Group group_freqparams_ffreqs3t;
        H5::Group group_freqparams_ffreqs3a2;
        H5::Group group_freqparams_ffreqs3p2;
        H5::Group group_freqparams_ffreqs3t2;
    if (!keep_existing_file) {
        group_freqparams = file_out.createGroup(FREQ_PARAMS);
        group_freqparams_ffreqs = group_freqparams.createGroup(FFREQS_LIST );
        group_freqparams_bfreqsa = group_freqparams.createGroup(BFREQS_LISTa );
        group_freqparams_bfreqsp = group_freqparams.createGroup(BFREQS_LISTp );
        group_freqparams_bfreqst = group_freqparams.createGroup(BFREQS_LISTt );
        group_freqparams_bfreqs2a = group_freqparams.createGroup(BFREQS2_LISTa);
        group_freqparams_bfreqs2p = group_freqparams.createGroup(BFREQS2_LISTp);
        group_freqparams_bfreqs2t = group_freqparams.createGroup(BFREQS2_LISTt);
        group_freqparams_bfreqs2ba = group_freqparams.createGroup(BFREQS2b_LISTa);
        group_freqparams_bfreqs2bp = group_freqparams.createGroup(BFREQS2b_LISTp);
        group_freqparams_bfreqs2bt = group_freqparams.createGroup(BFREQS2b_LISTt);
        group_freqparams_bfreqs3a = group_freqparams.createGroup(BFREQS3_LISTa);
        group_freqparams_bfreqs3p = group_freqparams.createGroup(BFREQS3_LISTp);
        group_freqparams_bfreqs3t = group_freqparams.createGroup(BFREQS3_LISTt);
        group_freqparams_ffreqs2a = group_freqparams.createGroup(FFREQS2_LISTa);
        group_freqparams_ffreqs2p = group_freqparams.createGroup(FFREQS2_LISTp);
        group_freqparams_ffreqs2t = group_freqparams.createGroup(FFREQS2_LISTt);
        group_freqparams_ffreqs2ba = group_freqparams.createGroup(FFREQS2b_LISTa);
        group_freqparams_ffreqs2bp = group_freqparams.createGroup(FFREQS2b_LISTp);
        group_freqparams_ffreqs2bt = group_freqparams.createGroup(FFREQS2b_LISTt);
        group_freqparams_ffreqs3a = group_freqparams.createGroup(FFREQS3_LISTa);
        group_freqparams_ffreqs3p = group_freqparams.createGroup(FFREQS3_LISTp);
        group_freqparams_ffreqs3t = group_freqparams.createGroup(FFREQS3_LISTt);
        group_freqparams_ffreqs3a2= group_freqparams.createGroup(FFREQS3_LISTa2);
        group_freqparams_ffreqs3p2= group_freqparams.createGroup(FFREQS3_LISTp2);
        group_freqparams_ffreqs3t2= group_freqparams.createGroup(FFREQS3_LISTt2);
    }
    else {
        group_freqparams = file_out.openGroup(FREQ_PARAMS);
        group_freqparams_ffreqs = group_freqparams.openGroup(FFREQS_LIST );
        group_freqparams_bfreqsa = group_freqparams.openGroup(BFREQS_LISTa );
        group_freqparams_bfreqsp = group_freqparams.openGroup(BFREQS_LISTp );
        group_freqparams_bfreqst = group_freqparams.openGroup(BFREQS_LISTt );
        group_freqparams_bfreqs2a = group_freqparams.openGroup(BFREQS2_LISTa);
        group_freqparams_bfreqs2p = group_freqparams.openGroup(BFREQS2_LISTp);
        group_freqparams_bfreqs2t = group_freqparams.openGroup(BFREQS2_LISTt);
        group_freqparams_bfreqs2ba = group_freqparams.openGroup(BFREQS2b_LISTa);
        group_freqparams_bfreqs2bp = group_freqparams.openGroup(BFREQS2b_LISTp);
        group_freqparams_bfreqs2bt = group_freqparams.openGroup(BFREQS2b_LISTt);
        group_freqparams_bfreqs3a = group_freqparams.openGroup(BFREQS3_LISTa);
        group_freqparams_bfreqs3p = group_freqparams.openGroup(BFREQS3_LISTp);
        group_freqparams_bfreqs3t = group_freqparams.openGroup(BFREQS3_LISTt);
        group_freqparams_ffreqs2a = group_freqparams.openGroup(FFREQS2_LISTa);
        group_freqparams_ffreqs2p = group_freqparams.openGroup(FFREQS2_LISTp);
        group_freqparams_ffreqs2t = group_freqparams.openGroup(FFREQS2_LISTt);
        group_freqparams_ffreqs2ba = group_freqparams.openGroup(FFREQS2b_LISTa);
        group_freqparams_ffreqs2bp = group_freqparams.openGroup(FFREQS2b_LISTp);
        group_freqparams_ffreqs2bt = group_freqparams.openGroup(FFREQS2b_LISTt);
        group_freqparams_ffreqs3a = group_freqparams.openGroup(FFREQS3_LISTa);
        group_freqparams_ffreqs3p = group_freqparams.openGroup(FFREQS3_LISTp);
        group_freqparams_ffreqs3t = group_freqparams.openGroup(FFREQS3_LISTt);
        group_freqparams_ffreqs3a2= group_freqparams.openGroup(FFREQS3_LISTa2);
        group_freqparams_ffreqs3p2= group_freqparams.openGroup(FFREQS3_LISTp2);
        group_freqparams_ffreqs3t2= group_freqparams.openGroup(FFREQS3_LISTt2);

    }
        /// Write frequency parameters
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs , state.selfenergy.Sigma.frequencies.get_freqGrid_b()      , Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqsa, state.vertex.avertex().K1.frequencies.get_freqGrid_b()  , Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqsp, state.vertex.pvertex().K1.frequencies.get_freqGrid_b()  , Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqst, state.vertex.tvertex().K1.frequencies.get_freqGrid_b()  , Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
#if MAX_DIAG_CLASS>1
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2a,state.vertex.avertex().K2.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2p,state.vertex.pvertex().K2.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2t,state.vertex.tvertex().K2.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2a,state.vertex.avertex().K2.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2p,state.vertex.pvertex().K2.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2t,state.vertex.tvertex().K2.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
#if DEBUG_SYMMETRIES
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2ba,state.vertex.avertex().K2b.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2bp,state.vertex.pvertex().K2b.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs2bt,state.vertex.tvertex().K2b.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2ba,state.vertex.avertex().K2b.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2bp,state.vertex.pvertex().K2b.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs2bt,state.vertex.tvertex().K2b.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
#endif
#endif
#if MAX_DIAG_CLASS>2
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs3a, state.vertex.avertex().K3.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs3p, state.vertex.pvertex().K3.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_bfreqs3t, state.vertex.tvertex().K3.frequencies.get_freqGrid_b(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3a, state.vertex.avertex().K3.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3p, state.vertex.pvertex().K3.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3t, state.vertex.tvertex().K3.frequencies.get_freqGrid_f(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3a2,state.vertex.avertex().K3.frequencies.get_freqGrid_3(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3p2,state.vertex.pvertex().K3.frequencies.get_freqGrid_3(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
        write_freqparams_to_hdf_LambdaLayer(group_freqparams_ffreqs3t2,state.vertex.tvertex().K3.frequencies.get_freqGrid_3(), Lambda_it, numberLambdaLayers, is_dataset_existent, verbose);
#endif

        file_out.close();
    }

}

/// Create file with fixed number of Lambda layers and save state to first Lambda layer
template <typename Q, bool diff>
void write_state_to_hdf(const H5std_string FILE_NAME, double Lambda_i, const int Lambda_size, const State<Q,diff>& state_in, const bool verbose=true) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        hdf5_impl::write_state_to_hdf_LambdaLayer(FILE_NAME, state_in, 0, Lambda_size, "w", false);
        if (verbose) {
            utils::print("Successfully saved in hdf5 file: ", FILE_NAME);
            utils::print_add(" in Lambda-layer ", 0, false);
            utils::print_add("", true);
        }
    }
}

/// Open file and save state to a specified Lambda layer
template <typename Q, bool diff>
void add_state_to_hdf(const H5std_string FILE_NAME, int Lambda_it, const State<Q,diff>& state_in, const bool is_converged=false, const bool verbose=true) {
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
            hdf5_impl::write_state_to_hdf_LambdaLayer(FILE_NAME, state_in, Lambda_it, Lambda_size, "rw", is_converged);
            if (verbose) {
                utils::print("Successfully saved in hdf5 file: ", FILE_NAME);
                utils::print_add(" in Lambda-layer ", Lambda_it, false);
                utils::print_add("", true);
            }
        } else {
            utils::print("\t\t  ERROR: Cannot write to file ", FILE_NAME, " since Lambda layer", Lambda_it,
                  " is out of range.", "\n");
        }
    }
}

/// Read state from specified Lambda layer of hdf file
State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it) ;
bool check_convergence_hdf(const H5std_string& filename, int& Lambda_it) ;


/// --- Functions for reading data from file --- ///

/**
 * Read Lambdas from an existing file
 * @param FILE_NAME   : File name
 * @return            : Lambdas in a vector of doubles
 */
rvec read_Lambdas_from_hdf(H5std_string FILE_NAME);



/// --- Test function --- ///
bool test_read_write_data_hdf(bool verbose);
bool test_read_write_state_hdf(bool verbose);

#endif //KELDYSH_MFRG_HDF5_ROUTINES_HPP