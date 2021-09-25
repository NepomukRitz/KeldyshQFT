/**
 * Functions to write/read a State object to/from an HDF5 file.
 */

#ifndef KELDYSH_MFRG_HDF5_ROUTINES_H
#define KELDYSH_MFRG_HDF5_ROUTINES_H

#include <cmath>
#include "util.h"               // printing text
#include "../parameters/master_parameters.h"         // system parameters (necessary for vector lengths etc.)
#include "../data_structures.h"    // comp data type, std::real/complex vector class
#include "../grids/frequency_grid.h"     // store frequency grid parameters
#include "H5Cpp.h"              // HDF5 functions

#ifdef MPI_FLAG
#include "mpi_setup.h"          // mpi routines: when using mpi, only the process with ID 0 writes into file
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

const int N_freq_params = 24; // number of parameters for frequency grid

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
const H5std_string	LAMBDA_LIST("lambdas");
const H5std_string  BFREQS_LIST ("bfreqs");
const H5std_string  BFREQS2_LIST("bfreqs2");
const H5std_string  BFREQS3_LIST("bfreqs3");
const H5std_string  FFREQS_LIST ("ffreqs");
const H5std_string  FFREQS2_LIST("ffreqs2");
const H5std_string  FFREQS3_LIST("ffreqs3");
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
H5::CompType def_mtype_comp() {
    H5::CompType mtype_comp(sizeof(h5_comp));
    mtype_comp.insertMember(RE, HOFFSET(h5_comp, re), H5::PredType::NATIVE_DOUBLE);
    mtype_comp.insertMember(IM, HOFFSET(h5_comp, im), H5::PredType::NATIVE_DOUBLE);
    return mtype_comp;
}

/// --- Helper classes: buffer, dimension arrays, data sets --- ///

/**
 * Class containing buffer lengths and arrays to buffer data for selfenergy and irreducible vertex
 * as well as K1, K2, K3 classes.
 * Constructor creates new arrays, destructor deletes them.
 */
class Buffer {
public:
    double * lambda;
    double * freq_params;
    double * bfreqs_buffer;
    double * ffreqs_buffer;
#ifdef KELDYSH_FORMALISM
    const int self_dim = 2 * nSE;                                     // length of self-energy buffer
    const int irred_dim = 16 * n_in;                                  // length of irreducible vertex buffer
#else
    const int self_dim = nSE;                                         // length of self-energy buffer
    const int irred_dim = n_in;                                       // length of irreducible vertex buffer
#endif
    h5_comp * selfenergy;
    h5_comp * irreducible_class;

#if MAX_DIAG_CLASS >= 1
    const int K1_dim = nK_K1 * nw1_t * n_in;                         // length of K1 buffer
    h5_comp * K1_class_a;
    h5_comp * K1_class_p;
    h5_comp * K1_class_t;
#endif
#if MAX_DIAG_CLASS >= 2
    double * bfreqs2_buffer;
    double * ffreqs2_buffer;
    const int K2_dim = nK_K2 * nw2_t * nv2_t * n_in;               // length of K2 buffer
    h5_comp * K2_class_a;
    h5_comp * K2_class_p;
    h5_comp * K2_class_t;
#endif
#if MAX_DIAG_CLASS >= 3
    double * bfreqs3_buffer;
    double * ffreqs3_buffer;
    const int K3_dim = nK_K3 * nw3_t * nv3_t * nv3_t * n_in;    // length of K3 buffer
    h5_comp * K3_class_a;
    h5_comp * K3_class_p;
    h5_comp * K3_class_t;
#endif

    Buffer() {
        freq_params = new double[N_freq_params];
        bfreqs_buffer = new double[nBOS];                        // create buffer for bosonic frequencies
        ffreqs_buffer = new double[nFER];                        // create buffer for fermionic frequencies
        selfenergy = new h5_comp[self_dim];                           // create buffer for self-energy
        irreducible_class = new h5_comp[irred_dim];              // create buffer for irreducible vertex
#if MAX_DIAG_CLASS >= 1
        K1_class_a = new h5_comp[K1_dim];                        // create buffer for K1_a
        K1_class_p = new h5_comp[K1_dim];                        // create buffer for K1_p
        K1_class_t = new h5_comp[K1_dim];                        // create buffer for K1_t
#endif
#if MAX_DIAG_CLASS >= 2
        bfreqs2_buffer = new double[nBOS2];                      // create buffer for bosonic frequencies for K2
        ffreqs2_buffer = new double[nFER2];                      // create buffer for fermionic frequencies for K2
        K2_class_a = new h5_comp[K2_dim];                        // create buffer for K2_a
        K2_class_p = new h5_comp[K2_dim];                        // create buffer for K2_p
        K2_class_t = new h5_comp[K2_dim];                        // create buffer for K2_t
#endif
#if MAX_DIAG_CLASS >= 3
        bfreqs3_buffer = new double[nBOS3];                      // create buffer for bosonic frequencies for K3
        ffreqs3_buffer = new double[nFER3];                      // create buffer for fermionic frequencies for K3
        K3_class_a = new h5_comp[K3_dim];                        // create buffer for K3_a
        K3_class_p = new h5_comp[K3_dim];                        // create buffer for K3_p
        K3_class_t = new h5_comp[K3_dim];                        // create buffer for K3_t
#endif
    }

    ~Buffer() {
        delete[] freq_params;
        delete[] bfreqs_buffer;
        delete[] ffreqs_buffer;
        delete[] selfenergy;
        delete[] irreducible_class;
#if MAX_DIAG_CLASS >= 1
        delete[] K1_class_a;
        delete[] K1_class_p;
        delete[] K1_class_t;
#endif
#if MAX_DIAG_CLASS >= 2
        delete[] bfreqs2_buffer;
        delete[] ffreqs2_buffer;
        delete[] K2_class_a;
        delete[] K2_class_p;
        delete[] K2_class_t;
#endif
#if MAX_DIAG_CLASS >= 3
        delete[] bfreqs3_buffer;
        delete[] ffreqs3_buffer;
        delete[] K3_class_a;
        delete[] K3_class_p;
        delete[] K3_class_t;
#endif
    }

template <typename Q>
    void initialize(const State<Q>& state_in) {
        print("Starting to copy to buffer...", true);
        FrequencyGrid bfreqs = state_in.vertex[0].avertex().frequencies.b_K1;
        FrequencyGrid ffreqs = state_in.selfenergy.frequencies;
        freq_params[0] = (double) bfreqs.N_w;
        freq_params[1] = bfreqs.w_upper;
        freq_params[2] = bfreqs.w_lower;
        freq_params[3] = bfreqs.W_scale;
        freq_params[4] = (double) ffreqs.N_w;
        freq_params[5] = ffreqs.w_upper;
        freq_params[6] = ffreqs.w_lower;
        freq_params[7] = ffreqs.W_scale;

        for (int i=0; i<nBOS; ++i) {
            bfreqs_buffer[i] = bfreqs.ws[i];
        }
        for (int i=0; i<nFER; ++i) {
            ffreqs_buffer[i] = ffreqs.ws[i];
        }
        for (int i=0; i<self_dim; ++i) {                        // write self-energy into buffer
#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            // in the particle-hole symmetric case in Matsubara we only save the imaginary part of the selfenergy
            selfenergy[i].re = glb_U/2.;
            selfenergy[i].im = state_in.selfenergy.acc(i);
#else
            selfenergy[i].re = std::real(state_in.selfenergy.acc(i));
            selfenergy[i].im = std::imag(state_in.selfenergy.acc(i));
#endif
        }
        for (int i=0; i<irred_dim; ++i) {                       // write irreducible vertex into buffer
            irreducible_class[i].re = std::real(state_in.vertex[0].irred().acc(i));
            irreducible_class[i].im = std::imag(state_in.vertex[0].irred().acc(i));
        }
#if MAX_DIAG_CLASS >= 1
        for(int i=0; i<K1_dim; ++i) {                                // write K1 into buffer
            K1_class_a[i].re = std::real(state_in.vertex[0].avertex().K1_acc(i));
            K1_class_a[i].im = std::imag(state_in.vertex[0].avertex().K1_acc(i));

            K1_class_p[i].re = std::real(state_in.vertex[0].pvertex().K1_acc(i));
            K1_class_p[i].im = std::imag(state_in.vertex[0].pvertex().K1_acc(i));

            K1_class_t[i].re = std::real(state_in.vertex[0].tvertex().K1_acc(i));
            K1_class_t[i].im = std::imag(state_in.vertex[0].tvertex().K1_acc(i));
        }
#endif
#if MAX_DIAG_CLASS >= 2
        FrequencyGrid bfreqs2 = state_in.vertex[0].avertex().frequencies.b_K2;
        FrequencyGrid ffreqs2 = state_in.vertex[0].avertex().frequencies.f_K2;
        freq_params[8]  = (double) bfreqs2.N_w;
        freq_params[9]  = bfreqs2.w_upper;
        freq_params[10] = bfreqs2.w_lower;
        freq_params[11] = bfreqs2.W_scale;
        freq_params[12] = (double) ffreqs2.N_w;
        freq_params[13] = ffreqs2.w_upper;
        freq_params[14] = ffreqs2.w_lower;
        freq_params[15] = ffreqs2.W_scale;

        for (int i=0; i<nBOS2; ++i) {
            bfreqs2_buffer[i] = bfreqs2.ws[i];
        }
        for (int i=0; i<nFER2; ++i) {
            ffreqs2_buffer[i] = ffreqs2.ws[i];
        }
        for(int i=0; i<K2_dim; ++i) {                                // write K2 into buffer
            K2_class_a[i].re = std::real(state_in.vertex[0].avertex().K2_acc(i));
            K2_class_a[i].im = std::imag(state_in.vertex[0].avertex().K2_acc(i));

            K2_class_p[i].re = std::real(state_in.vertex[0].pvertex().K2_acc(i));
            K2_class_p[i].im = std::imag(state_in.vertex[0].pvertex().K2_acc(i));

            K2_class_t[i].re = std::real(state_in.vertex[0].tvertex().K2_acc(i));
            K2_class_t[i].im = std::imag(state_in.vertex[0].tvertex().K2_acc(i));
        }
#endif
#if MAX_DIAG_CLASS >= 3
        FrequencyGrid bfreqs3 = state_in.vertex[0].avertex().frequencies.b_K3;
        FrequencyGrid ffreqs3 = state_in.vertex[0].avertex().frequencies.f_K3;
        freq_params[16] = (double) bfreqs3.N_w;
        freq_params[17] = bfreqs3.w_upper;
        freq_params[18] = bfreqs3.w_lower;
        freq_params[19] = bfreqs3.W_scale;
        freq_params[20] = (double) ffreqs3.N_w;
        freq_params[21] = ffreqs3.w_upper;
        freq_params[22] = ffreqs3.w_lower;
        freq_params[23] = ffreqs3.W_scale;

        for (int i=0; i<nBOS3; ++i) {
            bfreqs3_buffer[i] = bfreqs3.ws[i];
        }
        for (int i=0; i<nFER3; ++i) {
            ffreqs3_buffer[i] = ffreqs3.ws[i];
        }
        for(int i=0; i<K3_dim; ++i) {                                // write K3 into buffer
            K3_class_a[i].re = std::real(state_in.vertex[0].avertex().K3_acc(i));
            K3_class_a[i].im = std::imag(state_in.vertex[0].avertex().K3_acc(i));

            K3_class_p[i].re = std::real(state_in.vertex[0].pvertex().K3_acc(i));
            K3_class_p[i].im = std::imag(state_in.vertex[0].pvertex().K3_acc(i));

            K3_class_t[i].re = std::real(state_in.vertex[0].tvertex().K3_acc(i));
            K3_class_t[i].im = std::imag(state_in.vertex[0].tvertex().K3_acc(i));
        }
#endif
        print("Buffer ready. Preparing for saving into Hdf5 file...", true);
    }
};

// Wrapper for static cast from int to hsize_t, used in class Dims
hsize_t h5_cast(int dim) {
    return static_cast<hsize_t>(dim);
}

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
    H5::DataSet *bfreqs_p, *ffreqs_p, *params_p;
#if MAX_DIAG_CLASS >= 1
    H5::DataSet *K1_a_p, *K1_p_p, *K1_t_p;
#endif
#if MAX_DIAG_CLASS >= 2
    H5::DataSet *bfreqs2_p, *ffreqs2_p;
    H5::DataSet *K2_a_p, *K2_p_p, *K2_t_p;
#endif
#if MAX_DIAG_CLASS >= 3
    H5::DataSet *bfreqs3_p, *ffreqs3_p;
    H5::DataSet *K3_a_p, *K3_p_p, *K3_t_p;
#endif

    // Data sets from existing file: Objects
    H5::DataSet lambda, self, irred;
    H5::DataSet freq_params;
    H5::DataSet bfreqs_dataset, ffreqs_dataset;
#if MAX_DIAG_CLASS >= 1
    H5::DataSet K1_a, K1_p, K1_t;
#endif
#if MAX_DIAG_CLASS >= 2
    H5::DataSet bfreqs2_dataset, ffreqs2_dataset;
    H5::DataSet K2_a, K2_p, K2_t;
#endif
#if MAX_DIAG_CLASS >= 3
    H5::DataSet bfreqs3_dataset, ffreqs3_dataset;
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
             H5::DataSpace& dataSpaces_bfreqs,
             H5::DataSpace& dataSpaces_ffreqs,
             H5::DataSpace& dataSpaces_params,
#if MAX_DIAG_CLASS >= 1
             H5::DataSpace& dataSpaces_K1_a,
             H5::DataSpace& dataSpaces_K1_p,
             H5::DataSpace& dataSpaces_K1_t,
#endif
#if MAX_DIAG_CLASS >= 2
             H5::DataSpace& dataSpaces_bfreqs2,
             H5::DataSpace& dataSpaces_ffreqs2,
             H5::DataSpace& dataSpaces_K2_a,
             H5::DataSpace& dataSpaces_K2_p,
             H5::DataSpace& dataSpaces_K2_t,
#endif
#if MAX_DIAG_CLASS >= 3
             H5::DataSpace& dataSpaces_bfreqs3,
             H5::DataSpace& dataSpaces_ffreqs3,
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
            bfreqs_p = new H5::DataSet(
                    file->createDataSet(BFREQS_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs));
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
            bfreqs2_p = new H5::DataSet(
                    file->createDataSet(BFREQS2_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs2));
            ffreqs2_p = new H5::DataSet(
                    file->createDataSet(FFREQS2_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs2));
            K2_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_a, mtype_comp, dataSpaces_K2_a, plist_vert));
            K2_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_p, mtype_comp, dataSpaces_K2_p, plist_vert));
            K2_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_t, mtype_comp, dataSpaces_K2_t, plist_vert));
#endif
#if MAX_DIAG_CLASS >= 3
            // Create the datasets in file:
            bfreqs3_p = new H5::DataSet(
                    file->createDataSet(BFREQS3_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs3));
            ffreqs3_p = new H5::DataSet(
                    file->createDataSet(FFREQS3_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs3));
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
            bfreqs_dataset = file->openDataSet("bfreqs");
            ffreqs_dataset = file->openDataSet("ffreqs");
#if MAX_DIAG_CLASS >=1
            K1_a = file->openDataSet("K1_a");
            K1_p = file->openDataSet("K1_p");
            K1_t = file->openDataSet("K1_t");
#endif
#if MAX_DIAG_CLASS >=2
            bfreqs2_dataset = file->openDataSet("bfreqs2");
            ffreqs2_dataset = file->openDataSet("ffreqs2");
            K2_a = file->openDataSet("K2_a");
            K2_p = file->openDataSet("K2_p");
            K2_t = file->openDataSet("K2_t");
#endif
#if MAX_DIAG_CLASS >=3
            bfreqs3_dataset = file->openDataSet("bfreqs3");
            ffreqs3_dataset = file->openDataSet("ffreqs3");
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
        bfreqs_dataset = file->openDataSet("bfreqs");
        ffreqs_dataset = file->openDataSet("ffreqs");
#if MAX_DIAG_CLASS >=1
        K1_a = file->openDataSet("K1_a");
        K1_p = file->openDataSet("K1_p");
        K1_t = file->openDataSet("K1_t");
#endif
#if MAX_DIAG_CLASS >=2
        bfreqs2_dataset = file->openDataSet("bfreqs2");
        ffreqs2_dataset = file->openDataSet("ffreqs2");
        K2_a = file->openDataSet("K2_a");
        K2_p = file->openDataSet("K2_p");
        K2_t = file->openDataSet("K2_t");
#endif
#if MAX_DIAG_CLASS >=3
        bfreqs3_dataset = file->openDataSet("bfreqs3");
        ffreqs3_dataset = file->openDataSet("ffreqs3");
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
            bfreqs_p -> close();
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
            bfreqs2_p -> close();
            ffreqs2_p -> close();
            K2_a_p -> close();
            K2_p_p -> close();
            K2_t_p -> close();
#endif
#if MAX_DIAG_CLASS >= 3
            bfreqs3_p -> close();
            ffreqs3_p -> close();
            K3_a_p -> close();
            K3_p_p -> close();
            K3_t_p -> close();
#endif
        }
        else {
            self.close();
            irred.close();
            freq_params.close();
            bfreqs_dataset.close();
            ffreqs_dataset.close();

#if MAX_DIAG_CLASS >=1
            K1_a.close();
            K1_p.close();
            K1_t.close();
#endif
#if MAX_DIAG_CLASS >=2
            bfreqs2_dataset.close();
            ffreqs2_dataset.close();
            K2_a.close();
            K2_p.close();
            K2_t.close();
#endif
#if MAX_DIAG_CLASS >=3
            bfreqs3_dataset.close();
            ffreqs3_dataset.close();
            K3_a.close();
            K3_p.close();
            K3_t.close();
#endif
        }
    }
};

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
                 const State<Q>& state_in, rvec& Lambdas, bool file_exists) {
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
        H5::DataSpace dataSpaces_bfreqs(RANK_freqs, dims.bfreqs_dims);
        H5::DataSpace dataSpaces_ffreqs(RANK_freqs, dims.ffreqs_dims);
        H5::DataSpace dataSpaces_params(1, dims.params);
        H5::DataSpace dataSpaces_selfenergy(RANK_self, dims.selfenergy);
        H5::DataSpace dataSpaces_irreducible(RANK_irreducible, dims.irreducible);

        H5::DataSpace dataSpaces_freq_params_buffer(RANK_freqs-1, dims.freq_params_buffer_dims);
        H5::DataSpace dataSpaces_bfreqs_buffer(RANK_freqs-1, dims.bfreqs_buffer_dims);
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
        H5::DataSpace dataSpaces_bfreqs2(RANK_freqs, dims.bfreqs2_dims);
        H5::DataSpace dataSpaces_ffreqs2(RANK_freqs, dims.ffreqs2_dims);

        H5::DataSpace dataSpaces_bfreqs2_buffer(RANK_freqs-1, dims.bfreqs2_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs2_buffer(RANK_freqs-1, dims.ffreqs2_buffer_dims);

        H5::DataSpace dataSpaces_K2_a(RANK_K2, dims.K2);
        H5::DataSpace dataSpaces_K2_p(RANK_K2, dims.K2);
        H5::DataSpace dataSpaces_K2_t(RANK_K2, dims.K2);

        H5::DataSpace dataSpaces_K2_a_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_p_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_t_buffer(RANK_K2-1, dims.K2_buffer);
#endif
#if MAX_DIAG_CLASS >= 3
        H5::DataSpace dataSpaces_bfreqs3(RANK_freqs, dims.bfreqs3_dims);
        H5::DataSpace dataSpaces_ffreqs3(RANK_freqs, dims.ffreqs3_dims);

        H5::DataSpace dataSpaces_bfreqs3_buffer(RANK_freqs-1, dims.bfreqs3_buffer_dims);
        H5::DataSpace dataSpaces_ffreqs3_buffer(RANK_freqs-1, dims.ffreqs3_buffer_dims);

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
                          dataSpaces_bfreqs, dataSpaces_ffreqs, dataSpaces_params,
#if MAX_DIAG_CLASS >= 1
                          dataSpaces_K1_a, dataSpaces_K1_p, dataSpaces_K1_t,
#endif
#if MAX_DIAG_CLASS >= 2
                          dataSpaces_bfreqs2, dataSpaces_ffreqs2,
                          dataSpaces_K2_a, dataSpaces_K2_p, dataSpaces_K2_t,
#endif
#if MAX_DIAG_CLASS >= 3
                          dataSpaces_bfreqs3, dataSpaces_ffreqs3,
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
        dataSpaces_bfreqs.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs_p -> write(buffer.bfreqs_buffer, H5::PredType::NATIVE_DOUBLE,
                                       dataSpaces_bfreqs_buffer, dataSpaces_bfreqs);
        else
            dataSets.bfreqs_dataset.write(buffer.bfreqs_buffer, H5::PredType::NATIVE_DOUBLE,
                                          dataSpaces_bfreqs_buffer, dataSpaces_bfreqs);

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
        dataSpaces_bfreqs2.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs2_p -> write(buffer.bfreqs2_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs2_buffer, dataSpaces_bfreqs2);
        else
            dataSets.bfreqs2_dataset.write(buffer.bfreqs2_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs2_buffer, dataSpaces_bfreqs2);

        count[1] = nFER2;
        dataSpaces_ffreqs2.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs2_p -> write(buffer.ffreqs2_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs2_buffer, dataSpaces_ffreqs2);
        else
            dataSets.ffreqs2_dataset.write(buffer.ffreqs2_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs2_buffer, dataSpaces_ffreqs2);

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
        dataSpaces_bfreqs3.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.bfreqs3_p -> write(buffer.bfreqs3_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_bfreqs3_buffer, dataSpaces_bfreqs3);
        else
            dataSets.bfreqs3_dataset.write(buffer.bfreqs3_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_bfreqs3_buffer, dataSpaces_bfreqs3);

        count[1] = nFER3;
        dataSpaces_ffreqs3.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.ffreqs3_p -> write(buffer.ffreqs3_buffer, H5::PredType::NATIVE_DOUBLE,
                                        dataSpaces_ffreqs3_buffer, dataSpaces_ffreqs3);
        else
            dataSets.ffreqs3_dataset.write(buffer.ffreqs3_buffer, H5::PredType::NATIVE_DOUBLE,
                                           dataSpaces_ffreqs3_buffer, dataSpaces_ffreqs3);

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

        print("Successfully saved in hdf5 file: ", FILE_NAME);
        if (file_exists) print_add(" in Lambda-layer ", Lambda_it, false);
        print_add("", true);

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
#ifdef MPI_FLAG
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
void add_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size,
             State<Q>& state_in, rvec& Lambdas) {
#ifdef MPI_FLAG
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        // write data to file if Lambda iteration number is in allowed range, otherwise print error message
        if (Lambda_it < Lambda_size) {
            save_to_hdf(FILE_NAME, Lambda_it, Lambda_size, state_in, Lambdas, true);
        } else {
            print("Cannot write to file ", FILE_NAME, " since Lambda layer is out of range.", true);
        }
    }
}

/** overload of add_hdf for non-States, does not do anything */
template <typename Q>
void add_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size,
             Q& state_in, rvec& Lambdas) {}


/// --- Functions for reading data from file --- ///

/**
 * Initialize the frequency grids of the result using the parameters stored in buffer.
 * @param result : Empty State object into which result is copied.
 * @param buffer : Buffer from which result is read. Should contain data read from a file.
 */
template <typename Q>
void result_set_frequency_grids(State<Q>& result, Buffer& buffer) {
    // create new frequency grids
    FrequencyGrid bfreqs ('b', 1, Lambda_ini);
    FrequencyGrid ffreqs ('f', 1, Lambda_ini);
    // read grid parameters from buffer
    bfreqs.N_w = (int)buffer.freq_params[0];
    bfreqs.w_upper = buffer.freq_params[1];
    bfreqs.w_lower = buffer.freq_params[2];
    bfreqs.W_scale = buffer.freq_params[3];
    ffreqs.N_w = (int)buffer.freq_params[4];
    ffreqs.w_upper = buffer.freq_params[5];
    ffreqs.w_lower = buffer.freq_params[6];
    ffreqs.W_scale = buffer.freq_params[7];
    // initialize grids
    bfreqs.initialize_grid();
    ffreqs.initialize_grid();
    // copy grids to result
    result.selfenergy.frequencies = ffreqs;
    result.vertex[0].avertex().frequencies.b_K1 = bfreqs;
    result.vertex[0].pvertex().frequencies.b_K1 = bfreqs;
    result.vertex[0].tvertex().frequencies.b_K1 = bfreqs;
#if MAX_DIAG_CLASS >= 2
    FrequencyGrid bfreqs2 ('b', 2, Lambda_ini);
    FrequencyGrid ffreqs2 ('f', 2, Lambda_ini);
    bfreqs2.N_w = (int)buffer.freq_params[8];
    bfreqs2.w_upper = buffer.freq_params[9];
    bfreqs2.w_lower = buffer.freq_params[10];
    bfreqs2.W_scale = buffer.freq_params[11];
    ffreqs2.N_w = (int)buffer.freq_params[12];
    ffreqs2.w_upper = buffer.freq_params[13];
    ffreqs2.w_lower = buffer.freq_params[14];
    ffreqs2.W_scale = buffer.freq_params[15];
    bfreqs2.initialize_grid();
    ffreqs2.initialize_grid();
    result.vertex[0].avertex().frequencies.b_K2 = bfreqs2;
    result.vertex[0].pvertex().frequencies.b_K2 = bfreqs2;
    result.vertex[0].tvertex().frequencies.b_K2 = bfreqs2;
    result.vertex[0].avertex().frequencies.f_K2 = ffreqs2;
    result.vertex[0].pvertex().frequencies.f_K2 = ffreqs2;
    result.vertex[0].tvertex().frequencies.f_K2 = ffreqs2;
#endif
#if MAX_DIAG_CLASS >= 3
    FrequencyGrid bfreqs3 ('b', 3, Lambda_ini);
    FrequencyGrid ffreqs3 ('f', 3, Lambda_ini);
    bfreqs3.N_w = (int)buffer.freq_params[16];
    bfreqs3.w_upper = buffer.freq_params[17];
    bfreqs3.w_lower = buffer.freq_params[18];
    bfreqs3.W_scale = buffer.freq_params[19];
    ffreqs3.N_w = (int)buffer.freq_params[20];
    ffreqs3.w_upper = buffer.freq_params[21];
    ffreqs3.w_lower = buffer.freq_params[22];
    ffreqs3.W_scale = buffer.freq_params[23];
    bfreqs3.initialize_grid();
    ffreqs3.initialize_grid();
    result.vertex[0].avertex().frequencies.b_K3 = bfreqs3;
    result.vertex[0].pvertex().frequencies.b_K3 = bfreqs3;
    result.vertex[0].tvertex().frequencies.b_K3 = bfreqs3;
    result.vertex[0].avertex().frequencies.f_K3 = ffreqs3;
    result.vertex[0].pvertex().frequencies.f_K3 = ffreqs3;
    result.vertex[0].tvertex().frequencies.f_K3 = ffreqs3;
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
        result.vertex[0].irred().direct_set(i, val);
    }
#if MAX_DIAG_CLASS >= 1
    for (int i=0; i<buffer.K1_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K1_class_a[i].re, buffer.K1_class_a[i].im};
#else
        val = buffer.K1_class_a[i].re;
#endif
        result.vertex[0].avertex().K1_direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K1_class_p[i].re, buffer.K1_class_p[i].im};
#else
        val = buffer.K1_class_p[i].re;
#endif
        result.vertex[0].pvertex().K1_direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K1_class_t[i].re, buffer.K1_class_t[i].im};
#else
        val = buffer.K1_class_t[i].re;
#endif
        result.vertex[0].tvertex().K1_direct_set(i, val);
    }
#endif
#if MAX_DIAG_CLASS >= 2
    for (int i=0; i<buffer.K2_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K2_class_a[i].re, buffer.K2_class_a[i].im};
#else
        val = buffer.K2_class_a[i].re;
#endif
        result.vertex[0].avertex().K2_direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K2_class_p[i].re, buffer.K2_class_p[i].im};
#else
        val = buffer.K2_class_p[i].re;
#endif
        result.vertex[0].pvertex().K2_direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K2_class_t[i].re, buffer.K2_class_t[i].im};
#else
        val = buffer.K2_class_t[i].re;
#endif
        result.vertex[0].tvertex().K2_direct_set(i, val);
    }
#endif
#if MAX_DIAG_CLASS >= 3
    for (int i=0; i<buffer.K3_dim; ++i) {
#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K3_class_a[i].re, buffer.K3_class_a[i].im};
#else
        val = buffer.K3_class_a[i].re;
#endif
        result.vertex[0].avertex().K3_direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K3_class_p[i].re, buffer.K3_class_p[i].im};
#else
        val = buffer.K3_class_p[i].re;
#endif
        result.vertex[0].pvertex().K3_direct_set(i, val);

#if defined(KELDYSH_FORMALISM) or not defined(PARTICLE_HOLE_SYMM)
        val = {buffer.K3_class_t[i].re, buffer.K3_class_t[i].im};
#else
        val = buffer.K3_class_t[i].re;
#endif
        result.vertex[0].tvertex().K3_direct_set(i, val);
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
State<state_datatype> read_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size){
    State<state_datatype> result(Lambda_ini);
    if (Lambda_it < Lambda_size) {

        // Open the file. Access rights: read-only
        H5::H5File *file = 0;
        file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);

        // Prepare a buffer to which data from file is written. Buffer is copied to result eventually.
        Buffer buffer;

        // Create the memory data type for storing complex numbers in file
        H5::CompType mtype_comp = def_mtype_comp();

        // Create the dimension arrays for objects in file and in buffer
        Dims dims(buffer, Lambda_size);

        // Read the data sets from the file to copy their content into the buffer
        DataSets dataSets(file);

        // Create the data spaces for the data sets in file and for buffer objects
        H5::DataSpace dataSpaces_Lambda = dataSets.lambda.getSpace();
        H5::DataSpace dataSpaces_freq_params;
        try {   // storing frequency gri parameters was implemented later --> old files do not have it
            dataSpaces_freq_params = dataSets.freq_params.getSpace();
        }
        catch (H5::DataSetIException error) {
            error.printErrorStack();
        }
        H5::DataSpace dataSpaces_selfenergy = dataSets.self.getSpace();
        H5::DataSpace dataSpaces_irreducible = dataSets.irred.getSpace();

        H5::DataSpace dataSpaces_Lambda_buffer(1, dims.Lambda);
        H5::DataSpace dataSpaces_freq_params_buffer(RANK_freqs-1, dims.freq_params_buffer_dims);
        H5::DataSpace dataSpaces_selfenergy_buffer(RANK_self-1, dims.selfenergy_buffer);
        H5::DataSpace dataSpaces_irreducible_buffer(RANK_irreducible-1, dims.irreducible_buffer);

#if MAX_DIAG_CLASS >= 1
        H5::DataSpace dataSpaces_K1_a = dataSets.K1_a.getSpace();
        H5::DataSpace dataSpaces_K1_p = dataSets.K1_p.getSpace();
        H5::DataSpace dataSpaces_K1_t = dataSets.K1_t.getSpace();

        H5::DataSpace dataSpaces_K1_a_buffer(RANK_K1-1, dims.K1_buffer);
        H5::DataSpace dataSpaces_K1_p_buffer(RANK_K1-1, dims.K1_buffer);
        H5::DataSpace dataSpaces_K1_t_buffer(RANK_K1-1, dims.K1_buffer);
#endif

#if MAX_DIAG_CLASS >= 2
        H5::DataSpace dataSpaces_K2_a = dataSets.K2_a.getSpace();
        H5::DataSpace dataSpaces_K2_p = dataSets.K2_p.getSpace();
        H5::DataSpace dataSpaces_K2_t = dataSets.K2_t.getSpace();

        H5::DataSpace dataSpaces_K2_a_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_p_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_t_buffer(RANK_K2-1, dims.K2_buffer);
#endif

#if MAX_DIAG_CLASS >= 3
        H5::DataSpace dataSpaces_K3_a = dataSets.K3_a.getSpace();
        H5::DataSpace dataSpaces_K3_p = dataSets.K3_p.getSpace();
        H5::DataSpace dataSpaces_K3_t = dataSets.K3_t.getSpace();

        H5::DataSpace dataSpaces_K3_a_buffer(RANK_K3-1, dims.K3_buffer);
        H5::DataSpace dataSpaces_K3_p_buffer(RANK_K3-1, dims.K3_buffer);
        H5::DataSpace dataSpaces_K3_t_buffer(RANK_K3-1, dims.K3_buffer);
#endif


        //Select hyperslab in the file where the data should be located
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

        /// load data into buffer

        count[1] = 1;
        dataSpaces_Lambda.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSets.lambda.read(buffer.lambda, mtype_comp,
                           dataSpaces_Lambda_buffer, dataSpaces_Lambda);


        count[1] = N_freq_params;
        try {   // storing frequency gri parameters was implemented later --> old files do not have it
            dataSpaces_freq_params.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSets.freq_params.read(buffer.freq_params, H5::PredType::NATIVE_DOUBLE,
                                      dataSpaces_freq_params_buffer, dataSpaces_freq_params);
        }
        catch (H5::DataSpaceIException error) {
            error.printErrorStack();
        }

        count[1] = buffer.self_dim;
        dataSpaces_selfenergy.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSets.self.read(buffer.selfenergy, mtype_comp,
                           dataSpaces_selfenergy_buffer, dataSpaces_selfenergy);

        count[1] = buffer.irred_dim;
        dataSpaces_irreducible.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSets.irred.read(buffer.irreducible_class, mtype_comp,
                            dataSpaces_irreducible_buffer, dataSpaces_irreducible);


#if MAX_DIAG_CLASS >= 1
        count[1] = buffer.K1_dim;
        dataSpaces_K1_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K1_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K1_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataSets.K1_a.read(buffer.K1_class_a, mtype_comp, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
        dataSets.K1_p.read(buffer.K1_class_p, mtype_comp, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
        dataSets.K1_t.read(buffer.K1_class_t, mtype_comp, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
#endif
#if MAX_DIAG_CLASS >= 2
        count[1] = buffer.K2_dim;
        dataSpaces_K2_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K2_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K2_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataSets.K2_a.read(buffer.K2_class_a, mtype_comp, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
        dataSets.K2_p.read(buffer.K2_class_p, mtype_comp, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
        dataSets.K2_t.read(buffer.K2_class_t, mtype_comp, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
#endif
#if MAX_DIAG_CLASS >= 3
        count[1] = buffer.K3_dim;
        dataSpaces_K3_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K3_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K3_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataSets.K3_a.read(buffer.K3_class_a, mtype_comp, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
        dataSets.K3_p.read(buffer.K3_class_p, mtype_comp, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
        dataSets.K3_t.read(buffer.K3_class_t, mtype_comp, dataSpaces_K3_t_buffer, dataSpaces_K3_t);

#endif

        // Initialize the frequency grids of the result State using the parameters stored in buffer
        result_set_frequency_grids(result, buffer);

        // Copy the buffered result into State object
        copy_buffer_to_result(result, buffer);

        // Terminate
        dataSets.close(true);
        file->close();
        delete file;

        return result;

    } else {
        print("Cannot read from file ", FILE_NAME, " since Lambda layer out of range", true);
    }
}

/// --- Test function --- ///

void test_hdf5(H5std_string FILE_NAME, int i, State<state_datatype>& state) {
    // test hdf5: read files and compare to original file
    int cnt = 0;
    State<state_datatype> out = read_hdf(FILE_NAME, i, nODE);
    for (int iK=0; iK<2; ++iK) {
        for (int iSE = 0; iSE < nSE; ++iSE) {
            if (state.selfenergy.val(iK, iSE, 0) != out.selfenergy.val(iK, iSE, 0)) {
                std::cout << "Self-energy not equal, " << iK << ", " << iSE << std::endl;
                cnt += 1;
            }
        }
    }

    for (int iK=0; iK<nK_K1; ++iK) {
        for (int i_in=0; i_in<n_in; ++i_in) {
#if MAX_DIAG_CLASS >= 1
            for (int iw1=0; iw1<nBOS; ++iw1) {
                if (state.vertex[0].avertex().K1_val(iK, iw1, i_in) != out.vertex[0].avertex().K1_val(iK, iw1, i_in)) {
                    std::cout << "Vertex not equal, " << iK << ", " << iw1 << std::endl;
                    cnt += 1;
                }
                if (state.vertex[0].pvertex().K1_val(iK, iw1, i_in) != out.vertex[0].pvertex().K1_val(iK, iw1, i_in)) {
                    std::cout << "Vertex not equal, " << iK << ", " << iw1 << std::endl;
                    cnt += 1;
                }
                if (state.vertex[0].tvertex().K1_val(iK, iw1, i_in) != out.vertex[0].tvertex().K1_val(iK, iw1, i_in)) {
                    std::cout << "Vertex not equal, " << iK << ", " << iw1 << std::endl;
                    cnt += 1;
                }
#if MAX_DIAG_CLASS >= 2
                for (int iw2=0; iw2<nFER; ++iw2) {
                    if (state.vertex[0].avertex().K2_val(iK, iw1, iw2, i_in) != out.vertex[0].avertex().K2_val(iK, iw1, iw2, i_in)) {
                        std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << std::endl;
                        cnt += 1;
                    }
                    if (state.vertex[0].pvertex().K2_val(iK, iw1, iw2, i_in) != out.vertex[0].pvertex().K2_val(iK, iw1, iw2, i_in)) {
                        std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << std::endl;
                        cnt += 1;
                    }
                    if (state.vertex[0].tvertex().K2_val(iK, iw1, iw2, i_in) != out.vertex[0].tvertex().K2_val(iK, iw1, iw2, i_in)) {
                        std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << std::endl;
                        cnt += 1;
                    }
#if MAX_DIAG_CLASS == 3
                    for (int iw3=0; iw3<nFER; ++iw3) {
                        if (state.vertex[0].avertex().K3_val(iK, iw1, iw2, iw3, i_in) != out.vertex[0].avertex().K3_val(iK, iw1, iw2, iw3, i_in)) {
                            std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << std::endl;
                            cnt += 1;
                        }
                        if (state.vertex[0].pvertex().K3_val(iK, iw1, iw2, iw3, i_in) != out.vertex[0].pvertex().K3_val(iK, iw1, iw2, iw3, i_in)) {
                            std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << std::endl;
                            cnt += 1;
                        }
                        if (state.vertex[0].tvertex().K3_val(iK, iw1, iw2, iw3, i_in) != out.vertex[0].tvertex().K3_val(iK, iw1, iw2, iw3, i_in)) {
                            std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << std::endl;
                            cnt += 1;
                        }
                    }
#endif
                }
#endif
            }
#endif
        }
    }
    if (cnt == 0) print("HDF5 test successful.", true);
    else print("HDF5 test failed. Number of differences: ", cnt, true);
}

#endif //KELDYSH_MFRG_HDF5_ROUTINES_H