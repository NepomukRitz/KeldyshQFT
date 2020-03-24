// TODO: include mpi_world_rank check
// TODO: template Q ??

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <fstream>
#include <type_traits>
#include <string>
#include <cstdlib>
#include "H5Cpp.h"

#include "util.h"               // printing text
#include "parameters.h"         // system parameters (necessary for vector lengths etc.)
#include "data_structures.h"    // comp data type, real/complex vector class

#ifdef MPI_FLAG
#include "mpi_setup.h"
#endif

/********************************constants concerning HDF5 data format*************************/

// dataset dimensions
const int RANK_K1 = 2;
const int RANK_K2 = 2;
const int RANK_K3 = 2;
const int RANK_irreducible = 2;
const int RANK_self = 2;

// names of the individual datasets within the hdf5 file
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
const H5std_string  BFREQS_LIST("bfreqs");
const H5std_string  FFREQS_LIST("ffreqs");
const H5std_string  PARAM_LIST("parameters");
const H5std_string  MEMBER1( "spin_re" );
const H5std_string  MEMBER2( "spin_im" );
const H5std_string  MEMBER3( "dens_re" );
const H5std_string  MEMBER4( "dens_im" );
const H5std_string  MEMBER5( "re" );
const H5std_string  MEMBER6( "im" );


// define struct to save complex numbers with two spin components in hdf5 file
typedef struct h5_comp_spin {
    double spin_re; // real part for spin component
    double spin_im; // imag part for spin component
    double dens_re; // real part for dens component
    double dens_im; // imag part for dens component
} h5_comp_spin;

// define struct to save complex numbers in hdf5 file
typedef struct h5_comp {
    double re; // real part
    double im; // imaginary part
} h5_comp;

// Create the memory datatype
H5::CompType def_mtype_comp() {
    H5::CompType mtype_comp(sizeof(h5_comp));
    mtype_comp.insertMember(MEMBER5, HOFFSET(h5_comp, re), H5::PredType::NATIVE_DOUBLE);
    mtype_comp.insertMember(MEMBER6, HOFFSET(h5_comp, im), H5::PredType::NATIVE_DOUBLE);
    return mtype_comp;
}
H5::CompType def_mtype_comp_spin() {
    H5::CompType mtype_comp_spin(sizeof(h5_comp_spin));
    mtype_comp_spin.insertMember(MEMBER1, HOFFSET(h5_comp_spin, spin_re), H5::PredType::NATIVE_DOUBLE);
    mtype_comp_spin.insertMember(MEMBER2, HOFFSET(h5_comp_spin, spin_im), H5::PredType::NATIVE_DOUBLE);
    mtype_comp_spin.insertMember(MEMBER3, HOFFSET(h5_comp_spin, dens_re), H5::PredType::NATIVE_DOUBLE);
    mtype_comp_spin.insertMember(MEMBER4, HOFFSET(h5_comp_spin, dens_im), H5::PredType::NATIVE_DOUBLE);
    return mtype_comp_spin;
}

class Buffer {
public:
    const int self_dim = 2 * nSE;                                     // length of self-energy buffer
    h5_comp * selfenergy;
    const int irred_dim = 16 * n_in;                                  // length of irreducible vertex buffer
    h5_comp_spin * irreducible_class;
#if DIAG_CLASS >= 1
    const int K1_dim = nK_K1 * nw1_wt * n_in;                         // length of K1 buffer
    h5_comp_spin * K1_class_a;
    h5_comp_spin * K1_class_p;
    h5_comp_spin * K1_class_t;
#endif
#if DIAG_CLASS >= 2
    const int K2_dim = nK_K2 * nw2_wt * nw2_nut * n_in;               // length of K2 buffer
    h5_comp_spin * K2_class_a;
    h5_comp_spin * K2_class_p;
    h5_comp_spin * K2_class_t;
#endif
#if DIAG_CLASS >= 3
    const int K3_dim = nK_K3 * nw3_wt * nw3_nut * nw3_nutp * n_in;    // length of K3 buffer
    h5_comp_spin * K3_class_a;
    h5_comp_spin * K3_class_p;
    h5_comp_spin * K3_class_t;
#endif

    Buffer() {
        selfenergy = new h5_comp[self_dim];                           // create buffer for self-energy
        irreducible_class = new h5_comp_spin[irred_dim];              // create buffer for irreducible vertex
#if DIAG_CLASS >= 1
        K1_class_a = new h5_comp_spin[K1_dim];                        // create buffer for K1_a
        K1_class_p = new h5_comp_spin[K1_dim];                        // create buffer for K1_p
        K1_class_t = new h5_comp_spin[K1_dim];                        // create buffer for K1_t
#endif
#if DIAG_CLASS >= 2
        K2_class_a = new h5_comp_spin[K2_dim];                        // create buffer for K2_a
        K2_class_p = new h5_comp_spin[K2_dim];                        // create buffer for K2_p
        K2_class_t = new h5_comp_spin[K2_dim];                        // create buffer for K2_t
#endif
#if DIAG_CLASS >= 3
        K3_class_a = new h5_comp_spin[K3_dim];                        // create buffer for K3_a
        K3_class_p = new h5_comp_spin[K3_dim];                        // create buffer for K3_p
        K3_class_t = new h5_comp_spin[K3_dim];                        // create buffer for K3_t
#endif
    }

    ~Buffer() {
        delete[] selfenergy;
        delete[] irreducible_class;
#if DIAG_CLASS >= 1
        delete[] K1_class_a;
        delete[] K1_class_p;
        delete[] K1_class_t;
#endif
#if DIAG_CLASS >= 2
        delete[] K2_class_a;
        delete[] K2_class_p;
        delete[] K2_class_t;
#endif
#if DIAG_CLASS >= 3
        delete[] K3_class_a;
        delete[] K3_class_p;
        delete[] K3_class_t;
#endif
    }

    void initialize(State<comp>& state_in) {
        print("Starting to copy to buffer...", true);
        for (int i = 0; i < self_dim; ++i) {                        // write self-energy into buffer
            selfenergy[i].re = real(state_in.selfenergy.acc(i));
            selfenergy[i].im = imag(state_in.selfenergy.acc(i));
        }
        for (int i = 0; i < irred_dim; ++i) {                       // write irreducible vertex into buffer
            irreducible_class[i].spin_re = real(state_in.vertex.spinvertex.irred.acc(i));
            irreducible_class[i].dens_re = real(state_in.vertex.densvertex.irred.acc(i));
            irreducible_class[i].spin_im = imag(state_in.vertex.spinvertex.irred.acc(i));
            irreducible_class[i].dens_im = imag(state_in.vertex.densvertex.irred.acc(i));
        }
#if DIAG_CLASS >= 1
        for(int i=0; i<K1_dim; ++i){                                // write K1 into buffer
            K1_class_a[i].spin_re = real(state_in.vertex.spinvertex.avertex.K1_acc(i));
            K1_class_a[i].dens_re = real(state_in.vertex.densvertex.avertex.K1_acc(i));
            K1_class_a[i].spin_im = imag(state_in.vertex.spinvertex.avertex.K1_acc(i));
            K1_class_a[i].dens_im = imag(state_in.vertex.densvertex.avertex.K1_acc(i));

            K1_class_p[i].spin_re = real(state_in.vertex.spinvertex.pvertex.K1_acc(i));
            K1_class_p[i].dens_re = real(state_in.vertex.densvertex.pvertex.K1_acc(i));
            K1_class_p[i].spin_im = imag(state_in.vertex.spinvertex.pvertex.K1_acc(i));
            K1_class_p[i].dens_im = imag(state_in.vertex.densvertex.pvertex.K1_acc(i));

            K1_class_t[i].spin_re = real(state_in.vertex.spinvertex.tvertex.K1_acc(i));
            K1_class_t[i].dens_re = real(state_in.vertex.densvertex.tvertex.K1_acc(i));
            K1_class_t[i].spin_im = imag(state_in.vertex.spinvertex.tvertex.K1_acc(i));
            K1_class_t[i].dens_im = imag(state_in.vertex.densvertex.tvertex.K1_acc(i));
        }
#endif
#if DIAG_CLASS >= 2
        for(int i=0; i<K2_dim; ++i){                                // write K2 into buffer
            K2_class_a[i].spin_re = real(state_in.vertex.spinvertex.avertex.K2_acc(i));
            K2_class_a[i].dens_re = real(state_in.vertex.densvertex.avertex.K2_acc(i));
            K2_class_a[i].spin_im = imag(state_in.vertex.spinvertex.avertex.K2_acc(i));
            K2_class_a[i].dens_im = imag(state_in.vertex.densvertex.avertex.K2_acc(i));

            K2_class_p[i].spin_re = real(state_in.vertex.spinvertex.pvertex.K2_acc(i));
            K2_class_p[i].dens_re = real(state_in.vertex.densvertex.pvertex.K2_acc(i));
            K2_class_p[i].spin_im = imag(state_in.vertex.spinvertex.pvertex.K2_acc(i));
            K2_class_p[i].dens_im = imag(state_in.vertex.densvertex.pvertex.K2_acc(i));

            K2_class_t[i].spin_re = real(state_in.vertex.spinvertex.tvertex.K2_acc(i));
            K2_class_t[i].dens_re = real(state_in.vertex.densvertex.tvertex.K2_acc(i));
            K2_class_t[i].spin_im = imag(state_in.vertex.spinvertex.tvertex.K2_acc(i));
            K2_class_t[i].dens_im = imag(state_in.vertex.densvertex.tvertex.K2_acc(i));
        }
#endif
#if DIAG_CLASS >= 3
        for(int i=0; i<K3_dim; ++i){                                // write K3 into buffer
            K3_class_a[i].spin_re = real(state_in.vertex.spinvertex.avertex.K3_acc(i));
            K3_class_a[i].dens_re = real(state_in.vertex.densvertex.avertex.K3_acc(i));
            K3_class_a[i].spin_im = imag(state_in.vertex.spinvertex.avertex.K3_acc(i));
            K3_class_a[i].dens_im = imag(state_in.vertex.densvertex.avertex.K3_acc(i));

            K3_class_p[i].spin_re = real(state_in.vertex.spinvertex.pvertex.K3_acc(i));
            K3_class_p[i].dens_re = real(state_in.vertex.densvertex.pvertex.K3_acc(i));
            K3_class_p[i].spin_im = imag(state_in.vertex.spinvertex.pvertex.K3_acc(i));
            K3_class_p[i].dens_im = imag(state_in.vertex.densvertex.pvertex.K3_acc(i));

            K3_class_t[i].spin_re = real(state_in.vertex.spinvertex.tvertex.K3_acc(i));
            K3_class_t[i].dens_re = real(state_in.vertex.densvertex.tvertex.K3_acc(i));
            K3_class_t[i].spin_im = imag(state_in.vertex.spinvertex.tvertex.K3_acc(i));
            K3_class_t[i].dens_im = imag(state_in.vertex.densvertex.tvertex.K3_acc(i));
        }
#endif
        print("Buffer ready. Preparing for saving into Hdf5 file...", true);
    }
};

hsize_t h5_cast(int dim) {
    return static_cast<hsize_t>(dim);
}

// TODO: initializer list: "," for #if conditions??

class Dims {
public:
    hsize_t Lambda[1];
    hsize_t bfreqs[1];
    hsize_t ffreqs[1];
    hsize_t params[1];
    hsize_t selfenergy[2];
    hsize_t selfenergy_buffer[1];
    hsize_t irreducible[2];
    hsize_t irreducible_buffer[1];
#if DIAG_CLASS >= 1
    hsize_t K1[2];
    hsize_t K1_buffer[1];
#endif
#if DIAG_CLASS >= 2
    hsize_t K2[2];
    hsize_t K2_buffer[1];
#endif
#if DIAG_CLASS >= 3
    hsize_t K3[2];
    hsize_t K3_buffer[1];
#endif

    Dims (Buffer& buffer, long Lambda_size) :
        // Create the dimension arrays for objects in file and in buffer
        Lambda {h5_cast(Lambda_size)},
        bfreqs {h5_cast(nBOS)},
        ffreqs {h5_cast(nFER)},
        params {h5_cast(param_size)},
        selfenergy {h5_cast(Lambda_size), h5_cast(buffer.self_dim)},
        selfenergy_buffer {h5_cast(buffer.self_dim)},
        irreducible {h5_cast(Lambda_size), h5_cast(buffer.irred_dim)},
        irreducible_buffer {h5_cast(buffer.irred_dim)},
#if DIAG_CLASS >= 1
        K1 {h5_cast(Lambda_size), h5_cast(buffer.K1_dim)},
        K1_buffer {h5_cast(buffer.K1_dim)},
#endif
#if DIAG_CLASS >= 2
        K2 {h5_cast(Lambda_size), h5_cast(buffer.K2_dim)},
        K2_buffer {h5_cast(buffer.K2_dim)},
#endif
#if DIAG_CLASS >= 3
        K3 {h5_cast(Lambda_size), h5_cast(buffer.K3_dim)},
        K3_buffer {h5_cast(buffer.K3_dim)}
#endif
    {}
};

struct DataSpaces { // Create the data space for the dataset in file and for buffer objects
    H5::DataSpace Lambda;
    H5::DataSpace bfreqs;
    H5::DataSpace ffreqs;
    H5::DataSpace params;
    H5::DataSpace selfenergy, selfenergy_buffer;
    H5::DataSpace irreducible, irreducible_buffer;
#if DIAG_CLASS >= 1
    H5::DataSpace K1_a, K1_a_buffer;
    H5::DataSpace K1_p, K1_p_buffer;
    H5::DataSpace K1_t, K1_t_buffer;
#endif
#if DIAG_CLASS >= 2
    H5::DataSpace K2_a, K2_a_buffer;
    H5::DataSpace K2_p, K2_p_buffer;
    H5::DataSpace K2_t, K2_t_buffer;
#endif
#if DIAG_CLASS >= 3
    H5::DataSpace K3_a, K3_a_buffer;
    H5::DataSpace K3_p, K3_p_buffer;
    H5::DataSpace K3_t, K3_t_buffer;
#endif
};

DataSpaces initialize_DataSpaces(Dims& dims) {
    DataSpaces dataSpaces;

    dataSpaces.Lambda = H5::DataSpace (1, dims.Lambda);
    dataSpaces.bfreqs = H5::DataSpace (1, dims.bfreqs);
    dataSpaces.ffreqs = H5::DataSpace (1, dims.ffreqs);
    dataSpaces.params = H5::DataSpace (1, dims.params);
    dataSpaces.selfenergy = H5::DataSpace (RANK_self, dims.selfenergy);
    dataSpaces.selfenergy_buffer = H5::DataSpace (RANK_irreducible, dims.irreducible); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.irreducible = H5::DataSpace (RANK_self - 1, dims.selfenergy_buffer);
    dataSpaces.irreducible_buffer = H5::DataSpace (RANK_irreducible - 1, dims.irreducible_buffer); // data space for vertex with three dimensions (three independent frequencies)
#if DIAG_CLASS >= 1
    dataSpaces.K1_a = H5::DataSpace (RANK_K1, dims.K1); // data space for vertex with three dimensions (one independent frequencies)
    dataSpaces.K1_p = H5::DataSpace (RANK_K1, dims.K1); // data space for vertex with three dimensions (one independent frequencies)
    dataSpaces.K1_t = H5::DataSpace (RANK_K1, dims.K1); // data space for vertex with three dimensions (one independent frequencies)
    dataSpaces.K1_a_buffer = H5::DataSpace (RANK_K1-1, dims.K1_buffer); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K1_p_buffer = H5::DataSpace (RANK_K1-1, dims.K1_buffer); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K1_t_buffer = H5::DataSpace (RANK_K1-1, dims.K1_buffer); // data space for vertex with three dimensions (three independent frequencies)
#endif
#if DIAG_CLASS >= 2
    dataSpaces.K2_a = H5::DataSpace (RANK_K2, dims.K2); // data space for vertex with three dimensions (two independent frequencies)
    dataSpaces.K2_p = H5::DataSpace (RANK_K2, dims.K2); // data space for vertex with three dimensions (two independent frequencies)
    dataSpaces.K2_t = H5::DataSpace (RANK_K2, dims.K2); // data space for vertex with three dimensions (two independent frequencies)
    dataSpaces.K2_a_buffer = H5::DataSpace (RANK_K2-1, dims.K2_buffer); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K2_p_buffer = H5::DataSpace (RANK_K2-1, dims.K2_buffer); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K2_t_buffer = H5::DataSpace (RANK_K2-1, dims.K2_buffer); // data space for vertex with three dimensions (three independent frequencies)
#endif
#if DIAG_CLASS >= 3
    dataSpaces.K3_a = H5::DataSpace (RANK_K3, dims.K3); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K3_p = H5::DataSpace (RANK_K3, dims.K3); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K3_t = H5::DataSpace (RANK_K3, dims.K3); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K3_a_buffer = H5::DataSpace (RANK_K3-1, dims.K3_buffer); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K3_p_buffer = H5::DataSpace (RANK_K3-1, dims.K3_buffer); // data space for vertex with three dimensions (three independent frequencies)
    dataSpaces.K3_t_buffer = H5::DataSpace (RANK_K3-1, dims.K3_buffer); // data space for vertex with three dimensions (three independent frequencies)
#endif
    return dataSpaces;
}

//class DataSpaces { // Create the data space for the dataset in file and for buffer objects
//public:
//    H5::DataSpace Lambda;
//    H5::DataSpace bfreqs;
//    H5::DataSpace ffreqs;
//    H5::DataSpace params;
//    H5::DataSpace selfenergy, selfenergy_buffer;
//    H5::DataSpace irreducible, irreducible_buffer;
//#if DIAG_CLASS >= 1
//    H5::DataSpace K1_a, K1_a_buffer;
//    H5::DataSpace K1_p, K1_p_buffer;
//    H5::DataSpace K1_t, K1_t_buffer;
//#endif
//#if DIAG_CLASS >= 2
//    H5::DataSpace K2_a, K2_a_buffer;
//    H5::DataSpace K2_p, K2_p_buffer;
//    H5::DataSpace K2_t, K2_t_buffer;
//#endif
//#if DIAG_CLASS >= 3
//    H5::DataSpace K3_a, K3_a_buffer;
//    H5::DataSpace K3_p, K3_p_buffer;
//    H5::DataSpace K3_t, K3_t_buffer;
//#endif
//
//    DataSpaces(Dims& dims) //:
////        Lambda (1, dims.Lambda),
////        bfreqs (1, dims.bfreqs),
////        ffreqs (1, dims.ffreqs),
////        params (1, dims.params),
////        selfenergy (RANK_self, dims.selfenergy),
////        selfenergy_buffer (RANK_irreducible, dims.irreducible), // data space for vertex with three dimensions (three independent frequencies)
////        irreducible (RANK_self - 1, dims.selfenergy_buffer),
////        irreducible_buffer (RANK_irreducible - 1, dims.irreducible_buffer), // data space for vertex with three dimensions (three independent frequencies)
////#if DIAG_CLASS >= 1
////        K1_a (RANK_K1, dims.K1), // data space for vertex with three dimensions (one independent frequencies)
////        K1_p (RANK_K1, dims.K1), // data space for vertex with three dimensions (one independent frequencies)
////        K1_t (RANK_K1, dims.K1), // data space for vertex with three dimensions (one independent frequencies)
////        K1_a_buffer (RANK_K1-1, dims.K1_buffer), // data space for vertex with three dimensions (three independent frequencies)
////        K1_p_buffer (RANK_K1-1, dims.K1_buffer), // data space for vertex with three dimensions (three independent frequencies)
////        K1_t_buffer (RANK_K1-1, dims.K1_buffer), // data space for vertex with three dimensions (three independent frequencies)
////#endif
////#if DIAG_CLASS >= 2
////        K2_a (RANK_K2, dims.K2), // data space for vertex with three dimensions (two independent frequencies)
////        K2_p (RANK_K2, dims.K2), // data space for vertex with three dimensions (two independent frequencies)
////        K2_t (RANK_K2, dims.K2), // data space for vertex with three dimensions (two independent frequencies)
////        K2_a_buffer (RANK_K2-1, dims.K2_buffer), // data space for vertex with three dimensions (three independent frequencies)
////        K2_p_buffer (RANK_K2-1, dims.K2_buffer), // data space for vertex with three dimensions (three independent frequencies)
////        K2_t_buffer (RANK_K2-1, dims.K2_buffer), // data space for vertex with three dimensions (three independent frequencies)
////#endif
////#if DIAG_CLASS >= 3
////        K3_a (RANK_K3, dims.K3), // data space for vertex with three dimensions (three independent frequencies)
////        K3_p (RANK_K3, dims.K3), // data space for vertex with three dimensions (three independent frequencies)
////        K3_t (RANK_K3, dims.K3), // data space for vertex with three dimensions (three independent frequencies)
////        K3_a_buffer (RANK_K3-1, dims.K3_buffer), // data space for vertex with three dimensions (three independent frequencies)
////        K3_p_buffer (RANK_K3-1, dims.K3_buffer), // data space for vertex with three dimensions (three independent frequencies)
////        K3_t_buffer (RANK_K3-1, dims.K3_buffer) // data space for vertex with three dimensions (three independent frequencies)
////#endif
//    {
//        Lambda = H5::DataSpace (1, dims.Lambda);
//        bfreqs = H5::DataSpace (1, dims.bfreqs);
//        ffreqs = H5::DataSpace (1, dims.ffreqs);
//        params = H5::DataSpace (1, dims.params);
//        selfenergy = H5::DataSpace (RANK_self, dims.selfenergy);
//        selfenergy_buffer = H5::DataSpace (RANK_irreducible, dims.irreducible); // data space for vertex with three dimensions (three independent frequencies)
//        irreducible = H5::DataSpace (RANK_self - 1, dims.selfenergy_buffer);
//        irreducible_buffer = H5::DataSpace (RANK_irreducible - 1, dims.irreducible_buffer); // data space for vertex with three dimensions (three independent frequencies)
//#if DIAG_CLASS >= 1
//        K1_a = H5::DataSpace (RANK_K1, dims.K1); // data space for vertex with three dimensions (one independent frequencies)
//        K1_p = H5::DataSpace (RANK_K1, dims.K1); // data space for vertex with three dimensions (one independent frequencies)
//        K1_t = H5::DataSpace (RANK_K1, dims.K1); // data space for vertex with three dimensions (one independent frequencies)
//        K1_a_buffer = H5::DataSpace (RANK_K1-1, dims.K1_buffer); // data space for vertex with three dimensions (three independent frequencies)
//        K1_p_buffer = H5::DataSpace (RANK_K1-1, dims.K1_buffer); // data space for vertex with three dimensions (three independent frequencies)
//        K1_t_buffer = H5::DataSpace (RANK_K1-1, dims.K1_buffer); // data space for vertex with three dimensions (three independent frequencies)
//#endif
//#if DIAG_CLASS >= 2
//        K2_a = H5::DataSpace (RANK_K2, dims.K2); // data space for vertex with three dimensions (two independent frequencies)
//        K2_p = H5::DataSpace (RANK_K2, dims.K2); // data space for vertex with three dimensions (two independent frequencies)
//        K2_t = H5::DataSpace (RANK_K2, dims.K2); // data space for vertex with three dimensions (two independent frequencies)
//        K2_a_buffer = H5::DataSpace (RANK_K2-1, dims.K2_buffer); // data space for vertex with three dimensions (three independent frequencies)
//        K2_p_buffer = H5::DataSpace (RANK_K2-1, dims.K2_buffer); // data space for vertex with three dimensions (three independent frequencies)
//        K2_t_buffer = H5::DataSpace (RANK_K2-1, dims.K2_buffer); // data space for vertex with three dimensions (three independent frequencies)
//#endif
//#if DIAG_CLASS >= 3
//        K3_a = H5::DataSpace (RANK_K3, dims.K3); // data space for vertex with three dimensions (three independent frequencies)
//        K3_p = H5::DataSpace (RANK_K3, dims.K3); // data space for vertex with three dimensions (three independent frequencies)
//        K3_t = H5::DataSpace (RANK_K3, dims.K3); // data space for vertex with three dimensions (three independent frequencies)
//        K3_a_buffer = H5::DataSpace (RANK_K3-1, dims.K3_buffer); // data space for vertex with three dimensions (three independent frequencies)
//        K3_p_buffer = H5::DataSpace (RANK_K3-1, dims.K3_buffer); // data space for vertex with three dimensions (three independent frequencies)
//        K3_t_buffer = H5::DataSpace (RANK_K3-1, dims.K3_buffer); // data space for vertex with three dimensions (three independent frequencies)
//#endif
//    }
//};

class DataSets {
public:
    // Create the datasets in file:
    // For new file: pointers
    H5::DataSet *dataset_lambda_p, *dataset_self_p, *dataset_irred_p;
    H5::DataSet *dataset_bfreqs_p, *dataset_ffreqs_p, *dataset_params_p;
#if DIAG_CLASS >= 1
    H5::DataSet *dataset_K1_a_p, *dataset_K1_p_p, *dataset_K1_t_p;
#endif
#if DIAG_CLASS >= 2
    H5::DataSet *dataset_K2_a_p, *dataset_K2_p_p, *dataset_K2_t_p;
#endif
#if DIAG_CLASS >= 3
    H5::DataSet *dataset_K3_a_p, *dataset_K3_p_p, *dataset_K3_t_p;
#endif


    // For existing file: objects
    H5::DataSet dataset_lambda, dataset_self, dataset_irred;
#if DIAG_CLASS >= 1
    H5::DataSet dataset_K1_a, dataset_K1_p, dataset_K1_t;
#endif
#if DIAG_CLASS >= 2
    H5::DataSet dataset_K2_a, dataset_K2_p, dataset_K2_t;
#endif
#if DIAG_CLASS >= 3
    H5::DataSet dataset_K3_a, dataset_K3_p, dataset_K3_t;
#endif

    DataSets(H5::H5File* file, bool file_exists,
             H5::DataSpace& dataSpaces_Lambda,
             H5::DataSpace& dataSpaces_selfenergy,
             H5::DataSpace& dataSpaces_irreducible,
             H5::DataSpace& dataSpaces_bfreqs,
             H5::DataSpace& dataSpaces_ffreqs,
             H5::DataSpace& dataSpaces_params,
             H5::DataSpace& dataSpaces_K1_a,
             H5::DataSpace& dataSpaces_K1_p,
             H5::DataSpace& dataSpaces_K1_t,
             H5::DataSpace& dataSpaces_K2_a,
             H5::DataSpace& dataSpaces_K2_p,
             H5::DataSpace& dataSpaces_K2_t,
             H5::DataSpace& dataSpaces_K3_a,
             H5::DataSpace& dataSpaces_K3_p,
             H5::DataSpace& dataSpaces_K3_t,
             H5::CompType mtype_comp, H5::CompType mtype_comp_spin, H5::DSetCreatPropList plist_vert) {
        if (!file_exists) {
            dataset_lambda_p = new H5::DataSet(
                    file->createDataSet(LAMBDA_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_Lambda));
            dataset_self_p = new H5::DataSet(file->createDataSet(SELF_LIST, mtype_comp, dataSpaces_selfenergy));
            dataset_irred_p = new H5::DataSet(
                    file->createDataSet(DATASET_irred, mtype_comp_spin, dataSpaces_irreducible, plist_vert));
            dataset_bfreqs_p = new H5::DataSet(
                    file->createDataSet(BFREQS_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs));
            dataset_ffreqs_p = new H5::DataSet(
                    file->createDataSet(FFREQS_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs));
            dataset_params_p = new H5::DataSet(
                    file->createDataSet(PARAM_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_params));
#if DIAG_CLASS >= 1
            // Create the datasets in file:
            dataset_K1_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K1_a, mtype_comp_spin, dataSpaces_K1_a, plist_vert));
            dataset_K1_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K1_p, mtype_comp_spin, dataSpaces_K1_p, plist_vert));
            dataset_K1_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K1_t, mtype_comp_spin, dataSpaces_K1_t, plist_vert));
#endif
#if DIAG_CLASS >= 2
            // Create the datasets in file:
            dataset_K2_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_a, mtype_comp_spin, dataSpaces_K2_a, plist_vert));
            dataset_K2_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_p, mtype_comp_spin, dataSpaces_K2_p, plist_vert));
            dataset_K2_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K2_t, mtype_comp_spin, dataSpaces_K2_t, plist_vert));
#endif
#if DIAG_CLASS >= 3
            // Create the datasets in file:
            dataset_K3_a_p = new H5::DataSet(
                    file->createDataSet(DATASET_K3_a, mtype_comp_spin, dataSpaces_K3_a, plist_vert));
            dataset_K3_p_p = new H5::DataSet(
                    file->createDataSet(DATASET_K3_p, mtype_comp_spin, dataSpaces_K3_p, plist_vert));
            dataset_K3_t_p = new H5::DataSet(
                    file->createDataSet(DATASET_K3_t, mtype_comp_spin, dataSpaces_K3_t, plist_vert));
#endif
        }
        else {
            dataset_lambda = file->openDataSet("lambdas");
            dataset_self = file->openDataSet("selflist");
            dataset_irred = file->openDataSet("irred");
#if DIAG_CLASS >=1
            dataset_K1_a = file->openDataSet("K1_a");
            dataset_K1_p = file->openDataSet("K1_p");
            dataset_K1_t = file->openDataSet("K1_t");
#endif
#if DIAG_CLASS >=2
            dataset_K2_a = file->openDataSet("K2_a");
            dataset_K2_p = file->openDataSet("K2_p");
            dataset_K2_t = file->openDataSet("K2_t");
#endif
#if DIAG_CLASS >=3
            dataset_K3_a = file->openDataSet("K3_a");
            dataset_K3_p = file->openDataSet("K3_p");
            dataset_K3_t = file->openDataSet("K3_t");
#endif
        }
    }
    DataSets(H5::H5File* file) {
        dataset_lambda = file->openDataSet("lambdas");
        dataset_self = file->openDataSet("selflist");
        dataset_irred = file->openDataSet("irred");
#if DIAG_CLASS >=1
        dataset_K1_a = file->openDataSet("K1_a");
        dataset_K1_p = file->openDataSet("K1_p");
        dataset_K1_t = file->openDataSet("K1_t");
#endif
#if DIAG_CLASS >=2
        dataset_K2_a = file->openDataSet("K2_a");
        dataset_K2_p = file->openDataSet("K2_p");
        dataset_K2_t = file->openDataSet("K2_t");
#endif
#if DIAG_CLASS >=3
        dataset_K3_a = file->openDataSet("K3_a");
        dataset_K3_p = file->openDataSet("K3_p");
        dataset_K3_t = file->openDataSet("K3_t");
#endif
    }
};

void save_to_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size, State<complex<double> >& state_in, vector<double>& Lambdas, bool file_exists) {
    //Try block to detect exceptions raised by any of the calls inside it
    try {
        // Prepare a buffer for writing data into the file
        Buffer buffer;
        buffer.initialize(state_in);    // copy data from state_in into the buffer

        // Turn off the auto-printing when failure occurs so that we can handle the errors appropriately
        H5::Exception::dontPrint();

        H5::H5File* file = 0;
        if (!file_exists) {
            // Create a new file using the default property lists.
            file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
        }
        else {
            // Open an existing file and dataset. Access rights: read/write
            file = new H5::H5File(FILE_NAME, H5F_ACC_RDWR);
        }


        // Create the memory datatype
        H5::CompType mtype_comp = def_mtype_comp();
        H5::CompType mtype_comp_spin = def_mtype_comp_spin();

        h5_comp_spin fillvalue_vert;
        fillvalue_vert.spin_re = 0;
        fillvalue_vert.dens_re = 0;
        fillvalue_vert.spin_im = 0;
        fillvalue_vert.dens_im = 0;
        H5::DSetCreatPropList plist_vert;
        plist_vert.setFillValue(mtype_comp_spin, &fillvalue_vert);

        h5_comp fillvalue_self;
        fillvalue_self.re = 0;
        fillvalue_self.im = 0;

        H5::DSetCreatPropList plist_self;
        plist_self.setFillValue(mtype_comp, &fillvalue_self);

        // Create the dimension arrays for objects in file and in buffer
        Dims dims(buffer, Lambda_size);

        // Create the data space for the dataset in file and for buffer objects
        // DataSpaces dataSpaces(dims);
        DataSpaces dataSpaces = initialize_DataSpaces(dims);

        // Create the data space for the dataset in file and for buffer objects
        H5::DataSpace dataSpaces_Lambda(1, dims.Lambda);
        H5::DataSpace dataSpaces_bfreqs(1, dims.bfreqs);
        H5::DataSpace dataSpaces_ffreqs(1, dims.ffreqs);
        H5::DataSpace dataSpaces_params(1, dims.params);
        H5::DataSpace dataSpaces_selfenergy(RANK_self, dims.selfenergy);
        H5::DataSpace dataSpaces_irreducible(RANK_irreducible,
                                             dims.irreducible);//data space for vertex with three dimensions (three independent frequencies)

        H5::DataSpace dataSpaces_selfenergy_buffer(RANK_self - 1, dims.selfenergy_buffer);
        H5::DataSpace dataSpaces_irreducible_buffer(RANK_irreducible - 1,
                                                    dims.irreducible_buffer);//data space for vertex with three dimensions (three independent frequencies)

#if DIAG_CLASS >= 1
        H5::DataSpace dataSpaces_K1_a(RANK_K1, dims.K1);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataSpaces_K1_p(RANK_K1, dims.K1);//data space for vertex with three dimensions (one independent frequencies)
        H5::DataSpace dataSpaces_K1_t(RANK_K1, dims.K1);//data space for vertex with three dimensions (one independent frequencies)

        H5::DataSpace dataSpaces_K1_a_buffer(RANK_K1-1, dims.K1_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K1_p_buffer(RANK_K1-1, dims.K1_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K1_t_buffer(RANK_K1-1, dims.K1_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif
#if DIAG_CLASS >= 2
        H5::DataSpace dataSpaces_K2_a(RANK_K2, dims.K2);//data space for vertex with three dimensions (twoindependent frequencies)
        H5::DataSpace dataSpaces_K2_p(RANK_K2, dims.K2);//data space for vertex with three dimensions (twoindependent frequencies)
        H5::DataSpace dataSpaces_K2_t(RANK_K2, dims.K2);//data space for vertex with three dimensions (twoindependent frequencies)

        H5::DataSpace dataSpaces_K2_a_buffer(RANK_K2-1, dims.K2_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K2_p_buffer(RANK_K2-1, dims.K2_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K2_t_buffer(RANK_K2-1, dims.K2_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif
#if DIAG_CLASS >= 3
        H5::DataSpace dataSpaces_K3_a(RANK_K3, dims.K3);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K3_p(RANK_K3, dims.K3);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K3_t(RANK_K3, dims.K3);//data space for vertex with three dimensions (three independent frequencies)

        H5::DataSpace dataSpaces_K3_a_buffer(RANK_K3-1, dims.K3_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K3_p_buffer(RANK_K3-1, dims.K3_buffer);//data space for vertex with three dimensions (three independent frequencies)
        H5::DataSpace dataSpaces_K3_t_buffer(RANK_K3-1, dims.K3_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

        //DataSets dataSets(file, file_exists, dataSpaces, mtype_comp, mtype_comp_spin, plist_vert);
        DataSets dataSets(file, file_exists,
                          dataSpaces_Lambda,
                          dataSpaces_selfenergy,
                          dataSpaces_irreducible,
                          dataSpaces_bfreqs,
                          dataSpaces_ffreqs,
                          dataSpaces_params,
                          dataSpaces_K1_a,
                          dataSpaces_K1_p,
                          dataSpaces_K1_t,
                          dataSpaces_K2_a,
                          dataSpaces_K2_p,
                          dataSpaces_K2_t,
                          dataSpaces_K3_a,
                          dataSpaces_K3_p,
                          dataSpaces_K3_t,
                          mtype_comp, mtype_comp_spin, plist_vert);


//        // Create the datasets in file:
//        H5::DataSet *dataset_lambda_p, *dataset_self_p, *dataset_irred_p;
//        H5::DataSet *dataset_bfreqs_p, *dataset_ffreqs_p, *dataset_params_p;
//#if DIAG_CLASS >= 1
//        H5::DataSet *dataset_K1_a_p, *dataset_K1_p_p, *dataset_K1_t_p;
//#endif
//#if DIAG_CLASS >= 2
//        H5::DataSet *dataset_K2_a_p, *dataset_K2_p_p, *dataset_K2_t_p;
//#endif
//#if DIAG_CLASS >= 3
//        H5::DataSet *dataset_K3_a_p, *dataset_K3_p_p, *dataset_K3_t_p;
//#endif
//
//
//
//        H5::DataSet dataset_lambda, dataset_self, dataset_irred;
//#if DIAG_CLASS >= 1
//        H5::DataSet dataset_K1_a, dataset_K1_p, dataset_K1_t;
//#endif
//#if DIAG_CLASS >= 2
//        H5::DataSet dataset_K2_a, dataset_K2_p, dataset_K2_t;
//#endif
//#if DIAG_CLASS >= 3
//        H5::DataSet dataset_K3_a, dataset_K3_p, dataset_K3_t;
//#endif
//
//
//        if (!file_exists) {
//            dataset_lambda_p = new H5::DataSet(
//                    file->createDataSet(LAMBDA_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_Lambda));
//            dataset_self_p = new H5::DataSet(file->createDataSet(SELF_LIST, mtype_comp, dataSpaces_selfenergy));
//            dataset_irred_p = new H5::DataSet(
//                    file->createDataSet(DATASET_irred, mtype_comp_spin, dataSpaces_irreducible, plist_vert));
//            dataset_bfreqs_p = new H5::DataSet(
//                    file->createDataSet(BFREQS_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_bfreqs));
//            dataset_ffreqs_p = new H5::DataSet(
//                    file->createDataSet(FFREQS_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_ffreqs));
//            dataset_params_p = new H5::DataSet(
//                    file->createDataSet(PARAM_LIST, H5::PredType::NATIVE_DOUBLE, dataSpaces_params));
//#if DIAG_CLASS >= 1
//            // Create the datasets in file:
//            dataset_K1_a_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K1_a, mtype_comp_spin, dataSpaces_K1_a, plist_vert));
//            dataset_K1_p_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K1_p, mtype_comp_spin, dataSpaces_K1_p, plist_vert));
//            dataset_K1_t_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K1_t, mtype_comp_spin, dataSpaces_K1_t, plist_vert));
//#endif
//#if DIAG_CLASS >= 2
//            // Create the datasets in file:
//            dataset_K2_a_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K2_a, mtype_comp_spin, dataSpaces_K2_a, plist_vert));
//            dataset_K2_p_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K2_p, mtype_comp_spin, dataSpaces_K2_p, plist_vert));
//            dataset_K2_t_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K2_t, mtype_comp_spin, dataSpaces_K2_t, plist_vert));
//#endif
//#if DIAG_CLASS >= 3
//            // Create the datasets in file:
//            dataset_K3_a_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K3_a, mtype_comp_spin, dataSpaces_K3_a, plist_vert));
//            dataset_K3_p_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K3_p, mtype_comp_spin, dataSpaces_K3_p, plist_vert));
//            dataset_K3_t_p = new H5::DataSet(
//                    file->createDataSet(DATASET_K3_t, mtype_comp_spin, dataSpaces_K3_t, plist_vert));
//#endif
//        }
//        else {
//            dataset_lambda = file->openDataSet("lambdas");
//            dataset_self = file->openDataSet("selflist");
//            dataset_irred = file->openDataSet("irred");
//#if DIAG_CLASS >=1
//            dataset_K1_a = file->openDataSet("K1_a");
//            dataset_K1_p = file->openDataSet("K1_p");
//            dataset_K1_t = file->openDataSet("K1_t");
//#endif
//#if DIAG_CLASS >=2
//            dataset_K2_a = file->openDataSet("K2_a");
//            dataset_K2_p = file->openDataSet("K2_p");
//            dataset_K2_t = file->openDataSet("K2_t");
//#endif
//#if DIAG_CLASS >=3
//            dataset_K3_a = file->openDataSet("K3_a");
//            dataset_K3_p = file->openDataSet("K3_p");
//            dataset_K3_t = file->openDataSet("K3_t");
//#endif
//        }



        //Select hyperslabs:

        //Select hyperslab in the file where the data should be located and after write buffer data into file.
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
            dataSets.dataset_lambda_p->write(Lambdas.data(), H5::PredType::NATIVE_DOUBLE);
            dataSets.dataset_bfreqs_p->write(bfreqs.data(), H5::PredType::NATIVE_DOUBLE);
            dataSets.dataset_ffreqs_p->write(ffreqs.data(), H5::PredType::NATIVE_DOUBLE);
            dataSets.dataset_params_p->write(parameter_list, H5::PredType::NATIVE_DOUBLE);
        }
        else
            dataSets.dataset_lambda.write(Lambdas.data(), H5::PredType::NATIVE_DOUBLE);//overwrite vector containing all values for lambda

        count[1] = buffer.self_dim;
        dataSpaces_selfenergy.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.dataset_self_p->write(buffer.selfenergy, mtype_comp, dataSpaces_selfenergy_buffer, dataSpaces_selfenergy);
        else
            dataSets.dataset_self.write(buffer.selfenergy, mtype_comp, dataSpaces_selfenergy_buffer, dataSpaces_selfenergy);

        count[1] = buffer.irred_dim;
        dataSpaces_irreducible.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        if (!file_exists)
            dataSets.dataset_irred_p->write(buffer.irreducible_class, mtype_comp_spin, dataSpaces_irreducible_buffer, dataSpaces_irreducible);
        else
            dataSets.dataset_irred.write(buffer.irreducible_class, mtype_comp_spin, dataSpaces_irreducible_buffer, dataSpaces_irreducible);

#if DIAG_CLASS >= 1
        count[1]= buffer.K1_dim;

        dataSpaces_K1_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataSpaces_K1_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataSpaces_K1_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        if (!file_exists) {
            dataSets.dataset_K1_a_p->write(buffer.K1_class_a, mtype_comp_spin, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
            dataSets.dataset_K1_p_p->write(buffer.K1_class_p, mtype_comp_spin, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
            dataSets.dataset_K1_t_p->write(buffer.K1_class_t, mtype_comp_spin, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
        }
        else {
            dataSets.dataset_K1_a.write(buffer.K1_class_a, mtype_comp_spin, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
            dataSets.dataset_K1_p.write(buffer.K1_class_p, mtype_comp_spin, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
            dataSets.dataset_K1_t.write(buffer.K1_class_t, mtype_comp_spin, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
        }
#endif

#if DIAG_CLASS >= 2
        count[1]= buffer.K2_dim;

        dataSpaces_K2_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataSpaces_K2_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataSpaces_K2_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        if (!file_exists) {
            dataSets.dataset_K2_a_p->write(buffer.K2_class_a, mtype_comp_spin, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
            dataSets.dataset_K2_p_p->write(buffer.K2_class_p, mtype_comp_spin, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
            dataSets.dataset_K2_t_p->write(buffer.K2_class_t, mtype_comp_spin, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
        }
        else {
            dataSets.dataset_K2_a.write(buffer.K2_class_a, mtype_comp_spin, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
            dataSets.dataset_K2_p.write(buffer.K2_class_p, mtype_comp_spin, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
            dataSets.dataset_K2_t.write(buffer.K2_class_t, mtype_comp_spin, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
        }
#endif

#if DIAG_CLASS >= 3
        count[1]= buffer.K3_dim;

        dataSpaces_K3_a.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataSpaces_K3_p.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);
        dataSpaces_K3_t.selectHyperslab(H5S_SELECT_SET, count,start,stride,block);

        if (!file_exists) {
            dataSets.dataset_K3_a_p->write(buffer.K3_class_a, mtype_comp_spin, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
            dataSets.dataset_K3_p_p->write(buffer.K3_class_p, mtype_comp_spin, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
            dataSets.dataset_K3_t_p->write(buffer.K3_class_t, mtype_comp_spin, dataSpaces_K3_t_buffer, dataSpaces_K3_t);
        }
        else {
            dataSets.dataset_K3_a.write(buffer.K3_class_a, mtype_comp_spin, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
            dataSets.dataset_K3_p.write(buffer.K3_class_p, mtype_comp_spin, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
            dataSets.dataset_K3_t.write(buffer.K3_class_t, mtype_comp_spin, dataSpaces_K3_t_buffer, dataSpaces_K3_t);
        }
#endif

        cout << "Successfully saved in hdf5 file: " << FILE_NAME;
        if (file_exists)
            cout << " in Lambda-layer " << Lambda_it;
        cout << endl;

        // Terminate

        dataSets.close();

        if (!file_exists) {
            dataSets.dataset_lambda_p->close();
            dataSets.dataset_bfreqs_p->close();
            dataSets.dataset_ffreqs_p->close();
            dataSets.dataset_params_p->close();
            dataSets.dataset_irred_p->close();
            dataSets.dataset_self_p->close();

#if DIAG_CLASS >= 1
            dataSets.dataset_K1_a_p->close();
            dataSets.dataset_K1_p_p->close();
            dataSets.dataset_K1_t_p->close();
#endif
#if DIAG_CLASS >= 2
            dataSets.dataset_K2_a_p->close();
            dataSets.dataset_K2_p_p->close();
            dataSets.dataset_K2_t_p->close();
#endif
#if DIAG_CLASS >= 3
            dataSets.dataset_K3_a_p->close();
            dataSets.dataset_K3_p_p->close();
            dataSets.dataset_K3_t_p->close();
#endif
        }
        else {
            dataSets.dataset_self.close();
            dataSets.dataset_irred.close();

#if DIAG_CLASS >=1
            dataSets.dataset_K1_a.close();
            dataSets.dataset_K1_p.close();
            dataSets.dataset_K1_t.close();
#endif
#if DIAG_CLASS >=2
            dataSets.dataset_K2_a.close();
            dataSets.dataset_K2_p.close();
            dataSets.dataset_K2_t.close();
#endif
#if DIAG_CLASS >=3
            dataSets.dataset_K3_a.close();
            dataSets.dataset_K3_p.close();
            dataSets.dataset_K3_t.close();
#endif
        }

        file->close();
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

// Function that writes the inital state into hdf5 format.
// The second argument denotes the iteration number to which it is written.
// The third argument denotes the total number of saved iterations in the file.
/**
 * Function that writes the inital state into hdf5 format.
 *
 * @param FILE_NAME   : File name. Creates a new file if it does not exist.
 *                      If a file with this name already exists, overwrite the existing file.
 * @param Lambda_i    : Iteration number to which the state is written.
 * @param Lambda_size : Total number of iterations to be stored in the file.
 * @param state_in    : State to be written to file.
 */
void write_hdf(const H5std_string FILE_NAME, double Lambda_i, long Lambda_size, State<complex<double> >& state_in) {
#ifdef MPI_FLAG
    if (mpi_world_rank() == 0)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
    int Lambda_it = 0;

    // List with Lambda values where only the first one is non-zero -- TODO: do we need this?
    vector<double> Lambdas (Lambda_size);
    Lambdas[0] = Lambda_i;
    for (int i = 1; i < Lambda_size; i++) {
        Lambdas[i] = 0;
    }

    save_to_hdf(FILE_NAME, Lambda_it, Lambda_size, state_in, Lambdas, false);
    }
}

// TODO: add mpi_world_rank check
////writes state of new iteration into en EXISTING Hdf5 file. The second argument denotes the iteration number to which it is written. The thrid arguemnt denotes the total number of saved iterations in the file
void add_hdf(const H5std_string FILE_NAME, int Lambda_it, long Lambda_size, State<complex<double>>& state_in, vector<double>& Lambdas) {
    //Try block to detect exceptions raised by any of the calls inside it

    if(Lambda_it < Lambda_size){
        save_to_hdf(FILE_NAME, Lambda_it, Lambda_size, state_in, Lambdas, true);
    }
    else{
        cout << "Cannot write to file " << FILE_NAME << " since Lambda layer is out of range." << endl;
    }
}

template <typename Q>
void copy_buffer_to_result(State<comp>& result, Buffer& buffer) {

    Q val; //buffer value

    for(int i=0; i<buffer.self_dim; ++i) {
        val = {buffer.selfenergy[i].re, buffer.selfenergy[i].im};
        result.selfenergy.direct_set(i, val);
    }
    for(int i=0; i<buffer.irred_dim; ++i) {
        val = {buffer.irreducible_class[i].spin_re, buffer.irreducible_class[i].spin_im};
        result.vertex.spinvertex.irred.direct_set(i, val);
        val = {buffer.irreducible_class[i].dens_re, buffer.irreducible_class[i].dens_im};
        result.vertex.densvertex.irred.direct_set(i, val);
    }
#if DIAG_CLASS >= 1
    for (int i=0; i<buffer.K1_dim; ++i) {
        val = {buffer.K1_class_a[i].spin_re, buffer.K1_class_a[i].spin_im};
        result.vertex.spinvertex.avertex.K1_direct_set(i, val);
        val = {buffer.K1_class_a[i].dens_re, buffer.K1_class_a[i].dens_im};
        result.vertex.densvertex.avertex.K1_direct_set(i, val);

        val = {buffer.K1_class_p[i].spin_re, buffer.K1_class_p[i].spin_im};
        result.vertex.spinvertex.pvertex.K1_direct_set(i, val);
        val = {buffer.K1_class_p[i].dens_re, buffer.K1_class_p[i].dens_im};
        result.vertex.densvertex.pvertex.K1_direct_set(i, val);

        val = {buffer.K1_class_t[i].spin_re, buffer.K1_class_t[i].spin_im};
        result.vertex.spinvertex.tvertex.K1_direct_set(i, val);
        val = {buffer.K1_class_t[i].dens_re, buffer.K1_class_t[i].dens_im};
        result.vertex.densvertex.tvertex.K1_direct_set(i, val);
    }
#endif
#if DIAG_CLASS >= 2
    for (int i=0; i<buffer.K2_dim; ++i) {
        val = {buffer.K2_class_a[i].spin_re, buffer.K2_class_a[i].spin_im};
        result.vertex.spinvertex.avertex.K2_direct_set(i, val);
        val = {buffer.K2_class_a[i].dens_re, buffer.K2_class_a[i].dens_im};
        result.vertex.densvertex.avertex.K2_direct_set(i, val);

        val = {buffer.K2_class_p[i].spin_re, buffer.K2_class_p[i].spin_im};
        result.vertex.spinvertex.pvertex.K2_direct_set(i, val);
        val = {buffer.K2_class_p[i].dens_re, buffer.K2_class_p[i].dens_im};
        result.vertex.densvertex.pvertex.K2_direct_set(i, val);

        val = {buffer.K2_class_t[i].spin_re, buffer.K2_class_t[i].spin_im};
        result.vertex.spinvertex.tvertex.K2_direct_set(i, val);
        val = {buffer.K2_class_t[i].dens_re, buffer.K2_class_t[i].dens_im};
        result.vertex.densvertex.tvertex.K2_direct_set(i, val);
    }
#endif
#if DIAG_CLASS >= 3
    for (int i=0; i<buffer.K3_dim; ++i) {
        val = {buffer.K3_class_a[i].spin_re, buffer.K3_class_a[i].spin_im};
        result.vertex.spinvertex.avertex.K3_direct_set(i, val);
        val = {buffer.K3_class_a[i].dens_re, buffer.K3_class_a[i].dens_im};
        result.vertex.densvertex.avertex.K3_direct_set(i, val);

        val = {buffer.K3_class_p[i].spin_re, buffer.K3_class_p[i].spin_im};
        result.vertex.spinvertex.pvertex.K3_direct_set(i, val);
        val = {buffer.K3_class_p[i].dens_re, buffer.K3_class_p[i].dens_im};
        result.vertex.densvertex.pvertex.K3_direct_set(i, val);

        val = {buffer.K3_class_t[i].spin_re, buffer.K3_class_t[i].spin_im};
        result.vertex.spinvertex.tvertex.K3_direct_set(i, val);
        val = {buffer.K3_class_t[i].dens_re, buffer.K3_class_t[i].dens_im};
        result.vertex.densvertex.tvertex.K3_direct_set(i, val);
    }
#endif
}

//function to read out an exstiting Hdf5 file. Needed to resume computation if it has been interrupted during the flow
template<typename Q>
State<complex<double>> read_hdf(const H5std_string FILE_NAME,int Lambda_it, long Lambda_size, vector<double> &Lambdas){
    State<complex<double>> result;
//    try {
        if (Lambda_it < Lambda_size) {

            H5::H5File *file = 0;
            file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);

            Buffer buffer;

            auto Lambdas_buff = new double[Lambda_size];

            //Create memory data type
            H5::CompType mtype_comp = def_mtype_comp();
            H5::CompType mtype_comp_spin = def_mtype_comp_spin();

            // Create the dimension arrays for objects in file and in buffer
            Dims dims(buffer, Lambda_size);

            // Create the data space for the dataset in file and for buffer objects
            // DataSpaces dataSpaces(dims);
//            DataSpaces dataSpaces = initialize_DataSpaces(dims);

            DataSets dataSets(file);

//            H5::DataSet dataset_lambda = file->openDataSet("lambdas");
//            H5::DataSet dataset_self = file->openDataSet("selflist");
//            H5::DataSet dataset_irred = file->openDataSet("irred");


            H5::DataSpace dataSpaces_Lambda = dataSets.dataset_lambda.getSpace();
            H5::DataSpace dataSpaces_selfenergy = dataSets.dataset_self.getSpace();
            H5::DataSpace dataSpaces_irreducible = dataSets.dataset_irred.getSpace();//data space for vertex with three dimensions (three independent frequencies)

//            dataSpaces.Lambda = dataset_lambda.getSpace();
//            dataSpaces.selfenergy = dataset_self.getSpace();
//            dataSpaces.irreducible = dataset_irred.getSpace();

            // Create the data space for buffer objects.
            H5::DataSpace dataSpaces_selfenergy_buffer(RANK_self-1, dims.selfenergy_buffer);
            H5::DataSpace dataSpaces_irreducible_buffer(RANK_irreducible-1, dims.irreducible_buffer);//data space for vertex with three dimensions (three independent frequencies)

#if DIAG_CLASS >= 1
//            H5::DataSet dataset_K1_a = file->openDataSet("K1_a");
//            H5::DataSet dataset_K1_p = file->openDataSet("K1_p");
//            H5::DataSet dataset_K1_t = file->openDataSet("K1_t");

            H5::DataSpace dataSpaces_K1_a=dataSets.dataset_K1_a.getSpace();//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K1_p=dataSets.dataset_K1_p.getSpace();//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K1_t=dataSets.dataset_K1_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
//            dataSpaces.K1_a = dataset_K1_a.getSpace();
//            dataSpaces.K1_p = dataset_K1_p.getSpace();
//            dataSpaces.K1_t = dataset_K1_t.getSpace();


            // Create the data space for buffer objects.
            H5::DataSpace dataSpaces_K1_a_buffer(RANK_K1-1, dims.K1_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K1_p_buffer(RANK_K1-1, dims.K1_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K1_t_buffer(RANK_K1-1, dims.K1_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >= 2
//            H5::DataSet dataset_K2_a = file->openDataSet("K2_a");
//            H5::DataSet dataset_K2_p = file->openDataSet("K2_p");
//            H5::DataSet dataset_K2_t = file->openDataSet("K2_t");

            H5::DataSpace dataSpaces_K2_a=dataSets.dataset_K2_a.getSpace();//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K2_p=dataSets.dataset_K2_p.getSpace();//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K2_t=dataSets.dataset_K2_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
//            dataSpaces.K2_a = dataset_K2_a.getSpace();
//            dataSpaces.K2_p = dataset_K2_p.getSpace();
//            dataSpaces.K2_t = dataset_K2_t.getSpace();

            // Create the data space for buffer objects.
            H5::DataSpace dataSpaces_K2_a_buffer(RANK_K2-1, dims.K2_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K2_p_buffer(RANK_K2-1, dims.K2_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K2_t_buffer(RANK_K2-1, dims.K2_buffer);//data space for vertex with three dimensions (three independent frequencies)
#endif

#if DIAG_CLASS >= 3
//            H5::DataSet dataset_K3_a = file->openDataSet("K3_a");
//            H5::DataSet dataset_K3_p = file->openDataSet("K3_p");
//            H5::DataSet dataset_K3_t = file->openDataSet("K3_t");

            H5::DataSpace dataSpaces_K3_a=dataSets.dataset_K3_a.getSpace();//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K3_p=dataSets.dataset_K3_p.getSpace();//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K3_t=dataSets.dataset_K3_t.getSpace();//data space for vertex with three dimensions (three independent frequencies)
//            dataSpaces.K3_a = dataset_K3_a.getSpace();
//            dataSpaces.K3_p = dataset_K3_p.getSpace();
//            dataSpaces.K3_t = dataset_K3_t.getSpace();


            // Create the data space for buffer objects.
            H5::DataSpace dataSpaces_K3_a_buffer(RANK_K3-1, dims.K3_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K3_p_buffer(RANK_K3-1, dims.K3_buffer);//data space for vertex with three dimensions (three independent frequencies)
            H5::DataSpace dataSpaces_K3_t_buffer(RANK_K3-1, dims.K3_buffer);//data space for vertex with three dimensions (three independent frequencies)

#endif


            //Select hyperslabs:



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


            count[1] = buffer.self_dim;
            dataSpaces_selfenergy.selectHyperslab(H5S_SELECT_SET, count, start, stride, block); ///
            dataSets.dataset_self.read( buffer.selfenergy, mtype_comp, dataSpaces_selfenergy_buffer, dataSpaces_selfenergy); ///

            count[1] = buffer.irred_dim;
            dataSpaces_irreducible.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSets.dataset_irred.read(buffer.irreducible_class, mtype_comp_spin, dataSpaces_irreducible_buffer, dataSpaces_irreducible);


#if DIAG_CLASS >= 1
            count[1] = buffer.K1_dim;
            dataSpaces_K1_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSpaces_K1_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSpaces_K1_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

            dataSets.dataset_K1_a.read( buffer.K1_class_a, mtype_comp_spin, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
            dataSets.dataset_K1_p.read( buffer.K1_class_p, mtype_comp_spin, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
            dataSets.dataset_K1_t.read( buffer.K1_class_t, mtype_comp_spin, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
#endif
#if DIAG_CLASS >= 2
            count[1] = buffer.K2_dim;
            dataSpaces_K2_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSpaces_K2_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSpaces_K2_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

            dataSets.dataset_K2_a.read( buffer.K2_class_a, mtype_comp_spin, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
            dataSets.dataset_K2_p.read( buffer.K2_class_p, mtype_comp_spin, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
            dataSets.dataset_K2_t.read( buffer.K2_class_t, mtype_comp_spin, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
#endif
#if DIAG_CLASS >= 3
            count[1] = buffer.K3_dim;
            dataSpaces_K3_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSpaces_K3_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSpaces_K3_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

            dataSets.dataset_K3_a.read( buffer.K3_class_a, mtype_comp_spin, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
            dataSets.dataset_K3_p.read( buffer.K3_class_p, mtype_comp_spin, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
            dataSets.dataset_K3_t.read( buffer.K3_class_t, mtype_comp_spin, dataSpaces_K3_t_buffer, dataSpaces_K3_t);

#endif


            /* // never used --> remove?
            dataset_lambda.read(Lambdas_buff, H5::PredType::NATIVE_DOUBLE, dataSpaces_Lambda);

            for (int i = 0; i < Lambda_size; i++) {
                Lambdas[i] = Lambdas_buff[i];
            }
            */


            copy_buffer_to_result<Q>(result, buffer);

#if DIAG_CLASS >=1
            dataSets.dataset_K1_a.close();
            dataSets.dataset_K1_p.close();
            dataSets.dataset_K1_t.close();
#endif
#if DIAG_CLASS >=2
            dataSets.dataset_K2_a.close();
            dataSets.dataset_K2_p.close();
            dataSets.dataset_K2_t.close();
#endif
#if DIAG_CLASS >= 3
            dataSets.dataset_K3_a.close();
            dataSets.dataset_K3_p.close();
            dataSets.dataset_K3_t.close();
#endif

            dataSets.dataset_lambda.close();
            dataSets.dataset_self.close();
            dataSets.dataset_irred.close();

            file->close();
            delete file;

            cout << "File '" << FILE_NAME << "' successfully read out" << endl;
            return result;

        } else {
            cout << "Cannot read from file " << FILE_NAME << " since Lambda layer out of range" << endl;
        }
//    }
//    // catch failure caused by the H5File operations
//    catch (H5::FileIException error) {
//        error.printErrorStack();
//        return result;
//    }
//
//    // catch failure caused by the DataSet operations
//    catch (H5::DataSetIException error) {
//        error.printErrorStack();
//        return result;
//    }
//
//    // catch failure caused by the DataSpace operations
//    catch (H5::DataSpaceIException error) {
//        error.printErrorStack();
//        return result;
//    }
}

// TODO: include i_in!
// TODO: use print, revise
void test_hdf5(H5std_string FILE_NAME, int i, State<comp>& state) {
    // test hdf5: read files and compare to original file
    int cnt = 0;
    State<comp> out = read_hdf<comp>(FILE_NAME, i, nEVO, flow_grid);
    for (int iK=0; iK<2; ++iK) {
        for (int iSE = 0; iSE < nSE; ++iSE) {
            if (state.selfenergy.val(iK, iSE, 0) != out.selfenergy.val(iK, iSE, 0)) {
                cout << "Self-energy not equal, " << iK << ", " << iSE << endl;
                cnt += 1;
            }
        }
    }

    for (int iK=0; iK<nK_K1; ++iK) {
        for (int i_in=0; i_in<n_in; ++i_in) {
#if DIAG_CLASS >= 1
            for (int iw1=0; iw1<nBOS; ++iw1) {
                if (state.vertex.densvertex.avertex.K1_val(iK, iw1, i_in) != out.vertex.densvertex.avertex.K1_val(iK, iw1, i_in)) {
                    cout << "Vertex not equal, " << iK << ", " << iw1 << endl;
                    cnt += 1;
                }
                if (state.vertex.densvertex.pvertex.K1_val(iK, iw1, i_in) != out.vertex.densvertex.pvertex.K1_val(iK, iw1, i_in)) {
                    cout << "Vertex not equal, " << iK << ", " << iw1 << endl;
                    cnt += 1;
                }
                if (state.vertex.densvertex.tvertex.K1_val(iK, iw1, i_in) != out.vertex.densvertex.tvertex.K1_val(iK, iw1, i_in)) {
                    cout << "Vertex not equal, " << iK << ", " << iw1 << endl;
                    cnt += 1;
                }
#if DIAG_CLASS >= 2
                for (int iw2=0; iw2<nFER; ++iw2) {
                    if (state.vertex.densvertex.avertex.K2_val(iK, iw1, iw2, i_in) != out.vertex.densvertex.avertex.K2_val(iK, iw1, iw2, i_in)) {
                        cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << endl;
                        cnt += 1;
                    }
                    if (state.vertex.densvertex.pvertex.K2_val(iK, iw1, iw2, i_in) != out.vertex.densvertex.pvertex.K2_val(iK, iw1, iw2, i_in)) {
                        cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << endl;
                        cnt += 1;
                    }
                    if (state.vertex.densvertex.tvertex.K2_val(iK, iw1, iw2, i_in) != out.vertex.densvertex.tvertex.K2_val(iK, iw1, iw2, i_in)) {
                        cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << endl;
                        cnt += 1;
                    }
#if DIAG_CLASS == 3
                    for (int iw3=0; iw3<nFER; ++iw3) {
                        if (state.vertex.densvertex.avertex.K3_val(iK, iw1, iw2, iw3, i_in) != out.vertex.densvertex.avertex.K3_val(iK, iw1, iw2, iw3, i_in)) {
                            cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << endl;
                            cnt += 1;
                        }
                        if (state.vertex.densvertex.pvertex.K3_val(iK, iw1, iw2, iw3, i_in) != out.vertex.densvertex.pvertex.K3_val(iK, iw1, iw2, iw3, i_in)) {
                            cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << endl;
                            cnt += 1;
                        }
                        if (state.vertex.densvertex.tvertex.K3_val(iK, iw1, iw2, iw3, i_in) != out.vertex.densvertex.tvertex.K3_val(iK, iw1, iw2, iw3, i_in)) {
                            cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << endl;
                            cnt += 1;
                        }
                    }
#endif
                }
#endif
            }
#endif
            if (state.vertex.densvertex.irred.val(iK, i_in) != out.vertex.densvertex.irred.val(iK, i_in)) {
                cout << "Vertex not equal, " << iK << endl;
                cnt += 1;
            }
        }
    }
    if (cnt == 0) {
        print("Hdf5 test successful.", true);
    }
}


