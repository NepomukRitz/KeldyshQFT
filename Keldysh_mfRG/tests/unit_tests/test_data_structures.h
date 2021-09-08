#ifndef KELDYSH_MFRG_TESTING_TEST_DATA_STRUCTURES_H
#define KELDYSH_MFRG_TESTING_TEST_DATA_STRUCTURES_H

#include "../../data_structures.h"

TEST_CASE( "vector operations", "[data_structures]" ) {

    vec<comp> v (10);
    v[0] = 0. + 0.*glb_i;
    v[1] = 1.;
    v[2] = glb_i;
    v[3] = -1.;
    v[4] = -glb_i;
    v[5] = 1. + glb_i;
    v[6] = 1e-9 + glb_i;
    v[7] = 1. + 1e-9 * glb_i;
    v[8] = 1e-16 + glb_i;
    v[9] = 1. + 1e-16 * glb_i;

    SECTION( "real and imaginary parts" ) {
        vec<double> v0 (v.size());
        REQUIRE( v.real() - v.conj().real() == v0 );
        REQUIRE( v.imag() + v.conj().imag() == v0 );
    }

    SECTION( "inverse" ) {
        vec<comp> v0 = v * v.inv();

        for (int i=1; i<v0.size(); ++i) {
            REQUIRE( v0[i] == 1. );
        }
    }

    vec<comp> v1 = v;
    vec<comp> v2 = v;

    REQUIRE( v1 == v2 );

    SECTION( "multiplication" ) {
        v1 *= v;
        v2 = v * v;
        REQUIRE( v1 == v2 );
    }

    // TODO: implement more test cases
}

TEST_CASE( "multi-dimensional vector functions", "[multi-dimensional]" ) {
    const size_t dimensionality = 4;
    size_t dims[dimensionality] = {3, 5, 7, 1}; // last dimension is trivial "internal" index
    size_t dimflat = 1;
    for (int i=0; i<dimensionality; i++) {
        dimflat *= dims[i];
    }


    SECTION( "Get flat index" ) {
        int errorcount = 0;
        size_t dimstmp[dimensionality] = {3, 5, 7, 1};

        int flatidx_unrot = 0;
        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        size_t multiindex[dimensionality] = {i,j,k,l};
                        size_t test = ::getFlatIndex<dimensionality>(multiindex, dimstmp);
                        if (test != flatidx_unrot) errorcount++;
                        flatidx_unrot++;
                    }
                }
            }
        }

        REQUIRE( errorcount == 0);
    }

    SECTION( "Get rotated flat index" ) {
        int errorcount = 0;
        size_t dimstmp[dimensionality] = {5, 7, 3, 1};
        size_t permutation [dimensionality] = {2, 0, 1, 3};

        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        size_t multiindex[dimensionality] = {i,j,k,l};
                        size_t test = ::getFlatIndex<dimensionality>(multiindex, dimstmp, permutation);
                        size_t multiindex_permd[dimensionality] = {multiindex[permutation[0]], multiindex[permutation[1]], multiindex[permutation[2]], multiindex[permutation[3]]};
                        if (test != ::getFlatIndex<dimensionality>(multiindex_permd, dims)) errorcount++;
                    }
                }
            }
        }

        REQUIRE( errorcount == 0);
    }


    SECTION( "Rotate flat index" ) {
        int errorcount = 0;
        size_t dimstmp[dimensionality] = {5, 7, 3, 1};
        size_t permutation [dimensionality] = {2, 0, 1, 3};

        int flatidx_unrot = 0;
        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        size_t test = ::rotateFlatIndex<dimensionality>(flatidx_unrot, dimstmp, permutation);
                        size_t multiindex[dimensionality] = {i,j,k,l};
                        size_t multiindex_permd[dimensionality] = {multiindex[permutation[0]], multiindex[permutation[1]], multiindex[permutation[2]], multiindex[permutation[3]]};
                        if (test != ::getFlatIndex<dimensionality>(multiindex_permd, dims)) errorcount++;
                        flatidx_unrot++;
                    }
                }
            }
        }

        REQUIRE( errorcount == 0);
    }

}

namespace {
    double linfunc(size_t i, size_t j, size_t k) {
        return double(i) + 2.*double(j) + 3.*double(k);
    }
}
TEST_CASE( "Compute finite differences", "[finite_differences]") {
    const size_t dimensionality = 4;
    size_t dims[dimensionality] = {3, 5, 7, 1}; // last dimension is trivial "internal" index
    size_t dimflat = 1;
    for (int i=0; i<dimensionality; i++) {
        dimflat *= dims[i];
    }

    vec<double> linvals(dimflat);
    size_t multiindex[dimensionality] = {0,0,0,0};
    for (size_t i = 0; i<dimflat; i++) {
        ::getMultIndex<dimensionality>(multiindex, i, dims);
        linvals[i] = linfunc(multiindex[0], multiindex[1], multiindex[2]);
    }

    SECTION ( "Compute finite differences along last dimension" ) {
        vec<double> xs = {0,1,2,3,4,5,6};
        size_t permutation[dimensionality] = {1, 2, 3, 0};
        size_t dims_temp[dimensionality] = {1, 3, 5, 7};
        vec<double> dlinvals = ::get_finite_differences<double,4>(linvals, xs, dims_temp, permutation);

        int errorcount = 0;
        for (size_t i=1; i<dims[0]-1; i++) {
            for (size_t j=1; j<dims[1]-1; j++) {
                for (size_t k=1; k<dims[2]-1; k++) {
                    size_t multidx_temp[dimensionality] = {i,j,k,0};
                    if (dlinvals[::getFlatIndex<dimensionality>(multidx_temp, dims)] != 3.) errorcount ++;
                }
            }
        }


        REQUIRE( errorcount == 0);
    }


    SECTION ( "Compute finite differences along last dimension (with permutation)" ) {
        vec<double> xs = {0,1,2,3,4,5,6};
        size_t permutation[dimensionality] = {1, 2, 3, 0};
        size_t dims_temp[dimensionality] = {1, 3, 5, 7};
        vec<double> dlinvals = ::get_finite_differences(linvals, xs, dims_temp, permutation);


        int errorcount = 0;
        for (size_t i = 1; i < dims[0] - 1; i++) {
            for (size_t j = 1; j < dims[1] - 1; j++) {
                for (size_t k = 1; k < dims[2] - 1; k++) {
                    size_t multidx_temp[dimensionality] = {i, j, k, 0};
                    if (dlinvals[::getFlatIndex<dimensionality>(multidx_temp, dims)] != 3.) errorcount++;
                }
            }
        }


        REQUIRE(errorcount == 0);
    }


    SECTION ( "Compute finite differences along second dimension (with permutation)" ) {
        vec<double> xs = {0,1,2,3,4};
        size_t permutation[dimensionality] = {2, 3, 0, 1};
        size_t dims_temp[dimensionality] = {7, 1, 3, 5};
        vec<double> dlinvals = ::get_finite_differences(linvals, xs, dims_temp, permutation);


        int errorcount = 0;
        for (size_t i = 1; i < dims[0] - 1; i++) {
            for (size_t j = 1; j < dims[1] - 1; j++) {
                for (size_t k = 1; k < dims[2] - 1; k++) {
                    size_t multidx_temp[dimensionality] = {i, j, k, 0};
                    if (dlinvals[::getFlatIndex<dimensionality>(multidx_temp, dims)] != 2.) errorcount++;
                }
            }
        }


        REQUIRE(errorcount == 0);
    }


    SECTION ( "Compute finite differences along first dimension (with permutation)" ) {
        vec<double> xs = {0,1,2};
        size_t permutation[dimensionality] = {3, 0, 1, 2};
        size_t dims_temp[dimensionality] = {5, 7, 1, 3};
        vec<double> dlinvals = ::get_finite_differences(linvals, xs, dims_temp, permutation);


        int errorcount = 0;
        for (size_t i = 1; i < dims[0] - 1; i++) {
            for (size_t j = 1; j < dims[1] - 1; j++) {
                for (size_t k = 1; k < dims[2] - 1; k++) {
                    size_t multidx_temp[dimensionality] = {i, j, k, 0};
                    if (dlinvals[::getFlatIndex<dimensionality>(multidx_temp, dims)] != 1.) errorcount++;
                }
            }
        }


        REQUIRE(errorcount == 0);
    }

}

#endif //KELDYSH_MFRG_TESTING_TEST_DATA_STRUCTURES_H
