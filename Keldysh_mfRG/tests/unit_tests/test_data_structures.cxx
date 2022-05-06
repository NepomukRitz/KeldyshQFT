#include "catch.hpp"
#include "../../data_structures.hpp"
#include "../../utilities/math_utils.hpp"

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

        for (size_t i=1; i<v0.size(); ++i) {
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
    const size_t rank = 4;
    std::array<size_t,rank> dims = {3, 5, 7, 1}; // last dimension is trivial "internal" index
    size_t dimflat = 1;
    for (size_t i=0; i<rank; i++) {
        dimflat *= dims[i];
    }


    SECTION( "Get flat index" ) {
        int errorcount_flat = 0,  errorcount_multi = 0;
        std::array<size_t,rank> dimstmp = {3, 5, 7, 1};

        size_t flatidx_unrot = 0;
        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        std::array<size_t,rank> multiindex = {i,j,k,l};
                        size_t test = ::getFlatIndex<rank>(multiindex, dimstmp);
                        if (test != flatidx_unrot) errorcount_flat++;
                        size_t a,b,c,d;
                        getMultIndex<4,size_t,size_t,size_t,size_t>(a,b,c,d, test, dimstmp);
                        if (a != i or b != j or c != k or d != l) errorcount_multi++;
                        flatidx_unrot++;
                    }
                }
            }
        }

        REQUIRE( errorcount_flat == 0);
        REQUIRE( errorcount_multi == 0);
    }

    SECTION( "Get rotated flat index" ) {
        int errorcount = 0;
        std::array<size_t,rank>  dimstmp = {5, 7, 3, 1};
        std::array<size_t,rank>  permutation = {2, 0, 1, 3};

        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        std::array<size_t,rank> multiindex = {i,j,k,l};
                        size_t test = math_impl::getFlatIndex<rank>(multiindex, dimstmp, permutation);
                        std::array<size_t,rank> multiindex_permd = {multiindex[permutation[0]], multiindex[permutation[1]], multiindex[permutation[2]], multiindex[permutation[3]]};
                        if (test != ::getFlatIndex<rank>(multiindex_permd, dims)) errorcount++;
                    }
                }
            }
        }

        REQUIRE( errorcount == 0);
    }


    SECTION( "Rotate flat index" ) {
        int errorcount = 0;
        std::array<size_t,rank>  dimstmp = {5, 7, 3, 1};
        std::array<size_t,rank>  permutation = {2, 0, 1, 3};

        int flatidx_unrot = 0;
        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        size_t test = ::rotateFlatIndex<rank>(flatidx_unrot, dimstmp, permutation);
                        std::array<size_t,rank> multiindex = {i,j,k,l};
                        std::array<size_t,rank> multiindex_permd = {multiindex[permutation[0]], multiindex[permutation[1]], multiindex[permutation[2]], multiindex[permutation[3]]};
                        if (test != ::getFlatIndex<rank>(multiindex_permd, dims)) errorcount++;
                        flatidx_unrot++;
                    }
                }
            }
        }

        REQUIRE( errorcount == 0);
    }
    SECTION( "Rotate flat index cyclically" ) {
        int errorcount = 0;
        std::array<size_t,rank>  dimstmp = {5, 7, 1, 3};
        std::array<size_t,rank>  permutation  = {3, 0, 1, 2};

        size_t flatidx_unrot = 0;
        for (size_t i=0; i<dimstmp[0]; i++) {
            for (size_t j=0; j<dimstmp[1]; j++) {
                for (size_t k = 0; k < dimstmp[2]; k++) {
                    for (size_t l = 0; l < dimstmp[3]; l++) {
                        size_t test = ::rotateFlatIndex<rank>(flatidx_unrot, dims, 0);
                        std::array<size_t,rank>  multiindex = {i,j,k,l};
                        std::array<size_t,rank>  multiindex_permd = {multiindex[permutation[0]], multiindex[permutation[1]], multiindex[permutation[2]], multiindex[permutation[3]]};
                        if (test != ::getFlatIndex<rank>(multiindex_permd, dims)) errorcount++;
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
    const size_t rank = 4;
    std::array<size_t,rank> dims = {3, 5, 7, 1}; // last dimension is trivial "internal" index
    size_t dimflat = 1;
    for (size_t i=0; i<rank; i++) {
        dimflat *= dims[i];
    }

    vec<double> linvals(dimflat);
    std::array<size_t,rank>  multiindex = {0,0,0,0};
    for (size_t i = 0; i<dimflat; i++) {
        ::getMultIndex<rank>(multiindex, i, dims);
        linvals[i] = linfunc(multiindex[0], multiindex[1], multiindex[2]);
    }

    SECTION ( "Compute finite differences along last dimension" ) {
        vec<double> xs = {0,1,2,3,4,5,6};
        std::array<size_t,rank>  permutation = {1, 2, 3, 0};
        std::array<size_t,rank>  dims_temp = {1, 3, 5, 7};
        vec<double> dlinvals = ::get_finite_differences<double,4>(linvals, xs, dims_temp, permutation);

        int errorcount = 0;
        for (size_t i=0; i<dims[0]; i++) {
            for (size_t j=0; j<dims[1]; j++) {
                for (size_t k=0; k<dims[2]; k++) {
                    std::array<size_t,rank>  multidx_temp = {i,j,k,0};
                    if (dlinvals[::getFlatIndex<rank>(multidx_temp, dims)] != 3.) errorcount ++;
                }
            }
        }


        REQUIRE( errorcount == 0);
    }


    SECTION ( "Compute finite differences along last dimension (with permutation)" ) {
        vec<double> xs = {0,1,2,3,4,5,6};
        std::array<size_t,rank>  permutation = {1, 2, 3, 0};
        std::array<size_t,rank>  dims_temp = {1, 3, 5, 7};
        vec<double> dlinvals = ::get_finite_differences(linvals, xs, dims_temp, permutation);
        vec<double> dlinvals2 = ::get_finite_differences_v2(linvals, xs, dims_temp, permutation);
        vec<double> dlinvals3= ::partial_deriv(linvals, xs, dims, 2, 3);


        int errorcount = 0;
        for (size_t i = 0; i < dims[0]; i++) {
            for (size_t j = 0; j < dims[1]; j++) {
                for (size_t k = 0; k < dims[2]; k++) {
                    std::array<size_t,rank>  multidx_temp = {i, j, k, 0};
                    if (std::abs(dlinvals [::getFlatIndex<rank>(multidx_temp, dims)] - 3.) > 1e-10) errorcount++;
                    if (std::abs(dlinvals2[::getFlatIndex<rank>(multidx_temp, dims)] - 3.) > 1e-10) errorcount++;
                    if (std::abs(dlinvals3[::getFlatIndex<rank>(multidx_temp, dims)] - 3.) > 1e-10) errorcount++;
                }
            }
        }


        REQUIRE(errorcount == 0);
    }


    SECTION ( "Compute finite differences along second dimension (with permutation)" ) {
        vec<double> xs = {0,1,2,3,4};
        std::array<size_t,rank>  permutation = {2, 3, 0, 1};
        std::array<size_t,rank>  dims_temp = {7, 1, 3, 5};
        vec<double> dlinvals = ::get_finite_differences(linvals, xs, dims_temp, permutation);
        vec<double> dlinvals2= ::partial_deriv(linvals, xs, dims, 1);


        int errorcount = 0;
        for (size_t i = 0; i < dims[0] ; i++) {
            for (size_t j = 0; j < dims[1] ; j++) {
                for (size_t k = 0; k < dims[2] ; k++) {
                    std::array<size_t,rank>  multidx_temp = {i, j, k, 0};
                    if (std::abs(dlinvals [::getFlatIndex<rank>(multidx_temp, dims)] - 2.) > 1e-10) errorcount++;
                    if (std::abs(dlinvals2[::getFlatIndex<rank>(multidx_temp, dims)] - 2.) > 1e-10) errorcount++;
                }
            }
        }


        REQUIRE(errorcount == 0);
    }

/*
 *  // Current implementation of cubic interpolation scheme requires at least five grid points
    SECTION ( "Compute finite differences along first dimension (with permutation)" ) {
        vec<double> xs = {0,1,2};
        size_t permutation[rank] = {3, 0, 1, 2};
        size_t dims_temp[rank] = {5, 7, 1, 3};
        vec<double> dlinvals = ::get_finite_differences(linvals, xs, dims_temp, permutation);


        int errorcount = 0;
        for (size_t i = 1; i < dims[0] - 1; i++) {
            for (size_t j = 1; j < dims[1] - 1; j++) {
                for (size_t k = 1; k < dims[2] - 1; k++) {
                    size_t multidx_temp[rank] = {i, j, k, 0};
                    if (dlinvals[::getFlatIndex<rank>(multidx_temp, dims)] != 1.) errorcount++;
                }
            }
        }


        REQUIRE(errorcount == 0);
    }
*/
}
