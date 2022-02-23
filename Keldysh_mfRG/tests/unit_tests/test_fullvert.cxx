#include "catch.hpp"
#include "../../correlation_functions/four_point/vertex.hpp"
#include "../../utilities/math_utils.hpp"
#include "../../utilities/hdf5_routines.hpp"

// TODO(high): Write unit tests for cross projection functionality!



TEST_CASE("Does the vectorized interpolation work for the GeneralVertex?", "vectorized interpolation") {
    if (KELDYSH and MAX_DIAG_CLASS==3) {
        /// There was once a bug with automatic type deduction that produced wrong results when compiler optimization was switched on

        // construct vertex:
        std::vector<double> K1_raw_flat(K1_expanded_config.dims_flat, 1.);
        std::vector<double> K2_raw_flat(K2_expanded_config.dims_flat, 1.);
        std::vector<double> K3_raw_flat(K3_expanded_config.dims_flat, 1.);
        multidimensional::multiarray<double,K1_expanded_config.rank> K1_raw(K1_expanded_config.dims, K1_raw_flat);
        multidimensional::multiarray<double,K2_expanded_config.rank> K2_raw(K2_expanded_config.dims, K2_raw_flat);
        multidimensional::multiarray<double,K3_expanded_config.rank> K3_raw(K3_expanded_config.dims, K3_raw_flat);

        double Lambda = 100.;
        fullvert<double> fullvertex(Lambda, true);
        GeneralVertex<double,symmetric_full> vertex_symmetricfull(fullvertex);
        vertex_symmetricfull.symmetry_expand<'a',true>(); // here only to reserve space

        for (int i = 0; i < n_in*16; i++) vertex_symmetricfull.vertices_bubbleintegrand[0].irred.direct_set(i, 1.);
        vertex_symmetricfull.vertices_bubbleintegrand[0].avertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,1));
        vertex_symmetricfull.vertices_bubbleintegrand[0].avertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,2));
        vertex_symmetricfull.vertices_bubbleintegrand[0].avertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,3));
        vertex_symmetricfull.vertices_bubbleintegrand[0].avertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,4));
        vertex_symmetricfull.vertices_bubbleintegrand[0].pvertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,5));
        vertex_symmetricfull.vertices_bubbleintegrand[0].pvertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,6));
        vertex_symmetricfull.vertices_bubbleintegrand[0].pvertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,7));
        vertex_symmetricfull.vertices_bubbleintegrand[0].pvertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,8));
        vertex_symmetricfull.vertices_bubbleintegrand[0].tvertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,9));
        vertex_symmetricfull.vertices_bubbleintegrand[0].tvertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,10));
        vertex_symmetricfull.vertices_bubbleintegrand[0].tvertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,11));
        vertex_symmetricfull.vertices_bubbleintegrand[0].tvertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,12));

        GeneralVertex<double,symmetric_r_irred> vertex_irredfull(fullvertex);
        vertex_irredfull.symmetry_expand<'a',true>(); // here only to reserve space

        for (int i = 0; i < n_in*16; i++) vertex_irredfull.vertices_bubbleintegrand[0].irred.direct_set(i, 1.);
        vertex_irredfull.vertices_bubbleintegrand[0].avertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,1));
        vertex_irredfull.vertices_bubbleintegrand[0].avertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,2));
        vertex_irredfull.vertices_bubbleintegrand[0].avertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,3));
        vertex_irredfull.vertices_bubbleintegrand[0].avertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,4));
        vertex_irredfull.vertices_bubbleintegrand[0].pvertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,5));
        vertex_irredfull.vertices_bubbleintegrand[0].pvertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,6));
        vertex_irredfull.vertices_bubbleintegrand[0].pvertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,7));
        vertex_irredfull.vertices_bubbleintegrand[0].pvertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,8));
        vertex_irredfull.vertices_bubbleintegrand[0].tvertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,9));
        vertex_irredfull.vertices_bubbleintegrand[0].tvertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,10));
        vertex_irredfull.vertices_bubbleintegrand[0].tvertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,11));
        vertex_irredfull.vertices_bubbleintegrand[0].tvertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,12));

        GeneralVertex<double,non_symmetric_diffleft> vertex_asymmfull(fullvertex,fullvertex);
        vertex_asymmfull.symmetry_expand<'a',true>(); // here only to reserve space

        for (int i = 0; i < n_in*16; i++) vertex_asymmfull.vertices_bubbleintegrand[0].irred.direct_set(i, 1.);
        vertex_asymmfull.vertices_bubbleintegrand[0].avertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,1));
        vertex_asymmfull.vertices_bubbleintegrand[0].avertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,2));
        vertex_asymmfull.vertices_bubbleintegrand[0].avertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,3));
        vertex_asymmfull.vertices_bubbleintegrand[0].avertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,4));
        vertex_asymmfull.vertices_bubbleintegrand[0].pvertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,5));
        vertex_asymmfull.vertices_bubbleintegrand[0].pvertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,6));
        vertex_asymmfull.vertices_bubbleintegrand[0].pvertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,7));
        vertex_asymmfull.vertices_bubbleintegrand[0].pvertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,8));
        vertex_asymmfull.vertices_bubbleintegrand[0].tvertex.K1_symmetry_expanded .set_vec(K1_raw*pow(2.,9));
        vertex_asymmfull.vertices_bubbleintegrand[0].tvertex.K2_symmetry_expanded .set_vec(K2_raw*pow(2.,10));
        vertex_asymmfull.vertices_bubbleintegrand[0].tvertex.K2b_symmetry_expanded.set_vec(K2_raw*pow(2.,11));
        vertex_asymmfull.vertices_bubbleintegrand[0].tvertex.K3_symmetry_expanded .set_vec(K3_raw*pow(2.,12));

        int iK = 0;
        my_index_t spin = 0;
        double w = 0.; double v = 0.; double vp = 0.;
        my_index_t internal = 0;
        char channel = 'a';
        VertexInput input(iK, spin, w, v, vp, internal, channel);

        using result_type = Eigen::Matrix<state_datatype, 4, 1>;

        SECTION("symmetric_full spin 0") {
            result_type result = vertex_symmetricfull.left_same_bare_symmetry_expanded<0,'a',result_type>(input);
            double num = pow(2,0) + pow(2,1) + pow(2,3);
            result_type expected = {num, num, num, num};
            result_type deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_symmetricfull.right_same_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,0) + pow(2,1) + pow(2,2);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_symmetricfull.left_diff_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,2) + pow(2,4) + pow(2,5) + pow(2,6) + pow(2,7) + pow(2,8) + pow(2,9) + pow(2,10) + pow(2,11) + pow(2,12);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);

            result = vertex_symmetricfull.right_diff_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,3) + pow(2,4) + pow(2,5) + pow(2,6) + pow(2,7) + pow(2,8) + pow(2,9) + pow(2,10) + pow(2,11) + pow(2,12);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);



        }

        SECTION("symmetric_full spin 1") {
            result_type result = vertex_symmetricfull.left_same_bare_symmetry_expanded<1,'a',result_type>(input);
            double num = 0.;
            result_type expected = {num, num, num, num};
            result_type deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_symmetricfull.right_same_bare_symmetry_expanded<1,'a',result_type>(input);
            num = 0.;
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_symmetricfull.left_diff_bare_symmetry_expanded<1,'a',result_type>(input);
            num = 0.;
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);

            result = vertex_symmetricfull.right_diff_bare_symmetry_expanded<1,'a',result_type>(input);
            num = 0.;
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);



        }


        SECTION("irred_full spin 0") {
            result_type result = vertex_irredfull.left_same_bare_symmetry_expanded<0,'a',result_type>(input);
            double num = 1.;
            result_type expected = {num, num, num, num};
            result_type deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_irredfull.right_same_bare_symmetry_expanded<0,'a',result_type>(input);
            num = 1.;
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_irredfull.left_diff_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,5) + pow(2,6) + pow(2,7) + pow(2,8) + pow(2,9) + pow(2,10) + pow(2,11) + pow(2,12);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);

            result = vertex_irredfull.right_diff_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,5) + pow(2,6) + pow(2,7) + pow(2,8) + pow(2,9) + pow(2,10) + pow(2,11) + pow(2,12);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);



        }

        SECTION("asym_full spin 0") {
            result_type result = vertex_asymmfull.left_same_bare_symmetry_expanded<0,'a',result_type>(input);
            double num = pow(2,0) + pow(2,1) + pow(2,3);
            result_type expected = {num, num, num, num};
            result_type deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_asymmfull.right_same_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,0) + pow(2,1) + pow(2,2);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);


            result = vertex_asymmfull.left_diff_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,2) + pow(2,4);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);

            result = vertex_asymmfull.right_diff_bare_symmetry_expanded<0,'a',result_type>(input);
            num = pow(2,3) + pow(2,4);
            expected = {num, num, num, num};
            deviation = expected - result;

            REQUIRE(deviation.lpNorm<Eigen::Infinity>() < 1e-10);



        }

    }

    }
