/**
 * Auxiliary functions for conversions of Keldysh indices
 */
#ifndef KELDYSH_MFRG_KELDYSH_SYMMETRIES_HPP
#define KELDYSH_MFRG_KELDYSH_SYMMETRIES_HPP

#include <vector>     // standard std::vector
#include <algorithm>  // for find function in isInList
#include "../utilities/util.hpp"     // printing text output
#include "../grids/frequency_grid.hpp"
#include "../multidimensional/multiarray.hpp"

/// Keldysh index parameters ///
#ifdef KELDYSH_FORMALISM
#ifndef DEBUG_SYMMETRIES
// Number of independent Keldysh components for the respective diagrammatic class
const int nK_SE = 2;
const int nK_K1 = 2;        // For channels a and t, these two are components 1 and 3 (applies for K1 and K2),
                            // for channel p components 1 and 5
const int nK_K2 = 5;        // For channels a, p and t -channel separately
const int nK_K3 = 6;        // For all channels, these 6 components are 0, 1, 3, 5, 6, 7
                            // (independent components in order of appearance)
#else // DEBUG_SYMMETRIES
const int nK_SE = 2;
const int nK_K1 = 16;       // For channels a and t, these two are components 1 and 3 (applies for K1 and K2),
                            // for channel p components 1 and 5
const int nK_K2 = 16;       // For channels a, p and t -channel separately
const int nK_K3 = 16;       // For all channels, these 6 components are 0, 1, 3, 5, 6, 7
                            // (independent components in order of appearance)
#endif // DEBUG_SYMMETRIES
#else // KELDYSH_FORMALISM
const int nK_SE = 1;
const int nK_K1 = 1;
const int nK_K2 = 1;
const int nK_K3 = 1;
#endif // KELDYSH_FORMALISM

constexpr std::array<size_t,3> dimsSE = std::array<size_t,3>({nK_SE, nFER, n_in_K1});
constexpr std::array<size_t,4> dimsK1 = std::array<size_t,4>({nK_K1, n_spin, nBOS, n_in_K1});
constexpr std::array<size_t,5> dimsK2 = std::array<size_t,5>({nK_K2, n_spin, nBOS2, nFER2, n_in_K2});
constexpr std::array<size_t,6> dimsK3 = std::array<size_t,6>({nK_K3, n_spin, nBOS3, nFER3, nFER3, n_in_K3});
constexpr size_t dimsSE_flat = getFlatSize(dimsSE);
constexpr size_t dimsK1_flat = getFlatSize(dimsK1);
constexpr size_t dimsK2_flat = getFlatSize(dimsK2);
constexpr size_t dimsK3_flat = getFlatSize(dimsK3);

#ifdef KELDYSH_FORMALISM
// Vector of indices of the non-zero Keldysh components of the bubbles
const std::vector<int> glb_non_zero_Keldysh_bubble({3,6,7,9,11,12,13,14,15});
constexpr int glb_number_of_Keldysh_components_bubble = 9; // length of the previous vector

// Vector of indices of independent components of the diagrammatic classes, density channel
const std::vector<int> non_zero_Keldysh_K1a({1,3});
const std::vector<int> non_zero_Keldysh_K2a({0,1,2,3,11});
const std::vector<int> non_zero_Keldysh_K1p({1,5});
const std::vector<int> non_zero_Keldysh_K2p({0,1,4,5,13});
const std::vector<int> non_zero_Keldysh_K1t({1,3});
const std::vector<int> non_zero_Keldysh_K2t({0,1,2,3,7});
const std::vector<int> non_zero_Keldysh_K3({0,1,3,5,6,7});

// Vector of indices whose respective Keldysh indices add up to an odd number
const std::vector<int> odd_Keldysh({1, 2, 4, 7, 8, 11, 13, 14});
#else
const std::vector<int> glb_non_zero_Keldysh_bubble {0};
constexpr int glb_number_of_Keldysh_components_bubble = 1; // length of the previous vector


// Vector of indices of independent components of the diagrammatic classes, density channel
const std::vector<int> non_zero_Keldysh_K1a({0});
const std::vector<int> non_zero_Keldysh_K2a({0});
const std::vector<int> non_zero_Keldysh_K1p({0});
const std::vector<int> non_zero_Keldysh_K2p({0});
const std::vector<int> non_zero_Keldysh_K1t({0});
const std::vector<int> non_zero_Keldysh_K2t({0});
const std::vector<int> non_zero_Keldysh_K3({0});

const std::vector<int> odd_Keldysh({1});    // trivial, never used in Matsubara
#endif // KELDYSH_FORMALISM



// Checks if a given variable val is in a list passed by reference.
// Used to check if a Keldysh index is in a list of indices that should be equal.
template <typename T>
auto isInList (T val, const std::vector<T>& list) -> bool {
    return (std::find(list.begin(), list.end(), val) != list.end());
}

// This function converts indices in the range 0...5 to the actual Keldysh index they correspond to
// Rule: {0,1,3,5,6,7} -> {0,1,2,3,4,5}
// The components 0,1,3,5,6 and 7 are the chosen reference components, numerated in the 0...15 convention
auto convertToIndepIndex(int iK) -> int;

// This function returns the values of the 4 alphas for a given index in the 0...15 set
auto alphas(int index) -> std::vector<int>;

/**
 * Function that returns, for an input i0, i2 in 0...15, the two Keldysh indices of the left [0] and right [1] vertices
 * of a bubble in a given channel.
 * @param i0      : Keldysh index of the lhs of a derivative equation for the vertex
 * @param i2      : Keldysh index of the non-zero components of the bubble propagators (takes values in a set of size 9)
 * @param channel : channel of the bubble
 * @return        : Vector of two Keldysh indices for the left [0] and right [1] vertex in a bubble
 */
auto indices_sum(int i0, int i2, const char channel) -> std::vector<int>;


template<char ch, typename Q, typename VertexValue_l>
Q load_vertex_keldyshComponents_left(multidimensional::multiarray<Q,2>& values_vertex, const VertexValue_l& value_vertex, const VertexInput& input, const int spin_idx) {
    size_t len_1 = values_vertex.length()[0];
    assert(len_1 == glb_number_of_Keldysh_components_bubble);

    if (not KELDYSH) {
        values_vertex(0, spin_idx) =  value_vertex(input);
    }
    else {

        VertexInput input_tmp = input;
        //Q v11, v12, v21, v22;
        const std::array<size_t,2> dims_vtemp = {2,2};
        // v_temp contains values for different Keldysh indices:
        // i.e. in a-channel v_temp(i,j) = v_{1',i |j,2}
        //      in p-channel v_temp(i,j) = v_{1',2'|j,i}
        //      in t-channel v_temp(i,j) = v_{i ,2'|j,2}
        multidimensional::multiarray<Q,2> v_temp(dims_vtemp);
        std::vector<int> alpha = alphas(input.iK);

        // fill v_temp:
        const std::array<size_t,4> dims_K = {2,2, 2,2};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if constexpr(ch == 'a') {input_tmp.iK = getFlatIndex<4,int,int,int,int>(alpha[0]-1, i, j, alpha[3]-1, dims_K);}
                else if constexpr(ch == 'p') {input_tmp.iK = getFlatIndex<4,int,int,int,int>(alpha[0]-1, alpha[1]-1, j, i, dims_K);}
                else if constexpr(ch == 't') {input_tmp.iK = getFlatIndex<4,int,int,int,int>(i, alpha[1]-1, j, alpha[3]-1, dims_K);}
                else assert(false);

                v_temp.at(i,j) = value_vertex(input_tmp);
            }
        }

        // fill values_vertex:
        if constexpr(ch == 'a' or ch == 't') {
            values_vertex(1, spin_idx) =  v_temp(0,0);
            values_vertex(0, spin_idx) = values_vertex(2, spin_idx) = v_temp(1,0);
            values_vertex(5, spin_idx) = values_vertex(7, spin_idx) = v_temp(0,1);
            values_vertex(3, spin_idx) = values_vertex(4, spin_idx) = values_vertex(6, spin_idx) = values_vertex(8, spin_idx) = v_temp(1,1);
        }
        else if constexpr(ch == 'p') {
            values_vertex(0, spin_idx) =  v_temp(0,0);
            values_vertex(1, spin_idx) = values_vertex(2, spin_idx) = v_temp(1,0);
            values_vertex(3, spin_idx) = values_vertex(4, spin_idx) = v_temp(0,1);
            values_vertex(5, spin_idx) = values_vertex(6, spin_idx) = values_vertex(7, spin_idx) = values_vertex(8, spin_idx) = v_temp(1,1);
        }
        else assert(false);

    }

}

template<char ch, typename Q, typename VertexValue_r>
Q load_vertex_keldyshComponents_right(multidimensional::multiarray<Q,2>& values_vertex, const VertexValue_r& value_vertex, const VertexInput& input, const int spin_idx) {
    size_t len_1 = values_vertex.length()[0];
    assert(len_1 == glb_number_of_Keldysh_components_bubble);

    if (not KELDYSH) {
        values_vertex(0, spin_idx) =  value_vertex(input);
    }
    else {
        VertexInput input_tmp = input;
        //Q v11, v12, v21, v22;
        const std::array<size_t, 2> dims_vtemp = {2, 2};
        // v_temp contains values for different Keldysh indices:
        // i.e. in a-channel v_temp(i,j) = v_{i ,2'|1,j}
        //      in p-channel v_temp(i,j) = v_{j ,i |1,2}
        //      in t-channel v_temp(i,j) = v_{1',i |1,j}
        multidimensional::multiarray<Q, 2> v_temp(dims_vtemp);
        std::vector<int> alpha = alphas(input.iK);

        // fill v_temp:
        const std::array<size_t, 4> dims_K = {2, 2, 2, 2};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if constexpr(ch == 'a') { input_tmp.iK = getFlatIndex<4,int,int,int,int>(i, alpha[1] - 1, alpha[2] - 1, j, dims_K); }
                else if constexpr(ch == 'p') { input_tmp.iK = getFlatIndex<4,int,int,int,int>(j, i, alpha[2] - 1, alpha[3] - 1, dims_K); }
                else if constexpr(ch == 't') { input_tmp.iK = getFlatIndex<4,int,int,int,int>(alpha[0] - 1, i, alpha[2] - 1, j, dims_K); }
                else
                    assert(false);

                v_temp.at(i, j) = value_vertex(input_tmp);
            }
        }

        // fill values_vertex:
        if constexpr(ch == 'a' or ch == 't') {
            values_vertex(3, spin_idx) = v_temp(0, 0);
            values_vertex(5, spin_idx) = values_vertex(6, spin_idx) = v_temp(0, 1);
            values_vertex(0, spin_idx) = values_vertex(4, spin_idx) = v_temp(1, 0);
            values_vertex(1, spin_idx) = values_vertex(2, spin_idx) = values_vertex(7, spin_idx) = values_vertex(8,
                                                                                                                 spin_idx) = v_temp(
                    1, 1);
        } else if constexpr(ch == 'p') {
            values_vertex(5, spin_idx) = v_temp(0, 0);
            values_vertex(3, spin_idx) = values_vertex(6, spin_idx) = v_temp(1, 0);
            values_vertex(1, spin_idx) = values_vertex(7, spin_idx) = v_temp(0, 1);
            values_vertex(0, spin_idx) = values_vertex(2, spin_idx) = values_vertex(4, spin_idx) = values_vertex(8,
                                                                                                                 spin_idx) = v_temp(
                    1, 1);
        } else
            assert(false);
    }
}


template<char ch, typename Q, typename VertexValue_l, typename BubbleObj, typename VertexValue_r>
Q sum_over_internal(const VertexValue_l& value_vertex_l, const BubbleObj& Pi, const VertexValue_r& value_vertex_r, const VertexInput& input_external, const double vpp) {

    VertexInput input_l = input_external, input_r = input_external;
    input_l.v2 = vpp; input_r.v1 = vpp;

    // create multiarrays to store loaded values
    // in 1st dim: Keldysh indices for 9 nonzero values of bubble
    // in 2nd dim: spin components
    multidimensional::multiarray<Q,2> values_vertex_l(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2})), values_vertex_r(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2})), values_Pi(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 1}));

    // load vertex values:
    int spin_idx = 0;
    load_vertex_keldyshComponents_left <ch,Q>(values_vertex_l, value_vertex_l, input_l, spin_idx);
    load_vertex_keldyshComponents_right<ch,Q>(values_vertex_r, value_vertex_r, input_r, spin_idx);
    // load other spin component
    if ((ch == 't' and input_external.spin == 0)
        #ifdef DEBUG_SYMMETRIES
        or (ch == 'a' and input_external.spin == 1) or (ch == 'p')
        #endif
            ){
        input_l.spin = 1 - input_external.spin;
        input_r.spin = 1 - input_external.spin;
        spin_idx = 1;

        load_vertex_keldyshComponents_left <ch,Q>(values_vertex_l, value_vertex_l, input_l, spin_idx);
        load_vertex_keldyshComponents_right<ch,Q>(values_vertex_r, value_vertex_r, input_r, spin_idx);

    }

    for (int i = 0; i < glb_number_of_Keldysh_components_bubble; i++) {
        int i2 = glb_non_zero_Keldysh_bubble[i];
        values_Pi.at(i, 0) = Pi.value(i2, input_external.w, vpp, input_external.i_in, ch);
    }

    if (false) { // set to true for debugging the loading of vertex values

        multidimensional::multiarray<Q,2> values_vertex_l_alt(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2})), values_vertex_r_alt(std::array<size_t,2>({glb_number_of_Keldysh_components_bubble, 2}));

        // load values (basic inefficient version):
        for (int i = 0; i < glb_number_of_Keldysh_components_bubble; i++) {

            VertexInput input_l_alt = input_external, input_r_alt = input_external;
            input_l_alt.v2 = vpp; input_r_alt.v1 = vpp;

            int i2 = glb_non_zero_Keldysh_bubble[i];
            std::vector<int> indices_vertices = indices_sum(input_external.iK, i2, ch);
            input_l_alt.iK = indices_vertices[0];
            input_r_alt.iK = indices_vertices[1];

            values_vertex_l_alt.at(i,0) = value_vertex_l(input_l_alt);
            values_vertex_r_alt.at(i,0) = value_vertex_r(input_r_alt);

            // load other spin component
            if ((ch == 't' and input_external.spin == 0)
                #ifdef DEBUG_SYMMETRIES
                or (ch == 'a' and input_external.spin == 1) or (ch == 'p')
                #endif
                    ){
                if(ch == 't') assert(0 == input_external.spin);
                if(ch == 'a') assert(1 == input_external.spin);
                input_l_alt.spin = 1 - input_external.spin;
                input_r_alt.spin = 1 - input_external.spin;

                values_vertex_l_alt.at(i,1) = value_vertex_l(input_l_alt);
                values_vertex_r_alt.at(i,1) = value_vertex_r(input_r_alt);

            }
        }



        multidimensional::multiarray<Q,2> values_vertex_l_diff = values_vertex_l - values_vertex_l_alt;
        multidimensional::multiarray<Q,2> values_vertex_r_diff = values_vertex_r - values_vertex_r_alt;
        auto diff_l = values_vertex_l_diff.maxabs();
        auto diff_r = values_vertex_r_diff.maxabs();
        assert(std::abs(diff_l) < 1e-10);
        assert(std::abs(diff_r) < 1e-10);


    }

    // Assemble result
    Q result = 0;
    for (int i = 0; i < glb_number_of_Keldysh_components_bubble; i++) {
        Q result_tmp;
        if ((ch == 't' and input_external.spin == 0)
            #ifdef DEBUG_SYMMETRIES
            or (ch == 'a' and input_external.spin == 1)
            #endif
                )
            result_tmp = values_vertex_l(i,0) * values_Pi(i,0) * (values_vertex_r(i,0) + values_vertex_r(i,1)) + (values_vertex_l(i,0) + values_vertex_l(i,1)) * values_Pi(i,0) * values_vertex_r(i,0);
        else
            result_tmp = values_vertex_l(i,0) * values_Pi(i,0) * values_vertex_r(i,0);

        #ifdef DEBUG_SYMMETRIES
        if (ch == 'p' and input_external.spin == 0) {
            result_tmp = (values_vertex_l(i,0) * values_Pi(i,0) * values_vertex_r(i,0) +  values_vertex_l(i,1) * values_Pi(i,0) * values_vertex_r(i,1)) * 0.5;
        }
        else if (ch == 'p' and input_external.spin == 1) {
            result_tmp = (values_vertex_l(i,0) * values_Pi(i,0) * values_vertex_r(i,1) + values_vertex_l(i,1) * values_Pi(i,0) * values_vertex_r(i,0)) * 0.5;
        }
        #endif

        result += result_tmp;

    }

    return result;
}


#endif //KELDYSH_MFRG_KELDYSH_SYMMETRIES_HPP
