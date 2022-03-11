#ifndef KELDYSH_MFRG_DATA_BUFFER_H
#define KELDYSH_MFRG_DATA_BUFFER_H

#include "../../data_structures.hpp"
#include "../../symmetries/Keldysh_symmetries.hpp"
#include "data_container.hpp"
#include "../interpolations/InterpolatorSpline1D.hpp"
/*
template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename dataContainer_type, interpolMethod inter>
class Interpolator {
    explicit Interpolator(double Lambda) {
        assert(false);
    }
};*/

template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename dataContainer_type, interpolMethod inter>
class Interpolator : public dataContainer_type {
    using base_class = dataContainer_type;
    using this_class = Interpolator;

    using frequencies_type = std::array<double, numberFrequencyDims>;

    static constexpr my_index_t numSamples() {
        return my_integer_pow<numberFrequencyDims>(my_index_t(2));
    }
    static constexpr my_index_t numSamples_half() {
        return my_integer_pow<numberFrequencyDims-1>(my_index_t(2));
    }
    template <typename result_type>
    static constexpr my_index_t get_vecsize() {
        if constexpr(std::is_same_v<result_type, Q>) return 1;
        else { return result_type::RowsAtCompileTime; }
    }
public:
    using index_type = typename base_class::index_type;
    using dimensions_type = typename base_class::dimensions_type;

    mutable bool initialized = false;
    Interpolator() : initialized(false) {};
    explicit Interpolator (double Lambda, dimensions_type dims) : base_class(Lambda, dims) {};
    void initInterpolator() const {initialized = true;};
    void set_initializedInterpol(const bool is) const {initialized = is;}


    Eigen::Matrix<double, numSamples(), 1> get_weights(const frequencies_type& frequencies, index_type& idx_low) const {
        Eigen::Matrix<double, numSamples(), 1> weights;

        std::array<my_index_t,numberFrequencyDims> freq_idx;
        std::array<double,numberFrequencyDims> dw_normalized;
        if constexpr(inter == linear) base_class::frequencies.fconv(freq_idx, dw_normalized, frequencies);
        else if constexpr(inter == linear_on_aux) base_class::frequencies.fconv_on_aux(freq_idx, dw_normalized, frequencies);
        else assert(false);

        for (my_index_t i = 0; i < numberFrequencyDims; i++) {
            idx_low[pos_first_freqpoint+i] = freq_idx[i];
            assert(freq_idx[i] < 2000);
        }

        for (int j = 0; j < numSamples_half(); j++) {
            weights[j  ] = 1-dw_normalized[0];
            weights[numSamples_half()+j] = dw_normalized[0];
        }

        if constexpr (numberFrequencyDims == 2) {
            for (int j = 0; j < numSamples_half(); j++) {
                weights[2*j  ] *= 1-dw_normalized[1];
                weights[2*j+1] *= dw_normalized[1];
            }
        }
        if constexpr(numberFrequencyDims == 3) {
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 2; i++) {
                    weights[4*j+i]   *= 1-dw_normalized[1];
                    weights[4*j+i+2] *= dw_normalized[1];
                }
            }
            for (int j = 0; j < numSamples_half(); j++) {
                weights[j*2]   *= 1-dw_normalized[2];
                weights[1+j*2] *= dw_normalized[2];
            }
        }


        return weights;

    }
    /*
    template <typename result_type, my_index_t vecsize>
    Eigen::Matrix<Q, vecsize, numSamples()> get_values(const index_type& vertex_index) const {
        Eigen::Matrix<Q, vecsize, numSamples()> result;
        if constexpr(std::is_same_v<result_type, Q>) {
            if constexpr (numberFrequencyDims == 1) {
                for (int i = 0; i < 2; i++) {
                    index_type vertex_index_tmp = vertex_index;
                    vertex_index_tmp[pos_first_freqpoint] += i;
                    result[i] = base_class::val(vertex_index_tmp);
                }
            }
            else if constexpr(numberFrequencyDims == 2){ // k2 and k2b
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        index_type vertex_index_tmp = vertex_index;
                        vertex_index_tmp[pos_first_freqpoint] += i;
                        vertex_index_tmp[pos_first_freqpoint+1] += j;
                        result[i * 2 + j] = base_class::val(vertex_index_tmp);

                    }
                }
            }
            else if constexpr (numberFrequencyDims == 3) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int l = 0; l < 2; l++) {
                            index_type vertex_index_tmp = vertex_index;
                            vertex_index_tmp[pos_first_freqpoint] += i;
                            vertex_index_tmp[pos_first_freqpoint+1] += j;
                            vertex_index_tmp[pos_first_freqpoint+2] += l;
                            result[i * 4 + j * 2 + l] = base_class::val(vertex_index_tmp);
                        }
                    }
                }
            }
            else {
                assert(false); // numberFrequencyDims > 3 not supported
            }

        }
        else { // typ == vec_left or vec_right

            if constexpr (numberFrequencyDims == 1) {
                for (int i = 0; i < 2; i++) {
                    index_type vertex_index_tmp = vertex_index;
                    vertex_index_tmp[pos_first_freqpoint] += i;
                    auto res = base_class::template val_vectorized<numberFrequencyDims,vecsize>(vertex_index_tmp);
                    result.col(i) = res;
                }
            }
            else if constexpr(numberFrequencyDims == 2) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        index_type vertex_index_tmp = vertex_index;
                        vertex_index_tmp[pos_first_freqpoint] += i;
                        vertex_index_tmp[pos_first_freqpoint+1] += j;
                        result.col(i * 2 + j) = base_class::template val_vectorized<numberFrequencyDims,vecsize>(vertex_index_tmp);

                    }
                }
            }
            else if constexpr (numberFrequencyDims == 3) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int l = 0; l < 2; l++) {
                            index_type vertex_index_tmp = vertex_index;
                            vertex_index_tmp[pos_first_freqpoint] += i;
                            vertex_index_tmp[pos_first_freqpoint+1] += j;
                            vertex_index_tmp[pos_first_freqpoint+2] += l;
                            result.col(i * 4 + j * 2 + l) = base_class::template val_vectorized<numberFrequencyDims,vecsize>(vertex_index_tmp);
                        }
                    }
                }
            }
            else {
                assert(false); // numberFrequencyDims > 3 not supported
            }

        }



        return result;
    }

    */


    template <typename result_type = Q,
            typename std::enable_if_t<(pos_first_freqpoint+numberFrequencyDims < rank) and (numberFrequencyDims <= 3), bool> = true>
    auto interpolate_impl(const frequencies_type& frequencies, index_type indices) const -> result_type {

        constexpr my_index_t vecsize = get_vecsize<result_type>();
        using weights_type = Eigen::Matrix<double, numSamples(), 1>;
        using values_type = Eigen::Matrix<Q, vecsize, numSamples()>;

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (base_class::frequencies.is_in_box(frequencies))
        {

#ifdef DENSEGRID
            Q result;
            if constexpr(k == k1)
            {
                result = interpolate_nearest1D<Q>(indices.w, base_class::get_VertexFreqGrid().b,
                                                  [&](int i) -> Q {
                                                      return base_class::val(indices.spin, i,
                                                                                             indices.iK,
                                                                                             indices.i_in);
                                                  });
            }
            else if constexpr(k==k2){
                result =  interpolate_nearest2D<Q>(indices.w, indices.v1,
                                                   base_class::get_VertexFreqGrid().b,
                                                   base_class::get_VertexFreqGrid().f,
                                                   [&](int i, int j) -> Q {
                                                       return base_class::val(indices.spin, i,
                                                                                              j, indices.iK,
                                                                                              indices.i_in);
                                                   });
            }
            else if constexpr(k==k2b){
                result =  interpolate_nearest2D<Q>(indices.w, indices.v2,
                                                   base_class::get_VertexFreqGrid().b,
                                                   base_class::get_VertexFreqGrid().f,
                                                   [&](int i, int j) -> Q {
                                                       return base_class::val(indices.spin, i,
                                                                                              j, indices.iK,
                                                                                              indices.i_in);
                                                   });
            }
            else if constexpr(k == k3)
            {
                result =  interpolate_nearest3D<Q>(indices.w, indices.v1, indices.v2,
                                                   base_class::get_VertexFreqGrid().b,
                                                   base_class::get_VertexFreqGrid().f,
                                                   base_class::get_VertexFreqGrid().f,
                                                   [&](int i, int j, int l) -> Q {
                                                       return base_class::val(indices.spin, i,
                                                                                              j, l, indices.iK,
                                                                                              indices.i_in);
                                                   });
            }

            return result;
#else
            // get weights from frequency Grid
            index_type idx_tmp = indices;
            //index_type vertex_index;
            weights_type weights = get_weights(frequencies, indices);
            // fetch vertex values
            values_type values = base_class::template get_values<numberFrequencyDims, pos_first_freqpoint, vecsize, 2>(indices);
            assert(weights.allFinite());
            assert(values.allFinite());

            if constexpr(std::is_same_v<result_type, Q>) {
                Q result = values * weights;
                assert(isfinite(result));
                return result;
            }
            else if constexpr(std::is_same_v<result_type,Eigen::Matrix<Q,result_type::RowsAtCompileTime,1>>){
                Eigen::Matrix<Q, vecsize,1> result = values * weights;
                assert(result.allFinite());
                return result;
            }
            else {
                assert(false);
                result_type result;
                return result;
            }


#endif

            //assert(isfinite(result));
            //return result;

        } else { //asymptotic value
            if constexpr (std::is_same_v<result_type,Q>) return result_type{};
            else return result_type::Zero();

        }

    };



};


/*
template<typename Q, size_t rank, my_index_t pos_first_freqpoint, typename dataContainer_type>
class Interpolator<Q, rank, 1, pos_first_freqpoint, dataContainer_type, cubic> : public Spline1D<Q,rank,pos_first_freqpoint,dataContainer_type> {
    using base_class = Spline1D<Q,rank,pos_first_freqpoint,dataContainer_type>;
    using this_class = Interpolator<Q, rank, 1, pos_first_freq, dataContainer_type, cubic>;

    using frequencies_type = std::array<double, 1>;

    static constexpr my_index_t numSamples() {
        return my_integer_pow<1>(my_index_t(2));
    }
    static constexpr my_index_t numSamples_half() {
        return my_integer_pow<1-1>(my_index_t(2));
    }
    template <typename result_type>
    static constexpr my_index_t get_vecsize() {
        if constexpr(std::is_same_v<result_type, Q>) return 1;
        else { return 4; }
    }
public:
    using index_type = typename dataContainer_type::index_type;
    using dimensions_type = typename dataContainer_type::dimensions_type;

    //mutable bool initialized = false;
    Interpolator() = default;
    explicit Interpolator (double Lambda, dimensions_type dims) : base_class(Lambda, dims) {};
    //void initInterpolator() const {initialized = true;};
    //void set_initializedInterpol(const bool is) const {initialized = is;}



    template <typename result_type = Q,
            typename std::enable_if_t<(pos_first_freqpoint+1 < rank), bool> = true>
    auto interpolate_impl(const frequencies_type& frequencies, index_type indices) const -> result_type {

        constexpr my_index_t vecsize = get_vecsize<result_type>();
        using weights_type = Eigen::Matrix<double, numSamples(), 1>;
        using values_type = Eigen::Matrix<Q, vecsize, numSamples()>;

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (base_class::frequencies.is_in_box(frequencies))
        {

#ifdef DENSEGRID
            assert(false); // Use "linear" interpolation for dense grid
#else
            result_type result = base_class::template interpolate_spline<result_type>(frequencies, indices);
            return result;

#endif

            //assert(isfinite(result));
            //return result;

        } else { //asymptotic value
            if constexpr (std::is_same_v<result_type,Q>) return result_type{};
            else return result_type::Zero();

        }

    };



};
*/
template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename dataContainer_type>
class Interpolator<Q, rank, numberFrequencyDims, pos_first_freqpoint, dataContainer_type, cubic> : public Spline<Q,rank,numberFrequencyDims,pos_first_freqpoint,dataContainer_type> {
    using base_class = Spline<Q,rank,numberFrequencyDims,pos_first_freqpoint,dataContainer_type>;
    using this_class = Interpolator<Q, rank, numberFrequencyDims, pos_first_freq, dataContainer_type, cubic>;

    using frequencies_type = std::array<double, numberFrequencyDims>;

    static constexpr my_index_t numSamples() {
        return my_integer_pow<numberFrequencyDims>(my_index_t(2));
    }
    static constexpr my_index_t numSamples_half() {
        return my_integer_pow<numberFrequencyDims-1>(my_index_t(2));
    }
    template <typename result_type>
    static constexpr my_index_t get_vecsize() {
        if constexpr(std::is_same_v<result_type, Q>) return 1;
        else { return result_type::RowsAtCompileTime; }
    }
public:
    using index_type = typename dataContainer_type::index_type;
    using dimensions_type = typename dataContainer_type::dimensions_type;

    Interpolator() = default;
    explicit Interpolator (double Lambda, dimensions_type dims) : base_class(Lambda, dims) {};

    template <typename result_type = Q,
            typename std::enable_if_t<(pos_first_freqpoint+numberFrequencyDims < rank) and (numberFrequencyDims <= 3), bool> = true>
    auto interpolate_impl(const frequencies_type& frequencies, index_type indices) const -> result_type {
#ifdef DENSEGRID
        assert(false); // Use "linear" interpolation for dense grid
#endif
        constexpr my_index_t vecsize = get_vecsize<result_type>();
        using weights_type = Eigen::Matrix<double, numSamples(), 1>;
        using values_type = Eigen::Matrix<Q, vecsize, numSamples()>;

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (base_class::frequencies.is_in_box(frequencies))
        {

            result_type result = base_class::template interpolate_spline<result_type>(frequencies, indices);
            return result;


        } else { //asymptotic value
            if constexpr (std::is_same_v<result_type,Q>) return result_type{};
            else if constexpr(std::is_same_v<result_type,Eigen::Matrix<Q,result_type::RowsAtCompileTime,1>>){return result_type::Zero();}
            else {
                assert(false);
                result_type result;
                return result;
            }
        }


    };



};



/**
 * data buffers store the data and frequency grids and interpolates data points on the grid.
 */








template <typename Q, K_class k, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type, interpolMethod inter,
        typename std::enable_if_t<(pos_first_freqpoint+numberFrequencyDims < rank) and (numberFrequencyDims <= 3), bool> = true>
class dataBuffer: public Interpolator<Q, rank, numberFrequencyDims, pos_first_freqpoint, DataContainer<Q, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type>, inter> {
    using base_class = Interpolator<Q, rank, numberFrequencyDims, pos_first_freqpoint, DataContainer<Q, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type>, inter>;
    using this_class = dataBuffer<Q, k, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type, inter>;
    using frequencies_type = std::array<double, numberFrequencyDims>;

public:
    using index_type = typename base_class::index_type;
    using dimensions_type = typename base_class::dimensions_type;

    dataBuffer() : base_class() {};
    explicit dataBuffer (double Lambda, dimensions_type dims) : base_class(Lambda, dims) {};
    //void initInterpolator() const {initialized = true;};

    template<typename result_type=Q,
            typename std::enable_if_t<(pos_first_freqpoint+numberFrequencyDims < rank) and (numberFrequencyDims <= 3), bool> = true>
    result_type interpolate(const VertexInput& input) const {
        return base_class::template interpolate_impl<result_type>(input.template get_freqs<k>(), input.template get_indices<k>());
    }

    void update_grid(double Lambda) {
        frequencyGrid_type frequencies_new = base_class::get_VertexFreqGrid();  // new frequency grid
        frequencies_new.rescale_grid(Lambda);                     // rescale new frequency grid
        update_grid(frequencies_new, *this);
    }

    void update_grid(const frequencyGrid_type& frequencies_new, const this_class& buffer4data) {
        using buffer_type = multidimensional::multiarray<Q,rank>;
        base_class::initInterpolator();
        const dimensions_type dims = base_class::get_dims();
        const size_t flatsize = getFlatSize<rank>(dims);
        buffer_type data_new (dims);  // temporary data vector
        for (int iflat=0; iflat < flatsize; ++iflat) {
            index_type idx;
            std::array<my_index_t, numberFrequencyDims> i_freqs;
            frequencies_type freqs;
            getMultIndex<rank>(idx, iflat, dims);
            for (int i = 0; i < numberFrequencyDims; i++) {
                i_freqs[i] = idx[i+pos_first_freqpoint];
            }
            frequencies_new.get_freqs_w(freqs, i_freqs);
            // interpolate old values to new vector
            data_new.at(idx) = buffer4data.interpolate_impl(freqs, idx);
        }
        base_class::set_vec(std::move(data_new)); // update vertex to new interpolated values
        base_class::set_VertexFreqGrid(frequencies_new);
        base_class::set_initializedInterpol(false);
    }

    /*
    Eigen::Matrix<double, numSamples(), 1> get_weights(const frequencies_type& frequencies, index_type& idx_low) const {
        Eigen::Matrix<double, numSamples(), 1> weights;

        std::array<my_index_t,numberFrequencyDims> freq_idx;
        std::array<double,numberFrequencyDims> dw_normalized;
        base_class::frequencies.fconv(freq_idx, dw_normalized, frequencies);
        for (my_index_t i = 0; i < numberFrequencyDims; i++) {
            idx_low[pos_first_freqpoint+i] = freq_idx[i];
        }

        for (int j = 0; j < numSamples_half(); j++) {
            weights[j  ] = 1-dw_normalized[0];
            weights[numSamples_half()+j] = dw_normalized[0];
        }

        if constexpr (numberFrequencyDims == 2) {
            for (int j = 0; j < numSamples_half(); j++) {
                weights[2*j  ] *= 1-dw_normalized[1];
                weights[2*j+1] *= dw_normalized[1];
            }
        }
        if constexpr(numberFrequencyDims == 3) {
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 2; i++) {
                    weights[4*j+i]   *= 1-dw_normalized[1];
                    weights[4*j+i+2] *= dw_normalized[1];
                }
            }
            for (int j = 0; j < numSamples_half(); j++) {
                weights[j*2]   *= 1-dw_normalized[2];
                weights[1+j*2] *= dw_normalized[2];
            }
        }


        return weights;

    }

    template <typename result_type, my_index_t vecsize>
    Eigen::Matrix<Q, vecsize, numSamples()> get_values(const index_type& vertex_index) const {
        Eigen::Matrix<Q, vecsize, numSamples()> result;
        if constexpr(std::is_same_v<result_type, Q>) {
            if constexpr (numberFrequencyDims == 1) {
                for (int i = 0; i < 2; i++) {
                    index_type vertex_index_tmp = vertex_index;
                    vertex_index_tmp[pos_first_freqpoint] += i;
                    result[i] = base_class::val(vertex_index_tmp);
                }
            }
            else if constexpr(numberFrequencyDims == 2){ // k2 and k2b
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        index_type vertex_index_tmp = vertex_index;
                        vertex_index_tmp[pos_first_freqpoint] += i;
                        vertex_index_tmp[pos_first_freqpoint+1] += j;
                        result[i * 2 + j] = base_class::val(vertex_index_tmp);

                    }
                }
            }
            else if constexpr (numberFrequencyDims == 3) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int l = 0; l < 2; l++) {
                            index_type vertex_index_tmp = vertex_index;
                            vertex_index_tmp[pos_first_freqpoint] += i;
                            vertex_index_tmp[pos_first_freqpoint+1] += j;
                            vertex_index_tmp[pos_first_freqpoint+2] += l;
                            result[i * 4 + j * 2 + l] = base_class::val(vertex_index_tmp);
                        }
                    }
                }
            }
            else {
                assert(false); // numberFrequencyDims > 3 not supported
            }

        }
        else { // typ == vec_left or vec_right

            if constexpr (numberFrequencyDims == 1) {
                for (int i = 0; i < 2; i++) {
                    index_type vertex_index_tmp = vertex_index;
                    vertex_index_tmp[pos_first_freqpoint] += i;
                    auto res = base_class::template val_vectorized<numberFrequencyDims,vecsize>(vertex_index_tmp);
                    result.col(i) = res;
                }
            }
            else if constexpr(numberFrequencyDims == 2) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        index_type vertex_index_tmp = vertex_index;
                        vertex_index_tmp[pos_first_freqpoint] += i;
                        vertex_index_tmp[pos_first_freqpoint+1] += j;
                        result.col(i * 2 + j) = base_class::template val_vectorized<numberFrequencyDims,vecsize>(vertex_index_tmp);

                    }
                }
            }
            else if constexpr (numberFrequencyDims == 3) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int l = 0; l < 2; l++) {
                            index_type vertex_index_tmp = vertex_index;
                            vertex_index_tmp[pos_first_freqpoint] += i;
                            vertex_index_tmp[pos_first_freqpoint+1] += j;
                            vertex_index_tmp[pos_first_freqpoint+2] += l;
                            result.col(i * 4 + j * 2 + l) = base_class::template val_vectorized<numberFrequencyDims,vecsize>(vertex_index_tmp);
                        }
                    }
                }
            }
            else {
                assert(false); // numberFrequencyDims > 3 not supported
            }

        }



        return result;
    }




    template <typename result_type = Q,
            typename std::enable_if_t<(pos_first_freqpoint+numberFrequencyDims < rank) and (numberFrequencyDims <= 3), bool> = true>
    auto interpolate(const frequencies_type& frequencies, index_type indices) const -> result_type {

        constexpr my_index_t vecsize = get_vecsize<result_type>();
        using weights_type = Eigen::Matrix<double, numSamples(), 1>;
        using values_type = Eigen::Matrix<Q, vecsize, numSamples()>;

        // Check if the frequency runs out of the box; if yes: return asymptotic value
        if (base_class::frequencies.is_in_box(frequencies))
        {

#ifdef DENSEGRID
            Q result;
            if constexpr(k == k1)
            {
                result = interpolate_nearest1D<Q>(indices.w, base_class::get_VertexFreqGrid().b,
                                                  [&](int i) -> Q {
                                                      return base_class::val(indices.spin, i,
                                                                                             indices.iK,
                                                                                             indices.i_in);
                                                  });
            }
            else if constexpr(k==k2){
                result =  interpolate_nearest2D<Q>(indices.w, indices.v1,
                                                   base_class::get_VertexFreqGrid().b,
                                                   base_class::get_VertexFreqGrid().f,
                                                   [&](int i, int j) -> Q {
                                                       return base_class::val(indices.spin, i,
                                                                                              j, indices.iK,
                                                                                              indices.i_in);
                                                   });
            }
            else if constexpr(k==k2b){
                result =  interpolate_nearest2D<Q>(indices.w, indices.v2,
                                                   base_class::get_VertexFreqGrid().b,
                                                   base_class::get_VertexFreqGrid().f,
                                                   [&](int i, int j) -> Q {
                                                       return base_class::val(indices.spin, i,
                                                                                              j, indices.iK,
                                                                                              indices.i_in);
                                                   });
            }
            else if constexpr(k == k3)
            {
                result =  interpolate_nearest3D<Q>(indices.w, indices.v1, indices.v2,
                                                   base_class::get_VertexFreqGrid().b,
                                                   base_class::get_VertexFreqGrid().f,
                                                   base_class::get_VertexFreqGrid().f,
                                                   [&](int i, int j, int l) -> Q {
                                                       return base_class::val(indices.spin, i,
                                                                                              j, l, indices.iK,
                                                                                              indices.i_in);
                                                   });
            }

            return result;
#else
            // get weights from frequency Grid

            //index_type vertex_index;
            weights_type weights = get_weights(frequencies, indices);
            // fetch vertex values
            values_type values = get_values<result_type, vecsize>(indices);
            assert(weights.allFinite());
            assert(values.allFinite());

            if constexpr(std::is_same_v<result_type, Q>) {
                Q result = values * weights;
                assert(isfinite(result));
                return result;
            }
            else {
                Eigen::Matrix<Q, vecsize,1> result = values * weights;
                assert(result.allFinite());
                return result;
            }


#endif

            //assert(isfinite(result));
            //return result;

        } else { //asymptotic value
            if constexpr (std::is_same_v<result_type,Q>) return result_type{};
            else return result_type::Zero();

        }

    };
*/
    auto operator+= (const this_class& rhs) -> this_class {base_class::data += rhs.data; return *this;}
    auto operator-= (const this_class& rhs) -> this_class {base_class::data -= rhs.data; return *this;}
    auto operator*= (const this_class& rhs) -> this_class {base_class::data *= rhs.data; return *this;}
    auto operator/= (const this_class& rhs) -> this_class {base_class::data /= rhs.data; return *this;}
    friend this_class& operator+ (this_class& lhs, const this_class& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend this_class& operator- (this_class& lhs, const this_class& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend this_class& operator* (this_class& lhs, const this_class& rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend this_class& operator/ (this_class& lhs, const this_class& rhs) {
        lhs /= rhs;
        return lhs;
    }


    auto operator+= (const double rhs) -> this_class {base_class::data += rhs; return *this;}
    auto operator-= (const double rhs) -> this_class {base_class::data -= rhs; return *this;}
    auto operator*= (const double rhs) -> this_class {base_class::data *= rhs; return *this;}
    auto operator/= (const double rhs) -> this_class {base_class::data /= rhs; return *this;}
    friend this_class& operator+ (this_class& lhs, const double rhs) {
        lhs += rhs;
        return lhs;
    }
    friend this_class& operator- (this_class& lhs, const double rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend this_class& operator* (this_class& lhs, const double rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend this_class& operator+ (const double rhs, this_class& lhs) {
        lhs += rhs;
        return lhs;
    }
    friend this_class& operator* (const double rhs, this_class& lhs) {
        lhs *= rhs;
        return lhs;
    }
    friend this_class& operator/ (this_class& lhs, const double rhs) {
        lhs /= rhs;
        return lhs;
    }
};



#endif //KELDYSH_MFRG_DATA_BUFFER_H
