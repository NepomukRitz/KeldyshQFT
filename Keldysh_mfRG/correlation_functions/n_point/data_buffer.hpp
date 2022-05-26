#ifndef KELDYSH_MFRG_DATA_BUFFER_H
#define KELDYSH_MFRG_DATA_BUFFER_H

#include "../../data_structures.hpp"
#include "../../symmetries/Keldysh_symmetries.hpp"
#include "data_container.hpp"
#include "../../interpolations/InterpolatorLinOrSloppy.hpp"
#include "../../interpolations/InterpolatorSpline1D.hpp"
#include "../../interpolations/InterpolatorSpline2D.hpp"
#include "../../interpolations/InterpolatorSpline3D.hpp"

///  TODO: Profile and improve:
/// * get_weights()
/// * get_values()
/// * handling of frequency-/index-arrays (avoid copies)

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
        if constexpr(inter == linear) base_class::frequencies.get_grid_index(freq_idx, dw_normalized, frequencies);
        else if constexpr(inter == linear_on_aux) base_class::frequencies.get_auxgrid_index(freq_idx, dw_normalized, frequencies);
        else assert(false);

        for (my_index_t i = 0; i < numberFrequencyDims; i++) {
            idx_low[pos_first_freqpoint+i] = freq_idx[i];
            assert(freq_idx[i] < 2000);
        }

        for (my_index_t j = 0; j < numSamples_half(); j++) {
            weights[j  ] = 1-dw_normalized[0];
            weights[numSamples_half()+j] = dw_normalized[0];
        }

        if constexpr (numberFrequencyDims == 2) {
            for (my_index_t j = 0; j < numSamples_half(); j++) {
                weights[2*j  ] *= 1-dw_normalized[1];
                weights[2*j+1] *= dw_normalized[1];
            }
        }
        if constexpr(numberFrequencyDims == 3) {
            for (my_index_t j = 0; j < 2; j++) {
                for (my_index_t i = 0; i < 2; i++) {
                    weights[4*j+i]   *= 1-dw_normalized[1];
                    weights[4*j+i+2] *= dw_normalized[1];
                }
            }
            for (my_index_t j = 0; j < numSamples_half(); j++) {
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
            result_type result;
            if constexpr (std::is_same_v<result_type,Q>)
            {
                if constexpr(numberFrequencyDims == 1)
                {
                    result = interpolate_nearest1D<result_type>(frequencies[0],
                                                                base_class::get_VertexFreqGrid().primary_grid,
                                                                [&](int i) -> result_type {
                                                                    indices[pos_first_freqpoint] = i;
                                                                    return base_class::val(indices);
                                                                });
                }
                else if constexpr(numberFrequencyDims == 2)
                {
                    result = interpolate_nearest2D<result_type>(frequencies[0], frequencies[1],
                                                                base_class::get_VertexFreqGrid().primary_grid,
                                                                base_class::get_VertexFreqGrid().secondary_grid,
                                                                [&](int i, int j) -> result_type {
                                                                    indices[pos_first_freqpoint] = i;
                                                                    indices[pos_first_freqpoint + 1] = j;
                                                                    return base_class::val(indices);
                                                                });
                }
                else if constexpr(numberFrequencyDims == 3)
                {
                    result = interpolate_nearest3D<result_type>(frequencies[0], frequencies[1], frequencies[2],
                                                                base_class::get_VertexFreqGrid().primary_grid,
                                                                base_class::get_VertexFreqGrid().secondary_grid,
                                                                base_class::get_VertexFreqGrid().tertiary_grid,
                                                                [&](int i, int j, int l) -> result_type {
                                                                    indices[pos_first_freqpoint] = i;
                                                                    indices[pos_first_freqpoint + 1] = j;
                                                                    indices[pos_first_freqpoint + 2] = l;
                                                                    return base_class::val(indices);
                                                                });
                }
            }
            else {
                if constexpr(numberFrequencyDims == 1)
                {
                    result << interpolate_nearest1D<Q>(frequencies[0],
                                                                base_class::get_VertexFreqGrid().primary_grid,
                                                                [&](int i) -> Q {
                                                                    indices[pos_first_freqpoint] = i;
                                                                    return base_class::val(indices);
                                                                });
                }
                else if constexpr(numberFrequencyDims == 2)
                {
                    result << interpolate_nearest2D<Q>(frequencies[0], frequencies[1],
                                                                base_class::get_VertexFreqGrid().primary_grid,
                                                                base_class::get_VertexFreqGrid().secondary_grid,
                                                                [&](int i, int j) -> Q {
                                                                    indices[pos_first_freqpoint] = i;
                                                                    indices[pos_first_freqpoint + 1] = j;
                                                                    return base_class::val(indices);
                                                                });
                }
                else if constexpr(numberFrequencyDims == 3)
                {
                    result << interpolate_nearest3D<Q>(frequencies[0], frequencies[1], frequencies[2],
                                                                base_class::get_VertexFreqGrid().primary_grid,
                                                                base_class::get_VertexFreqGrid().secondary_grid,
                                                                base_class::get_VertexFreqGrid().tertiary_grid,
                                                                [&](int i, int j, int l) -> Q {
                                                                    indices[pos_first_freqpoint] = i;
                                                                    indices[pos_first_freqpoint + 1] = j;
                                                                    indices[pos_first_freqpoint + 2] = l;
                                                                    return base_class::val(indices);
                                                                });
                }
            }
            return result;
#else
            // get weights from frequency Grid
            const weights_type weights = get_weights(frequencies, indices);
            // fetch vertex values
            const values_type values = base_class::template get_values<numberFrequencyDims, pos_first_freqpoint, vecsize, 2>(indices);
            assert(weights.allFinite());
            assert(values.allFinite());

            if constexpr(std::is_same_v<result_type, Q>) {
                const Q result = values * weights;
                assert(isfinite(result));
                return result;
            }
            else if constexpr(std::is_same_v<result_type,Eigen::Matrix<Q,result_type::RowsAtCompileTime,1>>){
                const Eigen::Matrix<Q, vecsize,1> result = values * weights;
                assert(result.allFinite());
                return result;
            }
            else {
                assert(false);
                const result_type result;
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
        //constexpr my_index_t vecsize = get_vecsize<result_type>();

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



namespace {

    template<char direction, typename dataBuffer_type, typename freqGrid_type>
    class CostResolution {
        //using freqGrid_type = bufferFrequencyGrid<k1>;

        const dataBuffer_type buffer_backup;
        bool verbose;
    public:
        dataBuffer_type buffer;
        freqGrid_type frequencies = buffer.get_VertexFreqGrid();

        explicit CostResolution(const dataBuffer_type& buffer_in, bool verbose) : buffer_backup(buffer_in), verbose(verbose), buffer(buffer_in)
                                                                              {};

        auto operator()(double wscale_test) -> double {
            // cost function for resolution  =  curvature + magnitude in the tails

            if constexpr(std::is_same_v<freqGrid_type,FrequencyGrid<eliasGrid>>) {
                //double tail_height;
                if constexpr(direction == 'b')
                {
                    frequencies.primary_grid.update_Wscale(wscale_test);
                    //tail_height = buffer.template analyze_tails<0>(1);
                }
                else {
                    frequencies.secondary_grid.update_Wscale(wscale_test);
                    frequencies.tertiary_grid.update_Wscale(wscale_test);
                    //tail_height = buffer.template analyze_tails<1>(1);
                }
            }
            else if constexpr(std::is_same_v<freqGrid_type,FrequencyGrid<angularGrid>>) {
                if constexpr(direction == 'b')
                {
                    frequencies.primary_grid.update_power(wscale_test);
                    //tail_height = buffer.template analyze_tails<0>(1);
                }
                else {
                    frequencies.secondary_grid.update_power(wscale_test);
                    frequencies.tertiary_grid.update_power(wscale_test);
                    //tail_height = buffer.template analyze_tails<1>(1);
                }
            }

            //assert(false);
            buffer.update_grid(frequencies, buffer_backup);
            //rVert.K1.analyze_tails_K1();
            const double result = buffer.get_curvature_max();

            //if (verbose) {
            //    print( "max. Curvature in buffer:\t \t", result , "\t\t with wscale = " , wscale_test , "\n");
            //}

            return result; // std::max(result, tail_height*0);
            //}
        }

        auto operator()(const std::vector<double>& section_boundaries_test) -> double {
            // make sure that 0 <= section_boundaries_test[0] <= section_boundaries_test[1]:
            if (0. > section_boundaries_test[0] or section_boundaries_test[0] > section_boundaries_test[1]) {
                return 1.e3;
            }
            else {

                if constexpr(direction == 'b') {
                    frequencies.  primary_grid.update_pos_section_boundaries(std::array<double,2>({section_boundaries_test[0],section_boundaries_test[1]}));
                }
                else {
                    frequencies.secondary_grid.update_pos_section_boundaries(std::array<double,2>({section_boundaries_test[0],section_boundaries_test[1]}));
                    frequencies. tertiary_grid.update_pos_section_boundaries(std::array<double,2>({section_boundaries_test[0],section_boundaries_test[1]}));
                }



                buffer.update_grid(frequencies, buffer_backup);
                //rVert.K1.analyze_tails_K1();
                const double result = buffer.get_curvature_max();

                //if (verbose and mpi_world_rank() == 0) {
                //    std::cout << "max. Curvature in buffer"; // << std::endl;
                //    std::cout << "\t \t" << result << "\t\t with section boundaries = " << section_boundaries_test[0] << "\t" << section_boundaries_test[1] << std::endl;
                //}

                return result;

            }
        }
    };

}




/**
 * Stores and interpolates data
 * @tparam Q                    type of data
 * @tparam k                    K_class selfenergy / k1 / k2 / k2b / k3
 * @tparam rank                 number of dimensions
 * @tparam numberFrequencyDims  number of frequency indices
 * @tparam pos_first_freqpoint  position of first frequency index
 * @tparam frequencyGrid_type   frequency grid
 * @tparam inter                interpolation method
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
        frequencies_new.guess_essential_parameters(Lambda);                     // rescale new frequency grid
        update_grid(frequencies_new, *this);
    }

    void update_grid(const frequencyGrid_type& frequencies_new, const this_class& buffer4data) {
        using buffer_type = multidimensional::multiarray<Q,rank>;
        buffer4data.initInterpolator();
        const dimensions_type dims = base_class::get_dims();
        const size_t flatsize = getFlatSize<rank>(dims);
        buffer_type data_new (dims);  // temporary data vector
        for (my_index_t iflat=0; iflat < flatsize; ++iflat) {
            index_type idx;
            std::array<my_index_t, numberFrequencyDims> i_freqs;
            frequencies_type freqs;
            getMultIndex<rank>(idx, iflat, dims);
            for (my_index_t i = 0; i < numberFrequencyDims; i++) {
                i_freqs[i] = idx[i+pos_first_freqpoint];
            }
            frequencies_new.get_freqs_w(freqs, i_freqs);
            // interpolate old values to new vector
            Q value =  buffer4data.interpolate_impl(freqs, idx);
            data_new.at(idx) = value;
        }
        base_class::set_vec(std::move(data_new)); // update vertex to new interpolated values
        base_class::set_VertexFreqGrid(frequencies_new);
        buffer4data.set_initializedInterpol(false);
    }

    void optimize_grid(bool verbose) {

        if (verbose) utils::print("---> Now Optimize grid for ", k, " in direction w:\n");

        const bool superverbose = false;
        const double epsrel = 0.01;
        const double epsabs = 0.0;

        frequencyGrid_type frequencies_new = base_class::get_VertexFreqGrid();

        /// shrink grid
        const double rel_tail_threshold = 1e-6;
        frequencies_new = base_class::shrink_freq_box(rel_tail_threshold, verbose);
        update_grid(frequencies_new, *this);
        if (verbose and mpi_world_rank() == 0) {
            std::cout << "rel.tail height in direction w: " << base_class::template analyze_tails<0>() << std::endl;
        }


        /// Optimize grid parameters
        /// in first direction
        if constexpr(std::is_same_v<typename frequencyGrid_type::grid_type1, FrequencyGrid<eliasGrid>>) {
            double a_Wscale_b = base_class::get_VertexFreqGrid().  primary_grid.W_scale / 2.;
            double m_Wscale_b = base_class::get_VertexFreqGrid().  primary_grid.W_scale;
            double b_Wscale_b = base_class::get_VertexFreqGrid().  primary_grid.W_scale * 2;
            CostResolution<'b', this_class, frequencyGrid_type> cost_b(*this, verbose);
            minimizer(cost_b, a_Wscale_b, m_Wscale_b, b_Wscale_b, 20, verbose, superverbose, epsabs, epsrel);
            frequencies_new.  primary_grid.update_Wscale(m_Wscale_b);

        }
        else if constexpr(std::is_same_v<typename frequencyGrid_type::grid_type1, FrequencyGrid<hybridGrid>>) {

            std::array<double,2> section_boundaries_b = base_class::get_VertexFreqGrid().  primary_grid.pos_section_boundaries;
            vec<double> start_params_b = {section_boundaries_b[0], section_boundaries_b[1]};

            CostResolution<'b',this_class,frequencyGrid_type> cost_b(*this, verbose);
            double ini_stepsize_b = (section_boundaries_b[0] + section_boundaries_b[1])/10.;
            vec<double> result_b = minimizer_nD(cost_b, start_params_b, ini_stepsize_b, 100, verbose, superverbose, epsabs, epsrel);
            frequencies_new.  primary_grid.update_pos_section_boundaries(std::array<double,2>({result_b[0], result_b[1]}));
        }
        else if constexpr (std::is_same_v<typename frequencyGrid_type::grid_type1, FrequencyGrid<angularGrid>>) {
            double a_power_b = base_class::get_VertexFreqGrid().  primary_grid.power / 2.;
            double m_power_b = base_class::get_VertexFreqGrid().  primary_grid.power;
            double b_power_b = base_class::get_VertexFreqGrid().  primary_grid.power * 2;
            CostResolution<'b', this_class, frequencyGrid_type> cost_b(*this, verbose);
            minimizer(cost_b, a_power_b, m_power_b, b_power_b, 20, verbose, superverbose, epsabs, epsrel);
            frequencies_new.  primary_grid.update_power(m_power_b);

        }
        update_grid(frequencies_new, *this);

        /// in other direction(s)
        if constexpr (numberFrequencyDims > 1) { /// no optimization for angular grid
            if constexpr(std::is_same_v<typename frequencyGrid_type::grid_type2, FrequencyGrid<eliasGrid>>) {
                double a_Wscale_f = base_class::get_VertexFreqGrid().secondary_grid.W_scale / 2.;
                double m_Wscale_f = base_class::get_VertexFreqGrid().secondary_grid.W_scale;
                double b_Wscale_f = base_class::get_VertexFreqGrid().secondary_grid.W_scale * 2;
                CostResolution<'f', this_class, frequencyGrid_type> cost_f(*this, verbose);
                minimizer(cost_f, a_Wscale_f, m_Wscale_f, b_Wscale_f, 20, verbose, superverbose, epsabs, epsrel);
                frequencies_new.secondary_grid.update_Wscale(m_Wscale_f);
                frequencies_new. tertiary_grid.update_Wscale(m_Wscale_f);
            }
            else if constexpr(std::is_same_v<typename frequencyGrid_type::grid_type2, FrequencyGrid<hybridGrid>>) {
                std::array<double,2> section_boundaries_f = base_class::get_VertexFreqGrid().secondary_grid.pos_section_boundaries;
                vec<double> start_params_f = {section_boundaries_f[0], section_boundaries_f[1]};

                CostResolution<'f',this_class,frequencyGrid_type> cost_f(*this, verbose);
                double ini_stepsize_f = (section_boundaries_f[0] + section_boundaries_f[1])/10.;
                vec<double> result_f = minimizer_nD(cost_f, start_params_f, ini_stepsize_f, 100, verbose, superverbose, epsabs, epsrel);
                frequencies_new.secondary_grid.update_pos_section_boundaries(std::array<double,2>({result_f[0], result_f[1]}));
                frequencies_new. tertiary_grid.update_pos_section_boundaries(std::array<double,2>({result_f[0], result_f[1]}));
            }
            else if constexpr (std::is_same_v<typename frequencyGrid_type::grid_type2, FrequencyGrid<angularGrid>>) {
                double a_power_f = base_class::get_VertexFreqGrid().secondary_grid.power / 2.;
                double m_power_f = base_class::get_VertexFreqGrid().secondary_grid.power;
                double b_power_f = base_class::get_VertexFreqGrid().secondary_grid.power * 2;
                CostResolution<'f', this_class, frequencyGrid_type> cost_f(*this, verbose);
                minimizer(cost_f, a_power_f, m_power_f, b_power_f, 20, verbose, superverbose, epsabs, epsrel);
                frequencies_new.secondary_grid.update_power(m_power_f);
                frequencies_new. tertiary_grid.update_power(m_power_f);
            }
            update_grid(frequencies_new, *this);
        }



    }

    auto operator+= (const this_class& rhs) -> this_class {base_class::data += rhs.data; return *this;}
    auto operator-= (const this_class& rhs) -> this_class {base_class::data -= rhs.data; return *this;}
    auto operator*= (const this_class& rhs) -> this_class {base_class::data *= rhs.data; return *this;}
    auto operator/= (const this_class& rhs) -> this_class {base_class::data /= rhs.data; return *this;}
    friend this_class operator+ (const this_class& lhs, const this_class& rhs) {
        this_class lhs_temp = lhs;
        lhs_temp += rhs;
        return lhs_temp;
    }
    friend this_class operator- (const this_class& lhs, const this_class& rhs) {
        this_class lhs_temp = lhs;
        lhs_temp -= rhs;
        return lhs_temp;
    }
    friend this_class operator* (const this_class& lhs, const this_class& rhs) {
        this_class lhs_temp = lhs;
        lhs_temp *= rhs;
        return lhs_temp;
    }
    friend this_class operator/ (const this_class& lhs, const this_class& rhs) {
        this_class lhs_temp = lhs;
        lhs_temp /= rhs;
        return lhs_temp;
    }


    auto operator+= (const double rhs) -> this_class {base_class::data += rhs; return *this;}
    auto operator-= (const double rhs) -> this_class {base_class::data -= rhs; return *this;}
    auto operator*= (const double rhs) -> this_class {base_class::data *= rhs; return *this;}
    auto operator/= (const double rhs) -> this_class {base_class::data /= rhs; return *this;}
    friend this_class operator+ (const this_class& lhs, const double rhs) {
        this_class lhs_temp = lhs;
        lhs_temp += rhs;
        return lhs_temp;
    }
    friend this_class operator- (const this_class& lhs, const double rhs) {
        this_class lhs_temp = lhs;
        lhs_temp -= rhs;
        return lhs_temp;
    }
    friend this_class operator* (const this_class& lhs, const double rhs) {
        this_class lhs_temp = lhs;
        lhs_temp *= rhs;
        return lhs_temp;
    }
    friend this_class operator+ (const double rhs, const this_class& lhs) {
        this_class lhs_temp = lhs;
        lhs_temp += rhs;
        return lhs_temp;
    }
    friend this_class operator* (const double rhs, const this_class& lhs) {
        this_class lhs_temp = lhs;
        lhs_temp *= rhs;
        return lhs_temp;
    }
    friend this_class operator/ (const this_class& lhs, const double rhs) {
        this_class lhs_temp = lhs;
        lhs_temp /= rhs;
        return lhs_temp;
    }
};



#endif //KELDYSH_MFRG_DATA_BUFFER_H
