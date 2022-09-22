#ifndef FPP_MFRG_VERTEX_DATA_CONTAINER_H
#define FPP_MFRG_VERTEX_DATA_CONTAINER_H

/**
 * This header contains a class which is responsible for saving and retrieving vertex data
 */

#include "../../utilities/template_utils.hpp"
#include "../../utilities/math_utils.hpp"
#include "../../symmetries/Keldysh_symmetries.hpp"
#include "../../data_structures.hpp"          // real/complex vector classes
#include "../../parameters/master_parameters.hpp"               // system parameters (lengths of vectors etc.)
#include "../../parameters/frequency_parameters.hpp"
#include "../../grids/frequency_grid.hpp"            // functionality for the internal structure of the Hubbard model
#include "H5Cpp.h"

/// TODO: treat VertexFrequencyGrid as member of vertexContainerBase via template?

template <typename Q> class rvert; // forward declaration of rvert
//template <typename Q> class fullvert; // forward declaration of fullvert
template <typename Q, bool differentiated> class State; // forward declaration of State
template<typename Q, std::size_t depth, typename H5object>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const multidimensional::multiarray<Q, depth>& data, const bool data_set_exists);
template<typename Q, typename H5object, int nrows, int ncols>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const Eigen::Matrix<Q,nrows, ncols>& data, const bool data_set_exists);
template <typename Q, typename H5object>
void write_to_hdf(H5object& group, const H5std_string& dataset_name, const std::vector<Q>& data, const bool data_set_exists);
//template <typename Q, vertexType symm_type> class GeneralVertex;
////template <typename Q>class symmetric_full;
//template <typename Q>using Vertex = GeneralVertex<Q, symmetric_full>;
class Buffer;
//namespace n_point_fc {
/**
 * Offers basic functionality that is identical for all K_classes
 * @tparam Q        data type of vertex data
 * @tparam rank     rank of data tensor
 */
    template<typename Q, std::size_t rank>
    class dataContainerBase {
        friend class State<Q,false>;
        friend class State<Q,true>;
        friend State<state_datatype,false> read_state_from_hdf(const H5std_string& filename, const unsigned int Lambda_it);


    protected:
        using buffer_type = multidimensional::multiarray<Q, rank>;
        using index_type = typename buffer_type::index_type;
        using dimensions_type = typename buffer_type::dimensions_type;

        buffer_type data;

    public:
        /// constructor:
        dataContainerBase() = default;

        explicit dataContainerBase(const dimensions_type dims) : data(dims) {};

        template<typename... Types,
                typename std::enable_if_t<
                        (sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
        explicit
        dataContainerBase(const Types &... dims) : dataContainerBase(index_type({static_cast<size_t>(dims)...})) {};

        explicit dataContainerBase(const buffer_type &data_in) : buffer_type(data_in) {};

        /// Access data via flattened index.
        Q acc(const size_t flatIndex) const {
            assert(flatIndex < data.size());
            return data.flat_at(flatIndex);
        }

        void direct_set(const size_t flatIndex, Q value) {
            assert(flatIndex < data.size());
            data.flat_at(flatIndex) = value;
        }

        /// Returns value for a multiIndex
        template<typename... Types,
                typename std::enable_if_t<
                        (sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
        Q val(const Types &... i) const { return data(i...); }

        Q val(const index_type &idx) const { return data.at(idx); }

        /// Returns reference to a value for a multiIndex
        template<typename... Types,
                typename std::enable_if_t<
                        (sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true>
        const Q &at(const Types &... i) const { return data.at(i...); }

        template<std::size_t freqrank, std::size_t vecsize, typename... Types,
                typename std::enable_if_t<(sizeof...(Types) == freqrank + pos_first_freq + 1) and
                                          (are_all_integral<size_t, Types...>::value), bool> = true>
        auto val_vectorized(const Types &... i) const -> Eigen::Matrix<Q, vecsize, 1> {
            return data.template at_vectorized<vecsize>(i...);
        }
        template<std::size_t freqrank, std::size_t vecsize>
        auto val_vectorized(const index_type &idx) const -> Eigen::Matrix<Q, vecsize, 1> {
            return data.template at_vectorized<vecsize>(idx);
        }

        /// Sets a value at a multiIndex
        template<typename... Types, typename std::enable_if_t<
                (sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true
        >
        void setvert(const Q value, const Types &... i) { data.at(i...) = value; }

        void setvert(const Q value, const index_type &idx) { data.at(idx) = value; }

        template<std::size_t vecsize>
        void setvert_vectorized(const Eigen::Matrix<Q, vecsize, 1>& value, const index_type &idx) {
            data.template set_vectorized<vecsize>(value, idx);
        }

        auto get_dims() const { return data.length(); }

        /// Returns the buffer "data" containing the data
        const buffer_type &get_vec() const { return data; }

        /// Sets the the buffer "data"
        template<typename container,
                std::enable_if_t<std::is_same_v <
                                 typename container::value_type, Q> && !std::is_same_v <container, buffer_type>, bool> = true
        >

        void set_vec(const container &data_in) {
            assert(data.size() == data_in.size());
            data = buffer_type(data.length(), data_in);
        }

        void set_vec(const buffer_type &data_in) {
            assert(data.is_same_length(data_in));
            data = data_in;
        }

        void set_vec(const buffer_type &&data_in) {
            assert(data.is_same_length(data_in));
            data = data_in;
        }

        /// Adds a vector to the data
        template<typename container,
                std::enable_if_t<std::is_same_v <
                                 typename container::value_type, Q> && !std::is_same_v <container, buffer_type>, bool> = true
        >

        void add_vec(const container &summand) {
            assert(data.size() == summand.size()); /// Check that summand has the right length
            data += buffer_type(data.length(), summand);
        }

        void add_vec(const buffer_type &summand) {
            assert(data.is_same_length(summand)); /// Check that summand has the right length
            data += summand;
        }

        constexpr auto eigen_segment(const index_type &start, const index_type &end) {
            return data.eigen_segment(start, end);
        }

        constexpr auto eigen_segment(const index_type &start, const index_type &end) const {
            return data.eigen_segment(start, end);
        }

        template <my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, my_index_t vecsize, my_index_t sample_size>
        Eigen::Matrix<Q, vecsize, my_integer_pow<numberFrequencyDims>(sample_size)> get_values(const index_type& index) const {
            return data.template get_values<numberFrequencyDims, pos_first_freqpoint, vecsize, sample_size>(index);
        }

    };


    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type>
    class DataContainer : public dataContainerBase<Q, rank> {
        friend void test_PT4(double Lambda, bool write_flag);

        template<typename T>
        friend void test_PT_state(std::string outputFileName, double Lambda, bool write_flag);

        template<typename T>
        friend void result_set_frequency_grids(State<T,false> &result, Buffer &buffer);

        template<typename T>
        friend void check_FDTs(const State<T,false> &state, bool verbose);

        friend State<state_datatype,false> read_state_from_hdf(const H5std_string &filename, const int Lambda_it);

    protected:
        using base_class = dataContainerBase<Q, rank>;
    public:
        using index_type = typename base_class::index_type;
        using dimensions_type = typename base_class::dimensions_type;
        using buffer_type = typename base_class::buffer_type;

        frequencyGrid_type frequencies;    // frequency grid
        DataContainer() = default;

        explicit DataContainer(double Lambda, dimensions_type dims, const fRG_config& config) : base_class(dims), frequencies(Lambda, config) {};

        /// Functions for getting and setting the frequency grid and its members
        auto get_VertexFreqGrid() const -> const frequencyGrid_type &;

        void set_VertexFreqGrid(const frequencyGrid_type frequencyGrid);


        template<unsigned int idim,
                typename std::enable_if_t<(idim < numberFrequencyDims), bool> = true
                >
        double analyze_tails(const int index = 0) const {
            double maxabs_total = base_class::data.max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint+idim);

            const size_t number_of_frequency_points_along_w = base_class::get_dims()[pos_first_freqpoint+idim];
            return maxabs_along_w[number_of_frequency_points_along_w - 1 - index] / maxabs_total;
        };


        buffer_type get_deriv_x () const {
            buffer_type result = ::partial_deriv<Q,rank>(base_class::data, frequencies.  primary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(),  pos_first_freqpoint);
            return result;
        };
        buffer_type get_deriv_y () const {
            buffer_type result = ::partial_deriv<Q,rank>(base_class::data, frequencies.secondary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(),  pos_first_freqpoint+1);
            return result;
        };
        buffer_type get_deriv_z () const {
            buffer_type result = ::partial_deriv<Q,rank>(base_class::data, frequencies. tertiary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(),  pos_first_freqpoint+2);
            return result;
        };

        buffer_type get_deriv_xx() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.  primary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint  );
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.  primary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint  );
            return result;
        }
        buffer_type get_deriv_yy() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.secondary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.secondary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+1);
            return result;
        }
        buffer_type get_deriv_zz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.tertiary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+2);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.tertiary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+2);
            return result;
        }
        buffer_type get_deriv_xy() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.secondary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.  primary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint  );
            return result;
        }
        buffer_type get_deriv_xz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.tertiary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+2);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies. primary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint  );
            return result;
        }
        buffer_type get_deriv_yz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.secondary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies. tertiary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+2);
            return result;
        }

        buffer_type get_deriv_xyz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies. tertiary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+2);
            buffer_type inter_result2= ::partial_deriv<Q,rank>(inter_result,     frequencies.secondary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result2,    frequencies.  primary_grid.get_all_auxiliary_gridpoints(), base_class::data.length(), pos_first_freqpoint  );
            return result;
        }

        auto get_deriv_max() const -> double {
            if constexpr (numberFrequencyDims == 1) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().get_spacing_auxiliary_gridpoints();
                const buffer_type deriv = get_deriv_x();
                const vec<Q> deriv_vec = vec<Q>(deriv.begin(), deriv.end());
                double max = (::power2(deriv_vec * dtb * (1/Kmax))
                ).max_norm();
                return max;
            }
            else if constexpr (numberFrequencyDims == 2) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().get_spacing_auxiliary_gridpoints();
                double dtf = frequencies.get_freqGrid_f().get_spacing_auxiliary_gridpoints();
                double max_K2 = (  ::power2(get_deriv_y() * dtf * (1/Kmax))
                                 + ::power2(get_deriv_x() * dtb * (1/Kmax))
                ).max_norm();
                return max_K2;
            }
            else if constexpr (numberFrequencyDims == 3) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().get_spacing_auxiliary_gridpoints();
                double dtf = frequencies.get_freqGrid_f().get_spacing_auxiliary_gridpoints();
                double max_K3 = (  ::power2(get_deriv_z() * dtf * (1/Kmax))
                                 + ::power2(get_deriv_y() * dtf * (1/Kmax))
                                 + ::power2(get_deriv_x() * dtb * (1/Kmax))
                ).max_norm();
                return max_K3;
            }
            else assert(false);

        }
        auto get_curvature_max() const -> double {
            if constexpr (numberFrequencyDims == 1) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().get_spacing_auxiliary_gridpoints();
                buffer_type curvature = get_deriv_xx();
                double max = (     ::power2(curvature*dtb*dtb*(1/Kmax))
                ).max_norm();
                return max;
            }
            else if constexpr (numberFrequencyDims == 2) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().get_spacing_auxiliary_gridpoints();
                double dtf = frequencies.get_freqGrid_f().get_spacing_auxiliary_gridpoints();
                double max_K2 = (     ::power2(get_deriv_xx()*dtb*dtb*(1/Kmax))
                                    + ::power2(get_deriv_yy()*dtf*dtf*(1/Kmax))
                                    +(::power2(get_deriv_xy()*dtb*dtf*(1/Kmax)) )*2.
                ).max_norm();
                return max_K2;
            }
            else if constexpr (numberFrequencyDims == 3) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().get_spacing_auxiliary_gridpoints();
                double dtf = frequencies.get_freqGrid_f().get_spacing_auxiliary_gridpoints();
                double max_K3 = (     ::power2(get_deriv_xx()*dtb*dtb*(1/Kmax))
                                    + ::power2(get_deriv_yy()*dtf*dtf*(1/Kmax))
                                    + ::power2(get_deriv_zz()*dtf*dtf*(1/Kmax))
                                    +(::power2(get_deriv_xy()*dtb*dtf*(1/Kmax))
                                    + ::power2(get_deriv_xz()*dtb*dtf*(1/Kmax))
                                    + ::power2(get_deriv_yz()*dtf*dtf*(1/Kmax)) )*2.
                ).max_norm();
                return max_K3;
            }
            else assert(false);

        }

        template<typename gridType>
        auto shrink_freq_box_impl(const gridType &freqGrid, const double rel_tail_threshold,
                                  const vec<double> &maxabs_along_x, const bool verbose) const -> gridType {
            const size_t number_of_frequency_points_along_w = maxabs_along_x.size();
            assert(freqGrid.get_all_frequencies().size() == maxabs_along_x.size());

            gridType frequencies_new = freqGrid;

            if (maxabs_along_x.max_norm() < 1e-20) return frequencies_new; // don't shrink if there is no data yet

            int index = -1;
            while (true) {
                if (maxabs_along_x[number_of_frequency_points_along_w - 1 - (index + 1)] >= rel_tail_threshold) break;
                index++;
            }
            if (index > -1) { // if the frequency box is too big, shrink to appropriate size
                double t_belowthresh = freqGrid.get_auxiliary_gridpoint(
                        index); // auxiliary frequency point before passing threshold
                double t_abovethresh = freqGrid.get_auxiliary_gridpoint(
                        index + 1); // auxiliary frequency point after  passing threshold
                double h = (rel_tail_threshold - maxabs_along_x[number_of_frequency_points_along_w - 1 - index]) * (t_abovethresh - t_belowthresh) /
                           (maxabs_along_x[number_of_frequency_points_along_w - 1 - (index + 1)] - maxabs_along_x[number_of_frequency_points_along_w - 1 - index]);
                const double safety = 0.9;
                frequencies_new.set_w_upper(std::abs(freqGrid.frequency_from_t(t_belowthresh + h * safety)));
            } else if (index == -1) { // if data on outermost grid point is too big, then enlarge the box OR print warning
                //double t_upper_new = 1 - maxmax*rel_tail_threshold * (1-freqGrid.t_upper) / (maxabs_along_x[maxabs_along_x.size() -1]);
                //double w_upper_new = frequencies_new.frequency_from_t(t_upper_new);
                //frequencies_new.set_w_upper(w_upper_new);
            }
            frequencies_new.initialize_grid();

            return frequencies_new;
        }

        auto shrink_freq_box(const double rel_tail_threshold, bool verbose = true, double constant = 0.) const -> frequencyGrid_type;

    };


/************************************ MEMBER FUNCTIONS OF THE VERTEX Data Container************************************/



    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type>
    auto DataContainer<Q, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type>::get_VertexFreqGrid() const -> const frequencyGrid_type& {
        return frequencies;
    }


    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type>
    void DataContainer<Q, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type>::set_VertexFreqGrid(const frequencyGrid_type frequencyGrid) {
        frequencies = frequencyGrid;
    }





    template<typename Q, size_t rank, my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, typename frequencyGrid_type>
    auto DataContainer<Q, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type>::shrink_freq_box(const double rel_tail_threshold, const bool verbose, const double constant) const -> frequencyGrid_type {
        if constexpr(numberFrequencyDims == 1) {
            frequencyGrid_type frequencies_new = frequencies;
            typename base_class::buffer_type data_tmp = base_class::data;
            double maxmax = (data_tmp - constant).max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint) * (1 / maxmax);

            const double wmax_old = frequencies_new.  primary_grid.w_upper;
            frequencies_new.  primary_grid = shrink_freq_box_impl(frequencies.  primary_grid, rel_tail_threshold, maxabs_along_w, verbose);
        if (verbose and mpi_world_rank() == 0)
            std::cout << "Shrinking frequency box in direction w from " << wmax_old << " to " << frequencies_new.  primary_grid.w_upper << std::endl;

        return frequencies_new;
        }
        else if constexpr(numberFrequencyDims == 2) {
            frequencyGrid_type frequencies_new = frequencies;

            double maxmax = base_class::data.max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint  ) * (1 / maxmax);
            vec<double> maxabs_along_v = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint+1) * (1 / maxmax);

            const double wmax_old = frequencies_new.  primary_grid.w_upper;
            const double vmax_old = frequencies_new.secondary_grid.w_upper;
            frequencies_new.  primary_grid = shrink_freq_box_impl(frequencies.  primary_grid, rel_tail_threshold, maxabs_along_w, verbose);
            if (verbose and mpi_world_rank() == 0)
                std::cout << "Shrinking frequency box in direction w from " << wmax_old << " to " << frequencies_new.  primary_grid.w_upper << std::endl;
            if (GRID != 2) {
                frequencies_new.secondary_grid = shrink_freq_box_impl(frequencies.secondary_grid, rel_tail_threshold, maxabs_along_v, verbose);
                if (verbose and mpi_world_rank() == 0)
                    std::cout << "Shrinking frequency box in direction v from " << vmax_old << " to " << frequencies_new.secondary_grid.w_upper << std::endl;
            }

            return frequencies_new;
        }
        else if constexpr(numberFrequencyDims == 3) {
            frequencyGrid_type frequencies_new = frequencies;

            double maxmax = base_class::data.max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint  ) * (1 / maxmax);
            vec<double> maxabs_along_v = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint+1) * (1 / maxmax);
            vec<double> maxabs_along_vp= maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint+2) * (1 / maxmax);

            const double wmax_old = frequencies_new.  primary_grid.w_upper;
            const double vmax_old = frequencies_new.secondary_grid.w_upper;
            const double vpmax_old= frequencies_new. tertiary_grid.w_upper;
            frequencies_new.  primary_grid = shrink_freq_box_impl(frequencies.  primary_grid, rel_tail_threshold, maxabs_along_w, verbose);
            if (verbose and mpi_world_rank() == 0)
                std::cout << "Shrinking frequency box in direction w from " << wmax_old << " to " << frequencies_new.  primary_grid.w_upper << std::endl;
            if (GRID != 2) {
                frequencies_new.secondary_grid = shrink_freq_box_impl(frequencies.secondary_grid, rel_tail_threshold, maxabs_along_v, verbose);
                if (verbose and mpi_world_rank() == 0)
                    std::cout << "Shrinking frequency box in direction v from " << vmax_old << " to " << frequencies_new.secondary_grid.w_upper << std::endl;
                frequencies_new.tertiary_grid = shrink_freq_box_impl(frequencies.tertiary_grid, rel_tail_threshold, maxabs_along_vp, verbose);
                if (verbose and mpi_world_rank() == 0)
                    std::cout << "Shrinking frequency box in direction v from " << vmax_old << " to " << frequencies_new.tertiary_grid.w_upper << std::endl;
            }

            return frequencies_new;
        }
        else {
            assert(false);
            return frequencyGrid_type(0);
        }
    }

//}

#endif //FPP_MFRG_VERTEX_DATA_CONTAINER_H
