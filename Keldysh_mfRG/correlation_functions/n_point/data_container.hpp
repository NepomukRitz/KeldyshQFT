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
template <typename Q> class State; // forward declaration of State
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
        //friend class State<Q>;
        template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const rvert<T>& rhs);
        template<typename T> friend rvert<T> rvert<T>::operator+= (const rvert<T>& rhs);
        template<typename T> friend rvert<T> operator- (rvert<T> lhs, const rvert<T>& rhs);
        template<typename T> friend rvert<T> rvert<T>::operator-= (const rvert<T>& rhs);
        template<typename T> friend rvert<T> operator/ (rvert<T> lhs, const rvert<T>& rhs);
        template<typename T> friend rvert<T> rvert<T>::operator/= (const rvert<T>& rhs);
        template<typename T> friend rvert<T> operator* (rvert<T> lhs, const double& alpha);
        template<typename T> friend rvert<T> rvert<T>::operator*= (double alpha);
        template<typename T> friend rvert<T> operator+ (rvert<T> lhs, const double& alpha);
        template<typename T> friend rvert<T> rvert<T>::operator+= (double alpha);
        friend State<state_datatype> read_state_from_hdf(const H5std_string& filename, const int Lambda_it);


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
            return data.template at_vectorized<pos_first_freq, freqrank, vecsize>(i...);
        }
        template<std::size_t freqrank, std::size_t vecsize>
        auto val_vectorized(const index_type &idx) const -> Eigen::Matrix<Q, vecsize, 1> {
            return data.template at_vectorized<pos_first_freq, freqrank, vecsize>(idx);
        }

        /// Sets a value at a multiIndex
        template<typename... Types, typename std::enable_if_t<
                (sizeof...(Types) == rank) and (are_all_integral<size_t, Types...>::value), bool> = true
        >
        void setvert(const Q value, const Types &... i) { data.at(i...) = value; }

        void setvert(const Q value, const index_type &idx) { data.at(idx) = value; }

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
        friend void result_set_frequency_grids(State<T> &result, Buffer &buffer);

        template<typename T>
        friend void check_FDTs(const State<T> &state, bool verbose);

        friend State<state_datatype> read_state_from_hdf(const H5std_string &filename, const int Lambda_it);

    protected:
        using base_class = dataContainerBase<Q, rank>;
    public:
        using index_type = typename base_class::index_type;
        using dimensions_type = typename base_class::dimensions_type;
        using buffer_type = typename base_class::buffer_type;

        frequencyGrid_type frequencies;    // frequency grid
        DataContainer() = default;

        explicit DataContainer(double Lambda, dimensions_type dims) : frequencies(Lambda), base_class(dims) {};

        /// Functions for getting and setting the frequency grid and its members
        auto get_VertexFreqGrid() const -> const frequencyGrid_type &;

        void set_VertexFreqGrid(const frequencyGrid_type frequencyGrid);


        template<unsigned int idim,
                typename std::enable_if_t<(idim < numberFrequencyDims), bool> = true
                >
        double analyze_tails() const {
            double maxabs_total = base_class::data.max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint+idim);

            return maxabs_along_w[0] / maxabs_total;
        };


        buffer_type get_deriv_x () const {
            buffer_type result = ::partial_deriv<Q,rank>(base_class::data, frequencies.b.ts, base_class::data.length(),  pos_first_freqpoint);
            return result;
        };
        buffer_type get_deriv_y () const {
            buffer_type result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(),  pos_first_freqpoint+1);
            return result;
        };
        buffer_type get_deriv_z () const {
            buffer_type result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(),  pos_first_freqpoint+2);
            return result;
        };

        buffer_type get_deriv_xx() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.b.ts, base_class::data.length(), pos_first_freqpoint  );
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.b.ts, base_class::data.length(), pos_first_freqpoint  );
            return result;
        }
        buffer_type get_deriv_yy() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+1);
            return result;
        }
        buffer_type get_deriv_zz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+2);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+2);
            return result;
        }
        buffer_type get_deriv_xy() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.b.ts, base_class::data.length(), pos_first_freqpoint  );
            return result;
        }
        buffer_type get_deriv_xz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+2);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.b.ts, base_class::data.length(), pos_first_freqpoint  );
            return result;
        }
        buffer_type get_deriv_yz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result    , frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+2);
            return result;
        }

        buffer_type get_deriv_K3_xyz() const {
            buffer_type inter_result = ::partial_deriv<Q,rank>(base_class::data, frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+2);
            buffer_type inter_result2= ::partial_deriv<Q,rank>(inter_result,     frequencies.f.ts, base_class::data.length(), pos_first_freqpoint+1);
            buffer_type result       = ::partial_deriv<Q,rank>(inter_result2,    frequencies.b.ts, base_class::data.length(), pos_first_freqpoint  );
            return result;
        }

        auto get_deriv_max() const -> double {
            if constexpr (numberFrequencyDims == 1) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().dt;
                double max = (::power2(get_deriv_x() * dtb * (1/Kmax))
                ).max_norm();
                return max;
            }
            else if constexpr (numberFrequencyDims == 2) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().dt;
                double dtf = frequencies.get_freqGrid_f().dt;
                double max_K2 = (  ::power2(get_deriv_y() * dtf * (1/Kmax))
                                 + ::power2(get_deriv_x() * dtb * (1/Kmax))
                ).max_norm();
                return max_K2;
            }
            else if constexpr (numberFrequencyDims == 3) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().dt;
                double dtf = frequencies.get_freqGrid_f().dt;
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
                double dtb = frequencies.get_freqGrid_b().dt;
                double max = (     ::power2(get_deriv_xx()*dtb*dtb*(1/Kmax))
                ).max_norm();
                return max;
            }
            else if constexpr (numberFrequencyDims == 2) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().dt;
                double dtf = frequencies.get_freqGrid_f().dt;
                double max_K2 = (     ::power2(get_deriv_xx()*dtb*dtb*(1/Kmax))
                                    + ::power2(get_deriv_yy()*dtf*dtf*(1/Kmax))
                                    +(::power2(get_deriv_xy()*dtb*dtf*(1/Kmax)) )*2.
                ).max_norm();
                return max_K2;
            }
            else if constexpr (numberFrequencyDims == 3) {
                double Kmax = base_class::get_vec().max_norm();
                double dtb = frequencies.get_freqGrid_b().dt;
                double dtf = frequencies.get_freqGrid_f().dt;
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

        auto shrink_freq_box(const double rel_tail_threshold, bool verbose = true) const -> frequencyGrid_type;

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
    auto DataContainer<Q, rank, numberFrequencyDims, pos_first_freqpoint, frequencyGrid_type>::shrink_freq_box(const double rel_tail_threshold,
                                                         const bool verbose) const -> frequencyGrid_type {
        if constexpr(numberFrequencyDims == 1) {
            frequencyGrid_type frequencies_new = frequencies;
            double maxmax = base_class::data.max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint) * (1 / maxmax);

        frequencies_new.b = freqGrid::shrink_freq_box(frequencies.b, rel_tail_threshold, maxabs_along_w, verbose);
        if (verbose and mpi_world_rank() == 0)
            std::cout << "in direction w to " << frequencies_new.b.w_upper << std::endl;

        return frequencies_new;
        }
        else if constexpr(numberFrequencyDims == 2) {
            frequencyGrid_type frequencies_new = frequencies;

            double maxmax = base_class::data.max_norm();
            vec<double> maxabs_along_w = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint  ) * (1 / maxmax);
            vec<double> maxabs_along_v = maxabs(base_class::data, base_class::data.length(), pos_first_freqpoint+1) * (1 / maxmax);


            frequencies_new.b = freqGrid::shrink_freq_box(frequencies.b, rel_tail_threshold, maxabs_along_w, verbose);
            if (verbose and mpi_world_rank() == 0)
                std::cout << "in direction w to " << frequencies_new.b.w_upper << std::endl;
            frequencies_new.f = freqGrid::shrink_freq_box(frequencies.f, rel_tail_threshold, maxabs_along_v, verbose);
            if (verbose and mpi_world_rank() == 0)
                std::cout << "in direction v to " << frequencies_new.f.w_upper << std::endl;

            return frequencies_new;
        }
        else {
            assert(false);
        }
    }

//}

#endif //FPP_MFRG_VERTEX_DATA_CONTAINER_H
