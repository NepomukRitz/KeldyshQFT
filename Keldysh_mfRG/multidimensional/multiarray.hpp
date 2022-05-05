#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

// code originally from Marc Ritter
// slightly adapted by us

#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>

#include "Eigen/Dense"

#include "ranged_view.hpp"
#include "../utilities/template_utils.hpp"

#ifndef NDEBUG
#define MULTIARRAY_CHECK_BOUNDS
#endif

namespace multidimensional
{
    template <typename T1, typename T2 = T1>
    constexpr static bool abs_compare(T1 lhs, T2 rhs)
    {
        return (std::abs(lhs) < std::abs(rhs));
    }

    // Forward declarations for templates
    template <typename T, size_t depth>
    class multiarray;

    template <typename T, size_t depth>
    bool operator==(const multiarray<T, depth> &lhs, const multiarray<T, depth> &rhs);

    // template <typename Q, typename L, typename R, size_t depth>
    // auto elementwise(const Q &op, const multiarray<L, depth> &lhs, const multiarray<R, depth> &rhs);

    template <typename T, size_t depth>
    class multiarray
    {
    public:
        // === typedefs and constants ===
        using value_type = T;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using reference = value_type &;
        using const_reference = const value_type &;
        using pointer = value_type *;
        using const_pointer = const value_type *;
        using index_type = std::array<my_index_t, depth>;
        using dimensions_type = std::array<size_type, depth>;
        // iterator typedefs follow in the iterator section

    private:
        dimensions_type m_length;       // contains the number of points in each dimension
        dimensions_type get_cumul_length(const dimensions_type& length) const {
            dimensions_type result;
            result[depth-1] = 1;
            for (int i = (int)depth-2; i >= 0; i--) {
                result[i] = result[i+1] * length[i+1];
            }
            return result;
        }
        dimensions_type m_length_cumulative = get_cumul_length(m_length);

    public:
        // ensure we have at least 1 index
        static_assert(depth > 0, "A 0-dimensional multiarray is pointless.");

        // === members ===

        // using buffer_type = std::vector<T>;
        using buffer_type = Eigen::Array<T, Eigen::Dynamic, 1>;
        buffer_type elements;

        // === swap for copy-and-swap ===
        friend void swap(
                multiarray<value_type, depth> &lhs,
                multiarray<value_type, depth> &rhs) noexcept
        {
            using std::swap;
            swap(lhs.m_length, rhs.m_length);
            swap(lhs.elements, rhs.elements);
        }

        /// === constructors ===
        constexpr multiarray() noexcept
                : m_length{}, elements()
        {
        }

        explicit multiarray(dimensions_type length, const T &value = T())
                : m_length(std::move(length)), elements(_flat_size())
        {
            elements.setConstant(value);
        }

        multiarray(dimensions_type length, const buffer_type &elements)
                : m_length(std::move(length)), elements(elements)
        {
            check_size();
        }

        template <
                typename container,
                std::enable_if_t<
                        std::is_same_v<typename container::value_type, value_type> &&
                        !std::is_same_v<container, buffer_type>,
                        bool> = true>
        multiarray(dimensions_type length, const container &elements)
                : m_length(std::move(length)),
                  elements(Eigen::Map<const buffer_type>(elements.data(), elements.size()))
        {
            if ((size_t)elements.size() != this->size())
            {
                throw std::invalid_argument(
                        "Incosistent size of length description and data in "
                        "multiarray::multiarray.");
            }
        }

        multiarray(dimensions_type length, buffer_type &&elements)
                : m_length(std::move(length)), elements(std::move(elements))
        {
            check_size();
        }

        /// Move & copy constructors and assignment operators
        multiarray(const multiarray<T, depth> &) = default;
        multiarray(multiarray<T, depth> &&) = default;

        multiarray<T, depth> &operator=(const multiarray<T, depth> &) = default;
        multiarray<T, depth> &operator=(multiarray<T, depth> &&) = default;

        /// === iterators ===
        // using iterator = typename buffer_type::iterator;
        // using const_iterator = typename buffer_type::const_iterator;
        // using reverse_iterator = typename buffer_type::reverse_iterator;
        // using const_reverse_iterator = typename buffer_type::const_reverse_iterator;
        using iterator = pointer;
        using const_iterator = const_pointer;

        constexpr iterator begin() noexcept
        {
            return elements.data();
        }
        constexpr iterator end() noexcept
        {
            return elements.data() + elements.size();
        }

        constexpr const_iterator begin() const noexcept
        {
            return elements.data();
        }
        constexpr const_iterator end() const noexcept
        {
            return elements.data() + elements.size();
        }

        constexpr const_iterator cbegin() const noexcept
        {
            return elements.data();
        }
        constexpr const_iterator cend() const noexcept
        {
            return elements.data() + elements.size();
        }

        constexpr auto range(const index_type &start, const index_type &end)
        {
            if (check_bounds(start) && check_bounds_end(end))
            {
                return make_range(*this, flat_index(start), flat_index(end));
            }
            else
            {
                throw std::invalid_argument("Passed invalid indices during creation of a multiarray range");
            }
        }

        constexpr auto range(const index_type &start, const index_type &end) const
        {
            if (check_bounds(start) && check_bounds_end(end))
            {
                return make_const_range(*this, flat_index(start), flat_index(end));
            }
            else
            {
                throw std::invalid_argument("Passed invalid indices during creation of a multiarray range");
            }
        }

        /// return segment including start and end
        constexpr auto eigen_segment(const index_type &start, const index_type &end)
        {
            if (check_bounds(start) && check_bounds_end(end))
            {
                const auto flat_start = flat_index(start);
                const auto flat_end = flat_index(end);
                return elements.segment(flat_start, flat_end - flat_start);
            }
            else
            {
                throw std::invalid_argument("Passed invalid indices during creation of a multiarray range");
            }
        }

        constexpr auto eigen_segment(const index_type &start, const index_type &end) const
        {
            if (check_bounds(start) && check_bounds_end(end))
            {
                const auto flat_start = flat_index(start);
                const auto flat_end = flat_index(end);
                return elements.segment(flat_start, flat_end - flat_start);
            }
            else
            {
                throw std::invalid_argument("Passed invalid indices during creation of a multiarray range");
            }
        }


        bool is_same_length(const multiarray<value_type,depth>& other) const {
            return (m_length == other.m_length);
        }

        template <std::size_t pos_first_freq_index, std::size_t freqrank, std::size_t vecsize, typename... Types,
                typename std::enable_if_t<(sizeof...(Types) == depth) and (are_all_integral<size_t, Types...>::value) and (are_all_unsigned<size_t, Types...>::value) and (pos_first_freq_index + freqrank < depth), bool> = true>
        constexpr auto at_vectorized(Types & ...i  ) const -> Eigen::Matrix<T,vecsize,1>{
#ifndef NDEBUG
            my_index_t it = 0;
            ((assert((my_index_t)i < m_length[it]), it ++),...);
#endif
            const auto flat_start = flat_index(index_type({static_cast<size_t>(i)...}));;
            return elements.template segment<vecsize>(flat_start);
        }
        template <std::size_t pos_first_freq_index, std::size_t freqrank, std::size_t vecsize>
        constexpr auto at_vectorized(const index_type & idx  ) const -> Eigen::Matrix<T,vecsize,1>{
#ifndef NDEBUG
            for (size_t it = 0; it < depth; it++) {
                assert(idx[it] < m_length[it]);
            }
#endif
            const auto flat_start = flat_index(idx);
            return elements.template segment<vecsize>(flat_start);
        }
        template <std::size_t vecsize>
        constexpr auto at_vectorized(const std::size_t flat_idx  ) const -> Eigen::Matrix<T,vecsize,1>{
            assert(flat_idx+vecsize <= (size_t)elements.size());

            return elements.template segment<vecsize>(flat_idx);
        }

        template <my_index_t num_first_dims>
        size_t get_flatindex_ini(const index_type & index) const {
            size_t result = 0;
            for (my_index_t i = 0; i <= num_first_dims; i++) {
                result += index[i] * m_length_cumulative[i];
            }
            return result;
        }

        template <my_index_t numberFrequencyDims, my_index_t pos_first_freqpoint, my_index_t vecsize, my_index_t sample_size>
        auto get_values(const index_type& index) const -> Eigen::Matrix<T, vecsize, numberFrequencyDims == 1 ? sample_size : (numberFrequencyDims == 2 ? sample_size*sample_size : sample_size*sample_size*sample_size)>
                {
#ifndef NDEBUG
            assert(check_bounds(index));
#endif
            Eigen::Matrix<T, vecsize, numberFrequencyDims == 1 ? sample_size : (numberFrequencyDims == 2 ? sample_size*sample_size : sample_size*sample_size*sample_size)> result;
            const size_t flat_ini = get_flatindex_ini<pos_first_freqpoint+numberFrequencyDims>(index);

            if constexpr (numberFrequencyDims == 1) {
                for (my_index_t i = 0; i < sample_size; i++) {
                    const auto res = at_vectorized<vecsize>(flat_ini + m_length_cumulative[pos_first_freqpoint]*i);
                    result.col(i) = res;
                }
            }
            else if constexpr(numberFrequencyDims == 2) {
                for (my_index_t i = 0; i < sample_size; i++) {
                    for (my_index_t j = 0; j < sample_size; j++) {
                        result.col(i * sample_size + j) = at_vectorized<vecsize>(flat_ini + m_length_cumulative[pos_first_freqpoint]*i + m_length_cumulative[pos_first_freqpoint+1]*j);

                    }
                }
            }
            else if constexpr (numberFrequencyDims == 3) {
                for (my_index_t i = 0; i < sample_size; i++) {
                    for (my_index_t j = 0; j < sample_size; j++) {
                        for (my_index_t l = 0; l < sample_size; l++) {
                            result.col(i * sample_size*sample_size + j * sample_size + l) = at_vectorized<vecsize>( flat_ini + m_length_cumulative[pos_first_freqpoint] * i + m_length_cumulative[pos_first_freqpoint + 1] * j + m_length_cumulative[pos_first_freqpoint + 2] * l);
                        }
                    }
                }
            }
            else {
                assert(false); // numberFrequencyDims > 3 not supported
            }
            return result;
        }

        /// === public member functions ===

        /// random element access
        // The type index_type ensures we have exactly the right number of indices.
        T &at(const index_type &index);
        const T &at(const index_type &index) const;

        // The enable_if<> checks whether the index has the right length.
        // This is necessary as initialization of index_type with an initializer list of length < depth will just set the remaining indices to zero.
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        T &at(const Types &...i)
        {
            return at(index_type({static_cast<size_type>(i)...}));
        }
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        const T &at(const Types &...i) const
        {
            return at(index_type({static_cast<size_type>(i)...}));
        }
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        T &operator()(const Types &...i)
        {
            return at(i...);
        }
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        const T &operator()(const Types &...i) const
        {
            return at(i...);
        }

        T *data() noexcept
        {
            return elements.data();
        }

        const T *data() const noexcept
        {
            return elements.data();
        }

        buffer_type get_elements() const {
            return elements;
        }

        /// flat access

        T &flat_at(size_type i)
        {
            return elements.coeffRef(i);
        }

        const T &flat_at(size_type i) const
        {
            return elements.coeffRef(i);
        }

        T &operator[](size_type i)
        {
            return flat_at(i);
        }

        const T &operator[](size_type i) const
        {
            return flat_at(i);
        }

        /// elementwise arithmetics-assignment op's
        template <typename Q, typename R>
        multiarray<T, depth> &elementwise_map_assign(
                Q op,
                const multiarray<R, depth> &rhs)
        {
            if (rhs.length() != length())
            {
                throw std::length_error("Cannot perform pairwise operations on multiarrays of different length.");
            }
            for (size_type i = 0; i < elements.size(); i++)
            {
                elements[i] = op(elements[i], rhs.elements[i]);
            }
            return *this;
        }

        multiarray<T, depth> &operator+=(
                const multiarray<T, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements += rhs.elements;
            return *this;
        }

        multiarray<T, depth> &operator-=(
                const multiarray<T, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements -= rhs.elements;
            return *this;
        }

        multiarray<T, depth> &operator*=(
                const multiarray<T, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements *= rhs.elements;
            return *this;
        }

        multiarray<T, depth> &operator/=(
                const multiarray<T, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements /= rhs.elements;
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator+=(
                const multiarray<R, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements += rhs.elements.template cast<T>();
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator-=(
                const multiarray<R, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements -= rhs.elements.template cast<T>();
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator*=(
                const multiarray<R, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements *= rhs.elements.template cast<T>();
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator/=(
                const multiarray<R, depth> &rhs)
        {
            assert(is_same_length(rhs));
            elements /= rhs.elements.template cast<T>();
            return *this;
        }

        /// scalar arithmetic assignment
        template <typename Q, typename R>
        multiarray<T, depth> &scalar_map_assign(Q op, const R &rhs) noexcept(noexcept(op(T(), R())))
        {
            for (auto &e : elements)
            {
                e = op(e, rhs);
            }
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator+=(const R &rhs) noexcept(noexcept(T() + R()))
        {
            elements += rhs;
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator-=(const R &rhs) noexcept(noexcept(T() - R()))
        {
            elements -= rhs;
            return *this;
            ;
        }

        template <typename R>
        multiarray<T, depth> &operator*=(const R &rhs) noexcept(noexcept(T() * R()))
        {
            elements *= rhs;
            return *this;
        }

        template <typename R>
        multiarray<T, depth> &operator/=(const R &rhs) noexcept(noexcept(T() - R()))
        {
            elements /= rhs;
            return *this;
        }

        /// other function related to arithmetic

        multiarray<T,depth> abs() const {
            return transform([] (const T& x) {return static_cast<T>(std::abs(x));}, *this);
        }

        T max() const
        {
            return elements.maxCoeff();
        }

        T min() const
        {
            return elements.minCoeff();
        }

        T maxabs() const
        {
            return std::abs(*std::max_element(begin(), end(), abs_compare<value_type>));
        }

        T minabs() const
        {
            return std::abs(*std::min_element(begin(), end(), abs_compare<value_type>));
        }

        double max_norm() const {
                return std::abs(maxabs());
            };

        template<int p> double lpNorm() const {
            return elements.template lpNorm<p>();
        }


        bool isfinite() const
        {
            for (const auto &c : *this)
            {
                if (!std::isfinite(c))
                {
                    return false;
                }
            }
            return true;
        }

        /// non-member op's
        friend bool operator==<>(const multiarray<T, depth> &, const multiarray<T, depth> &);

        size_type size() const noexcept
        {
            return elements.size();
        }

        constexpr const dimensions_type &length() const noexcept
        {
            return m_length;
        }

        static constexpr size_type get_depth() noexcept
        {
            return depth;
        }

    private:
        // === private members ===

        // === private member functions (mostly helpers) ===
        size_type _flat_size() const noexcept
        {
            size_type res = 1;
            for (auto &l : m_length)
            {
                res *= l;
            }
            return res;
        }

        void check_size() const
        {
#if DEBUG
            if (_flat_size() != elements.size())
            {
                std::cerr << "Flat size: " << _flat_size() << std::endl
                          << "elements.size(): " << elements.size() << std::endl;
                throw std::logic_error("multiarrray::operator= got an invalid parameter and cannot construct a consistent object.");
            }
#endif
        }

        size_type flat_index(const index_type &index) const noexcept;
        bool check_bounds(const index_type &index) const noexcept;
        bool check_bounds_end(const index_type &index) const noexcept;
    };

    // === template definitions ===

    template <typename T, size_t depth>
    typename multiarray<T, depth>::size_type // <- This should be size_t. Compiler is too stupid to figure that out himself.
    multiarray<T, depth>::flat_index(const index_type &index) const noexcept
    {
        size_type res = index[0];
        for (size_t i = 1; i < depth; i++)
        {
            res *= m_length[i];
            res += index[i];
        }

        return res;
    }

    template <typename T, size_t depth>
    bool multiarray<T, depth>::check_bounds(const index_type &index) const noexcept
    {
#ifdef MULTIARRAY_CHECK_BOUNDS
        for (size_t i = 0; i < depth; i++)
        {
            if (index[i] >= m_length[i] || index[i] < 0)
            {
                return false;
            }
        }
#endif
        return true;
    }

    template <typename T, size_t depth>
    bool multiarray<T, depth>::check_bounds_end(const index_type &index) const noexcept
    {
#ifdef MULTIARRAY_CHECK_BOUNDS
        for (size_t i = 0; i < depth; i++)
        {
            if (index[i] > m_length[i] || index[i] < 0)
            {
                return false;
            }
        }
#endif
        return true;
    }

    template <typename T, size_t depth>
    T &multiarray<T, depth>::at(const index_type &index)
    {
        if (check_bounds(index))
        {
            return flat_at(flat_index(index));
        }
        else
        {
            assert(false);
            throw std::out_of_range("Attempted to access multiarray element at invalid index.");
        }
    }

    template <typename T, size_t depth>
    const T &multiarray<T, depth>::at(const index_type &index) const
    {
        if (check_bounds(index))
        {
            return flat_at(flat_index(index));
        }
        else
        {
            assert(false);
            throw std::out_of_range("Attempted to access multiarray element at invalid index.");
        }
    }

    // === operator implementation ===
    template <typename T, size_t depth>
    bool operator==(const multiarray<T, depth> &lhs, const multiarray<T, depth> &rhs)
    {
        return (lhs.m_length == rhs.m_length) && (lhs.elements == rhs.elements).all();
    }

    template <typename T, size_t depth>
    bool operator!=(const multiarray<T, depth> &lhs, const multiarray<T, depth> &rhs)
    {
        return !(lhs == rhs);
    }

    // === arithmetic operators ===

    // unary op
    template <
            typename Q, typename L, size_t depth,
            std::enable_if_t<std::is_same_v<std::result_of_t<Q && (L &&)>, L>, bool> = true>
    auto transform(Q op, multiarray<L, depth> lhs) noexcept(noexcept(op(std::declval<L>())))
    {
    for (size_t i = 0; i < lhs.size(); i++)
    {
    lhs.elements[i] = op(lhs.elements[i]);
    }
    return lhs;
    }

    // binary op for two multiarrays
    template <
            typename Q,
            typename L, typename R,
            size_t depth,
            std::enable_if_t<std::is_same_v<std::result_of_t<Q && (L &&, R &&)>, L>, bool> = true>
    auto elementwise(Q op, multiarray<L, depth> lhs, const multiarray<R, depth> &rhs)
    {
        if (lhs.length() != rhs.length())
        {
            throw std::length_error("Cannot perform pairwise operations on multiarrays of different length.");
        }
        for (size_t i = 0; i < lhs.size(); i++)
        {
            lhs.elements[i] = op(lhs.elements[i], rhs.elements[i]);
        }
        return lhs;
    }

    // binary op for two multiarrays
    template <
            typename Q,
            typename L, typename R,
            size_t depth,
            std::enable_if_t<!std::is_same_v<std::result_of_t<Q && (L &&, R &&)>, L>, bool> = true>
    auto elementwise(Q op, const multiarray<L, depth> &lhs, const multiarray<R, depth> &rhs)
    {
        if (lhs.length() != rhs.length())
        {
            throw std::length_error("Cannot perform pairwise operations on multiarrays of different length.");
        }
        multiarray<std::result_of_t<Q && (L &&, R &&)>, depth> res(lhs.length());
        for (size_t i = 0; i < res.size(); i++)
        {
            res.elements[i] = op(lhs.elements[i], rhs.elements[i]);
        }
        return res;
    }

    // binary op for a multiarray and a scalar
    template <
            typename Q, typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<std::result_of_t<Q && (L &&, R &&)>, L>, bool> = true>
    auto scalar_map(Q op, multiarray<L, depth> lhs, const R &rhs) noexcept(noexcept(op(std::declval<L>(), std::declval<R>())))
    {
        for (size_t i = 0; i < lhs.size(); i++)
        {
            lhs.elements[i] = op(lhs.elements[i], rhs);
        }
        return lhs;
    }

    // binary op for a multiarray and a scalar
    template <
            typename Q, typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<std::result_of_t<Q && (L &&, R &&)>, L>, bool> = true>
    auto scalar_map(Q op, const multiarray<L, depth> &lhs, const R &rhs)
    {
        multiarray<std::result_of_t<Q && (L &&, R &&)>, depth> res(
                lhs.length());
        for (size_t i = 0; i < res.size(); i++)
        {
            res.elements[i] = op(lhs.elements[i], rhs);
        }
        return res;
    }

    // Binary operators
    // Addition

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator+(const multiarray<L, depth> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return l + r; },
                          lhs, rhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator+(multiarray<L, depth> lhs, const R &rhs)
    {
        return lhs += rhs;
    }

    template <typename L, typename R, size_t depth>
    auto operator+(const R &lhs, const multiarray<L, depth> &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return r + l; },
                          rhs, lhs);
    }

    template <
            typename L, typename R,
            size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator+(
            const multiarray<L, depth> &lhs,
            const multiarray<R, depth> &rhs)
    {
        return elementwise([](const L &l, const R &r)
                           { return l + r; },
                           lhs, rhs);
    }

    template <
            typename L, typename R,
            size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator+(
            multiarray<L, depth> lhs,
            const multiarray<R, depth> &rhs)
    {
        return lhs += rhs;
    }

    // Subtraction

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator-(const multiarray<L, depth> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return l - r; },
                          lhs, rhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() - std::declval<R>()), L>, bool> = true>
    auto operator-(multiarray<L, depth> lhs, const R &rhs)
    {
        return lhs -= rhs;
    }

    template <typename L, typename R, size_t depth>
    auto operator-(const R &lhs, const multiarray<L, depth> &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return r - l; },
                          rhs, lhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator-(
            const multiarray<L, depth> &lhs,
            const multiarray<R, depth> &rhs)
    {
        return elementwise([](const L &l, const R &r)
                           { return l - r; },
                           lhs, rhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() - std::declval<R>()), L>, bool> = true>
    auto operator-(
            multiarray<L, depth> lhs,
            const multiarray<R, depth> &rhs)
    {
        return lhs -= rhs;
    }

    // Multiplication

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator*(const multiarray<L, depth> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return l * r; },
                          lhs, rhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() * std::declval<R>()), L>, bool> = true>
    auto operator*(multiarray<L, depth> lhs, const R &rhs)
    {
        return lhs *= rhs;
    }

    template <typename L, typename R, size_t depth>
    auto operator*(const R &lhs, const multiarray<L, depth> &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return r * l; },
                          rhs, lhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator*(
            const multiarray<L, depth> &lhs,
            const multiarray<R, depth> &rhs)
    {
        return elementwise([](const L &l, const R &r)
                           { return l * r; },
                           lhs, rhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() * std::declval<R>()), L>, bool> = true>
    auto operator*(multiarray<L, depth> lhs, const multiarray<R, depth> &rhs)
    {
        return lhs *= rhs;
    }

    // Division

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator/(const multiarray<L, depth> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return l / r; },
                          lhs, rhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() / std::declval<R>()), L>, bool> = true>
    auto operator/(multiarray<L, depth> lhs, const R &rhs)
    {
        return lhs /= rhs;
    }

    template <typename L, typename R, size_t depth>
    auto operator/(const R &lhs, const multiarray<L, depth> &rhs)
    {
        return scalar_map([](const L &l, const R &r)
                          { return r / l; },
                          rhs, lhs);
    }

    template <
            typename L, typename R, size_t depth,
            std::enable_if_t<!std::is_same_v<decltype(std::declval<L>() + std::declval<R>()), L>, bool> = true>
    auto operator/(
            const multiarray<L, depth> &lhs,
            const multiarray<R, depth> &rhs)
    {
        return elementwise([](const L &l, const R &r)
                           { return l / r; },
                           lhs, rhs);
    }

    template <
            typename L, typename R,
            size_t depth,
            std::enable_if_t<std::is_same_v<decltype(std::declval<L>() / std::declval<R>()), L>, bool> = true>
    auto operator/(
            multiarray<L, depth> lhs,
            const multiarray<R, depth> &rhs)
    {
        return lhs /= rhs;
    }

    // unary minus
    template <typename value_type, size_t depth>
    auto operator-(multiarray<value_type, depth> &&rhs)
    {
        return (-1) * std::move(rhs);
    }

    template <typename value_type, size_t depth>
    auto operator-(const multiarray<value_type, depth> &rhs)
    {
        return (-1) * rhs;
    }
} // namespace multidimensional

#endif
