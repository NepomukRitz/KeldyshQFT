#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

// code originally from Marc Ritter
// slightly adapted by us

#include <array>
#include <exception>
#include <functional>
#include <iostream>
#include <vector>

#include "ranged_view.hpp"
#include "../data_structures.h"


#if DEBUG >= 1
#define MULTIARRAY_CHECK_BOUNDS
#endif

namespace multidimensional
{

    // Forward declarations for templates
    template <typename T, size_t depth, typename Allocator>
    class multiarray;

    template <typename T, size_t depth, typename Allocator>
    bool operator==(const multiarray<T, depth, Allocator> &lhs, const multiarray<T, depth, Allocator> &rhs);

    template <typename Q, typename L, typename R, size_t depth, typename Allocator>
    auto elementwise(const Q &op, const multiarray<L, depth, Allocator> &lhs, const multiarray<R, depth, Allocator> &rhs);

    template <typename Q, typename L, size_t depth, typename Allocator>
    auto transform_multiarray(const Q &op, const multiarray<L, depth, Allocator> &arr);

    /**
     * tensor of rank = depth
     * @tparam T            datatype of tensor
     * @tparam depth        rank of tensor
     * @tparam Allocator
     */
    template <typename T, size_t depth, typename Allocator=std::allocator<T>>
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
        using index_type = std::array<size_type, depth>;
        using dimensions_type = std::array<size_type, depth>;
        // iterator typedefs follow in the iterator section

    private:
        dimensions_type m_length;       // contains the number of points in each dimension

    public:
        // ensure we have at least 1 index
        static_assert(depth > 0, "A 0-dimensional multiarray is pointless.");

        // === members ===

        using buffer_type = vec<T>; //std::vector<T, Allocator>;
        buffer_type elements;           // contains data; stored in a vector; accessed via a flattened index

        /// === constructors ===
        multiarray() = delete;

        // constructs a multiarray, initializes elements with value everywhere (standard value = zero)
        // lengths of dimensions handed over in an std::array
        explicit multiarray(const dimensions_type &length, const T &value = T())
                : m_length(length), elements(_flat_size(), value)
        {
        }

        // constructs a multiarray with zeros,
        // lengths of dimensions is handed over as a list of integers
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        explicit multiarray(const Types &... lengths)
                : multiarray(dimensions_type{static_cast<size_t>(lengths)...})
        {
        }

        // constructs a multiarray, copies input vector to elements (if input is of buffer_type)
        // lengths of dimensions handed over in an std::array
        multiarray(dimensions_type length, const buffer_type &elements)
                : m_length(std::move(length)), elements(elements)
        {
            check_size();
        }

        // constructs a multiarray, copies input vector to elements (if input is NOT of buffer_type)
        template <typename container>
        multiarray(dimensions_type length, const container &elements)
                : m_length(std::move(length)), elements(elements.begin(), elements.end())
        {
            check_size();
        }

        // constructs a multiarray, copies input vector to elements (if input is of buffer_type)
        multiarray(dimensions_type length, buffer_type &&elements)
                : m_length(std::move(length)), elements(std::move(elements))
        {
            check_size();
        }

        // Move & copy constructors and assignment operators
        multiarray(const multiarray<T, depth, Allocator> &other)
                : m_length(other.m_length), elements(other.elements)
        {
        }

        template <typename Other_Alloc>
        multiarray(const multiarray<T, depth, Other_Alloc> &other)
                : m_length(other.m_length), elements()
        {
            elements = other.elements;
            check_size();
        }

        multiarray(multiarray<T, depth, Allocator> &&other)
                : m_length(std::move(other.m_length)), elements(std::move(other.elements))
        {
            other.m_length = {0};
            check_size();
        }

        template <typename Other_Alloc>
        multiarray(multiarray<T, depth, Other_Alloc> &&other)
                : m_length(std::move(other.m_length)),
                  elements()
        {
            elements = std::move(other.elements);
            check_size();
        }

        multiarray<T, depth, Allocator> &operator=(const multiarray<T, depth, Allocator> &other)
        {
            return do_copy_assign(other);
        }

        template <typename Other_Alloc>
        multiarray<T, depth, Allocator> &operator=(const multiarray<T, depth, Other_Alloc> &other)
        {
            return do_copy_assign(other);
        }

        multiarray<T, depth, Allocator> &operator=(multiarray<T, depth, Allocator> &&other)
        {
            return do_move_assign(std::move(other));
        }

        template <typename Other_Alloc>
        multiarray<T, depth, Allocator> &operator=(multiarray<T, depth, Other_Alloc> &&other)
        {
            return do_move_assign(std::move(other));
        }

        // === iterators ===
        using iterator = typename buffer_type::iterator;
        using const_iterator = typename buffer_type::const_iterator;
        using reverse_iterator = typename buffer_type::reverse_iterator;
        using const_reverse_iterator = typename buffer_type::const_reverse_iterator;

        constexpr iterator begin() noexcept
        {
            return elements.begin();
        }
        constexpr iterator end() noexcept
        {
            return elements.end();
        }

        constexpr const_iterator begin() const noexcept
        {
            return elements.begin();
        }
        constexpr const_iterator end() const noexcept
        {
            return elements.end();
        }

        constexpr const_iterator cbegin() const noexcept
        {
            return elements.cbegin();
        }
        constexpr const_iterator cend() const noexcept
        {
            return elements.cend();
        }

        constexpr reverse_iterator rbegin() noexcept
        {
            return elements.rbegin();
        }
        constexpr reverse_iterator rend() noexcept
        {
            return elements.rend();
        }

        constexpr const_reverse_iterator rbegin() const noexcept
        {
            return elements.rbegin();
        }
        constexpr const_reverse_iterator rend() const noexcept
        {
            return elements.rend();
        }

        constexpr const_reverse_iterator crbegin() const noexcept
        {
            return elements.crbegin();
        }
        constexpr const_reverse_iterator crend() const noexcept
        {
            return elements.crend();
        }

        constexpr auto range(const index_type &start, const index_type &end)
        {
            if(check_bounds(start) && check_bounds_end(end))
            {
                return make_range(elements, flat_index(start), flat_index(end));
            }
            else
            {
                throw std::invalid_argument("Passed invalid indices during creation of a multiarray range");
            }
        }

        // === public member functions ===

        // random element access
        // The type index_type ensures we have exactly the right number of indices.
        T &at(const index_type &index);
        const T &at(const index_type &index) const;

        // The enable_if<> checks whether the index has the right length.
        // This is necessary as initialization of index_type with an initializer list of length < depth will just set the remaining indices to zero.
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        T &at(const Types &... i)
        {
            return at(index_type({static_cast<size_t>(i)...}));
        }
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        const T &at(const Types &... i) const
        {
            return at(index_type({static_cast<size_t>(i)...}));
        }
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        T &operator()(const Types &... i)
        {
            return at(i...);
        }
        template <typename... Types,
                typename std::enable_if_t<sizeof...(Types) == depth, bool> = true>
        const T &operator()(const Types &... i) const
        {
            return at(i...);
        }

        T *data()
        {
            return elements.data();
        }

        // flat access

        T &flat_at(size_type i)
        {
            return elements.at(i);
        }

        const T &flat_at(size_type i) const
        {
            return elements.at(i);
        }

        T &operator[](size_type i)
        {
            return flat_at(i);
        }

        const T &operator[](size_type i) const
        {
            return flat_at(i);
        }


        const auto get_elements(const int i1, const int i2)
        {
            size_type flat_size = _flat_size();
            if(std::abs(i1) <= flat_size && std::abs(i2) <= flat_size)
            {
                auto it1 = (i1 >= 0) ? elements.begin() : elements.end();
                auto it2 = (i2 >= 0) ? elements.begin() : elements.end();
                buffer_type subvector (it1 + i1, it2 + i2 + 1);
                return subvector;
            }
            else
            {
                throw std::invalid_argument("Passed invalid indices during creation of a subarray");
            }
        }

        // elementwise arithmetics-assignment op's
        template <typename Q, typename R>
        multiarray<T, depth, Allocator> &elementwise_map_assign(
                const Q &op,
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

        template <typename R>
        multiarray<T, depth, Allocator> &operator+=(const multiarray<R, depth, Allocator> &rhs)
        {
            return elementwise_map_assign([](const T &l, const R &r) { return l + r; }, rhs);
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator-=(const multiarray<R, depth, Allocator> &rhs)
        {
            return elementwise_map_assign([](const T &l, const R &r) { return l - r; }, rhs);
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator*=(const multiarray<R, depth, Allocator> &rhs)
        {
            return elementwise_map_assign([](const T &l, const R &r) { return l * r; }, rhs);
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator/=(const multiarray<R, depth, Allocator> &rhs)
        {
            return elementwise_map_assign([](const T &l, const R &r) { return l / r; }, rhs);
        }

        // scalar arithmetic assignment
        template <typename Q, typename R>
        multiarray<T, depth, Allocator> &scalar_map_assign(const Q &op, const R &rhs)
        {
            for (auto &e : elements)
            {
                e = op(e, rhs);
            }
            return *this;
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator+=(const R &rhs)
        {
            return scalar_map_assign([](const T &l, const R &r) { return l + r; }, rhs);
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator-=(const R &rhs)
        {
            return scalar_map_assign([](const T &l, const R &r) { return l - r; }, rhs);
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator*=(const R &rhs)
        {
            return scalar_map_assign([](const T &l, const R &r) { return l * r; }, rhs);
        }

        template <typename R>
        multiarray<T, depth, Allocator> &operator/=(const R &rhs)
        {
            return scalar_map_assign([](const T &l, const R &r) { return l / r; }, rhs);
        }

        multiarray<T, depth, Allocator> &conj()
        {
            return transform_multiarray(myconj, *this);
        }
        multiarray<T, depth, Allocator> &real()
        {
            return transform_multiarray(myreal, *this);
        }
        multiarray<T, depth, Allocator> &imag()
        {
            return transform_multiarray(myimag, *this);
        }

        /// partial derivative in i_dim-th dimension
        multiarray<T, depth, Allocator> &partial_deriv(const  vec<double>& xs, const size_t i_dim)
        {
            return ::partial_deriv(elements, xs, m_length, i_dim);
        }

        double max_norm() {
            return elements.max_norm();
        }


        // non-member op's
        friend bool operator==<>(const multiarray<T, depth, Allocator> &, const multiarray<T, depth, Allocator> &);

        template <typename Q, typename L, typename R, size_t depthLR, typename AllocatorLR>
        friend auto elementwise(const Q &op, const multiarray<L, depthLR, AllocatorLR> &lhs, const multiarray<R, depthLR, AllocatorLR> &rhs);

        size_type size() const noexcept
        {
            return elements.size();
        }

        dimensions_type length() const noexcept
        {
            return m_length;
        }

        static constexpr size_type get_depth()
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

        template <typename Other_Alloc>
        multiarray<T, depth, Allocator> &do_copy_assign(const multiarray<T, depth, Other_Alloc> &other)
        {
            if (this != &other)
            {
                m_length = other.m_length;
                elements = other.elements;
                check_size();
            }
            return *this;
        }

        template <typename Other_Alloc>
        multiarray<T, depth, Allocator> &do_move_assign(multiarray<T, depth, Other_Alloc> &&other)
        {
            if (this != &other)
            {
                m_length = std::move(other.m_length);
                elements = std::move(other.elements);
                other.m_length = {0};
                check_size();
            }
            return *this;
        }
    };

    // === template definitions ===

    template <typename T, size_t depth, typename Allocator>
    typename multiarray<T, depth, Allocator>::size_type // <- This should be size_t. Compiler is too stupid to figure that out himself.
    multiarray<T, depth, Allocator>::flat_index(const index_type &index) const noexcept
    {
        //size_type res = index[0];
        //for (size_t i = 1; i < depth; i++)
        //{
        //    res *= m_length[i];
        //    res += index[i];
        //}

        return getFlatIndex<depth>(index, m_length);
    }

    template <typename T, size_t depth, typename Allocator>
    bool multiarray<T, depth, Allocator>::check_bounds(const index_type &index) const noexcept
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

    template <typename T, size_t depth, typename Allocator>
    bool multiarray<T, depth, Allocator>::check_bounds_end(const index_type &index) const noexcept
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

    template <typename T, size_t depth, typename Allocator>
    T &multiarray<T, depth, Allocator>::at(const index_type &index)
    {
        if (check_bounds(index))
        {
            return flat_at(flat_index(index));
        }
        else
        {
            throw std::out_of_range("Attempted to access multiarray element at invalid index.");
        }
    }

    template <typename T, size_t depth, typename Allocator>
    const T &multiarray<T, depth, Allocator>::at(const index_type &index) const
    {
        if (check_bounds(index))
        {
            return flat_at(flat_index(index));
        }
        else
        {
            throw std::out_of_range("Attempted to access multiarray element at invalid index.");
        }
    }

    // === operator implementation ===
    template <typename T, size_t depth, typename Allocator>
    bool operator==(const multiarray<T, depth, Allocator> &lhs, const multiarray<T, depth, Allocator> &rhs)
    {
        return lhs.elements == rhs.elements;
    }

    template <typename T, size_t depth, typename Allocator>
    bool operator!=(const multiarray<T, depth, Allocator> &lhs, const multiarray<T, depth, Allocator> &rhs)
    {
        return !(lhs == rhs);
    }

    // === arithmetic operators ===
    template <typename Q, typename L, typename R, size_t depth, typename Allocator>
    auto elementwise(const Q &op, const multiarray<L, depth, Allocator> &lhs, const multiarray<R, depth, Allocator> &rhs)
    {
        if (lhs.length() != rhs.length())
        {
            throw std::length_error("Cannot perform pairwise operations on multiarrays of different length.");
        }
        multiarray<decltype(op(std::declval<L>(), std::declval<R>())), depth> res(lhs.length());
        for (size_t i = 0; i < res.size(); i++)
        {
            res.elements[i] = op(lhs.elements[i], rhs.elements[i]);
        }
        return res;
    }

    template <typename Q, typename L, size_t depth, typename Allocator>
    auto transform_multiarray(const Q &op, const multiarray<L, depth, Allocator> &arr)
    {

        multiarray<decltype(op(std::declval<L>())), depth> res(arr.length());
        for (size_t i = 0; i < res.size(); i++)
        {
            res.elements[i] = op(arr.elements[i]);
        }
        return res;
    }

    template <typename Q, typename L, typename R, size_t depth, typename Allocator>
    auto scalar_map(const Q &op, const multiarray<L, depth, Allocator> &lhs, const R &rhs)
    {
        multiarray<decltype(op(std::declval<L>(), std::declval<R>())), depth> res(lhs.length());
        for (size_t i = 0; i < res.size(); i++)
        {
            res.elements[i] = op(lhs.elements[i], rhs);
        }
        return res;
    }

    /// Unary operators


    /// Binary operators
    /// Addition

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator+(const multiarray<L, depth, Allocator> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return l + r; }, lhs, rhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator+(const R &lhs, const multiarray<L, depth, Allocator> &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return r + l; }, rhs, lhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator+(const multiarray<L, depth, Allocator> &lhs, const multiarray<R, depth, Allocator> &rhs)
    {
        return elementwise([](const L &l, const R &r) { return l + r; }, lhs, rhs);
    }

    // Subtraction

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator-(const multiarray<L, depth, Allocator> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return l - r; }, lhs, rhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator-(const R &lhs, const multiarray<L, depth, Allocator> &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return r - l; }, rhs, lhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator-(const multiarray<L, depth, Allocator> &lhs, const multiarray<R, depth, Allocator> &rhs)
    {
        return elementwise([](const L &l, const R &r) { return l - r; }, lhs, rhs);
    }

    // Multiplication

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator*(const multiarray<L, depth, Allocator> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return l * r; }, lhs, rhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator*(const R &lhs, const multiarray<L, depth, Allocator> &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return r * l; }, rhs, lhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator*(const multiarray<L, depth, Allocator> &lhs, const multiarray<R, depth, Allocator> &rhs)
    {
        return elementwise([](const L &l, const R &r) { return l * r; }, lhs, rhs);
    }

    // Division

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator/(const multiarray<L, depth, Allocator> &lhs, const R &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return l / r; }, lhs, rhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator/(const R &lhs, const multiarray<L, depth, Allocator> &rhs)
    {
        return scalar_map([](const L &l, const R &r) { return r / l; }, rhs, lhs);
    }

    template <typename L, typename R, size_t depth, typename Allocator>
    auto operator/(const multiarray<L, depth, Allocator> &lhs, const multiarray<R, depth, Allocator> &rhs)
    {
        return elementwise([](const L &l, const R &r) { return l / r; }, lhs, rhs);
    }

} // namespace multidimensional

#endif

