#ifndef RANGED_VIEW_HPP
#define RANGED_VIEW_HPP

#include <stdexcept>
#include <string>
#include "../utilities/math_utils.hpp"
#include "../data_structures.hpp"
#include "multiarray.hpp"

namespace multidimensional
{

    template <class Iter>
    class ranged_view
    {
        Iter b;
        Iter e;

    public:
        ranged_view(Iter b, Iter e) : b(b), e(e) {}

        Iter begin() { return b; }
        Iter end() { return e; }

        auto &operator[] (size_t i)
        {
#if DEBUG
            if (b + i < e)
        {
#endif
            return b[i];
#if DEBUG
            }
        else
        {
            throw std::out_of_range("Attempted to access ranged_view beyond last element (at " + std::to_string(i) + ").");
        }
#endif
        }
    };

    template <class Container>
    ranged_view<typename Container::iterator>
    make_range(Container &c, size_t b, size_t e)
    {
        return ranged_view<typename Container::iterator>(c.begin() + b, c.begin() + e);
    }

    template <typename T, size_t depth>
    class BlockView {
        private:
            using index_type = std::array<std::size_t,depth>;
        public:
            BlockView(const multiarray<T,depth>& data, const index_type& start, const index_type length)
                : m_data(data), m_start(start), m_length(length), m_stride(data.length()) {
            // check lengths
            const index_type& length_base = data.length();
            bool is_sensible = true;
            for (int i = 0; i < depth; i++) {
                if (m_start[i] + m_length[i] > length_base[i]) is_sensible = false;
            }
            assert(is_sensible);
        };

        private:
            const index_type m_start;
            const index_type m_length;
            const index_type m_stride;
            const multidimensional::multiarray<T,depth>& m_data;




        bool check_bounds(const index_type &index) const
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

        size_t flat_index(const index_type &index) const noexcept
        {
            size_t res = m_start[0] + index[0];
            for (size_t i = 1; i < depth; i++)
            {
                res *= m_stride[i];
                res += m_start[i] + index[i];
            }

            return res;
        }

    public:
        const T &at(const index_type &index) const {
            if (check_bounds(index))
            {
                const size_t flat_idx = flat_index(index);
                return m_data.flat_at(flat_idx);
            }
            else
            {
                assert(false);
                throw std::out_of_range("Attempted to access multiarray element at invalid index.");
            }
        }


        vec<T> get_vec() const {
            const size_t flat_size = getFlatSize(m_length);
            vec<T> res(flat_size);
            for (size_t i = 0; i < flat_size; i++) {
                index_type idx;
                getMultIndex<depth>(idx, i, m_length);
                res[i] = at(idx);
            }
            return res;
        }


    };

} // namespace multi

#endif
