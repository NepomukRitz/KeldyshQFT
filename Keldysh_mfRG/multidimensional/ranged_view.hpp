#ifndef RANGED_VIEW_HPP
#define RANGED_VIEW_HPP

#include <stdexcept>
#include <string>

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

} // namespace multi

#endif
