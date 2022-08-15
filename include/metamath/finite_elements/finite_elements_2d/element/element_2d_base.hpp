#ifndef FINITE_ELEMENT_2D_BASE_HPP
#define FINITE_ELEMENT_2D_BASE_HPP

#include "finite_element_base.hpp"
#include "side_2d.hpp"

#include <array>

namespace metamath::finite_element {

template<class T>
class element_2d_base : public element_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

public:
    ~element_2d_base() override = default;

    virtual const std::array<T, 2>& node(const size_t i) const = 0;

    virtual T N   (const size_t i, const std::array<T, 2>& x) const = 0;
    virtual T Nxi (const size_t i, const std::array<T, 2>& x) const = 0;
    virtual T Neta(const size_t i, const std::array<T, 2>& x) const = 0;

    virtual T boundary(const side_2d bound, const T x) const = 0;
};

}

#endif