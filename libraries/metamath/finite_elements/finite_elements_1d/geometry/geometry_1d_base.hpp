#ifndef FINITE_ELEMENT_GEOMETRY_1D_BASE_HPP
#define FINITE_ELEMENT_GEOMETRY_1D_BASE_HPP

#include "side_1d.hpp"

namespace metamath::finite_element {

template<class T>
class geometry_1d_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    static inline constexpr symbolic::variable<0> x{};

public:
    virtual ~geometry_1d_base() noexcept = default;

    virtual T boundary(const side_1d bound) const = 0;
};

}

#endif