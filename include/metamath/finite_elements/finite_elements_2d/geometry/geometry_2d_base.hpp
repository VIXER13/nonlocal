#ifndef FINITE_ELEMENT_GEOMETRY_2D_BASE_HPP
#define FINITE_ELEMENT_GEOMETRY_2D_BASE_HPP

#include "variable.hpp"
#include "side_2d.hpp"

#include <type_traits>

namespace metamath::finite_element {

template<class T>
class geometry_2d_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    static inline constexpr symbolic::variable<0> x{};
    static inline constexpr symbolic::variable<1> y{};

public:
    virtual ~geometry_2d_base() noexcept = default;
    virtual T boundary(const side_2d bound, const T x) const = 0;
};

}

#endif