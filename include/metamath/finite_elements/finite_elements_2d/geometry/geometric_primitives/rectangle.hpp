#ifndef FINITE_ELEMENT_RECTANGLE_HPP
#define FINITE_ELEMENT_RECTANGLE_HPP

#include <array>
#include <functional>

namespace metamath::finite_element {

template<class T>
class rectangle_element_geometry {
protected:
    static inline const std::array<std::function<T(const T)>, 4>
        boundary = { [](const T x) constexpr noexcept { return T{-1}; },
                     [](const T x) constexpr noexcept { return T{ 1}; },
                     [](const T y) constexpr noexcept { return T{-1}; },
                     [](const T y) constexpr noexcept { return T{ 1}; } };

    explicit rectangle_element_geometry() = default;

public:
    virtual ~rectangle_element_geometry() noexcept = default;
};

}


#endif