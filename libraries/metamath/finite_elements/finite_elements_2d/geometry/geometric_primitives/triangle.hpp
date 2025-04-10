#pragma once

#include <array>
#include <functional>

namespace metamath::finite_element {

template<class T>
class triangle_element_geometry {
protected:
    static inline const std::array<std::function<T(const T)>, 4>
        boundary = { [](const T x) constexpr noexcept { return T{0};     },
                     [](const T x) constexpr noexcept { return T{1};     },
                     [](const T y) constexpr noexcept { return T{0};     },
                     [](const T y) constexpr noexcept { return T{1} - y; } };

    explicit triangle_element_geometry() = default;

public:
    virtual ~triangle_element_geometry() noexcept = default;
};

}