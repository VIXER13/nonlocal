#pragma once

#include <metamath/symbolic/base/variable.hpp>

#include <array>

namespace metamath::finite_element {

template<class T>
class standart_segment_geometry {
protected:
    static inline constexpr std::array<T, 2> boundary = { T{-1}, T{1} };

    constexpr explicit standart_segment_geometry() noexcept = default;

public:
    virtual ~standart_segment_geometry() noexcept = default;
};

}