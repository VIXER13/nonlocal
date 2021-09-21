#ifndef FINITE_ELEMENT_1D_BASIS_QUBIC_ELEMENT_HPP
#define FINITE_ELEMENT_1D_BASIS_QUBIC_ELEMENT_HPP

#include "geometry_1d.hpp"
#include <functional>

namespace metamath::finite_element {

template<class T>
class qubic : protected geometry_1d<T, standart_segment_geometry> {
protected:
    using geometry_1d<T, standart_segment_geometry>::xi;

    explicit qubic() noexcept = default;
    ~qubic() override = default;

    // Нумерация узлов на кубическом элементе: 0--1--2--3
    static constexpr std::array<T, 4> nodes = { T{-1}, T{-1}/T{3}, T{1}/T{3}, T{1} };

    static constexpr auto basis = std::make_tuple(
        T{-0.5625} * (xi -      T{1}) * (xi * xi - T{1}/T{9}),
        T{ 1.6875} * (xi - T{1}/T{3}) * (xi * xi -      T{1}),
        T{-1.6875} * (xi + T{1}/T{3}) * (xi * xi -      T{1}),
        T{ 0.5625} * (xi +      T{1}) * (xi * xi - T{1}/T{9})
    );
};

}

#endif