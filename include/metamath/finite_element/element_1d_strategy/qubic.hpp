#ifndef FINITE_ELEMENT_QUBIC_ELEMENT_HPP
#define FINITE_ELEMENT_QUBIC_ELEMENT_HPP

#include <functional>
#include "geometry_1d.hpp"

namespace metamath::finite_element {

template<class T>
class qubic : protected geometry_1d<T, standart_segment_geometry> {
protected:
    using geometry_1d<T, standart_segment_geometry>::xi;

    explicit qubic() noexcept = default;

    // Нумерация узлов на кубическом элементе: 0--1--2--3
    static constexpr std::array<T, 4> nodes = { -1., -1./3., 1./3., 1. };

    static constexpr auto basis = std::make_tuple(
        -0.5625 * (xi -    1.) * (xi*xi - 1./9.),
         1.6875 * (xi - 1./3.) * (xi*xi -    1.),
        -1.6875 * (xi + 1./3.) * (xi*xi -    1.),
         0.5625 * (xi +    1.) * (xi*xi - 1./9.)
    );

    static inline const std::array<std::function<T(const std::array<T, 1>&)>, 4>
        N   = symdiff::to_function<T, 1>(basis),
        Nxi = symdiff::to_function<T, 1>(symdiff::derivative<xi>(basis));
};

}

#endif