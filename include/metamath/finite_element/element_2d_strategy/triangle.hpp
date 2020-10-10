#ifndef FINITE_ELEMENT_TRIANGLE_ELEMENT_HPP
#define FINITE_ELEMENT_TRIANGLE_ELEMENT_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class T>
class triangle : public barycentric<T> {
protected:
    using barycentric<T>::xi;
    using barycentric<T>::eta;
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    explicit triangle() = default;
    ~triangle() override = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          2--0
    */
    static constexpr std::array<std::array<T, 2>, 3> nodes = { 1., 0.,
                                                               0., 1.,
                                                               0., 0. };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i, i = 1..3
    static constexpr auto basis = std::make_tuple(L1, L2, L3);
};

}

#endif