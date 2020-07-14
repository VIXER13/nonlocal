#ifndef FINITE_ELEMENT_QUBIC_TRIANGLE_ELEMENT_HPP
#define FINITE_ELEMENT_QUBIC_TRIANGLE_ELEMENT_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class T>
class qubic_triangle : public barycentric<T> {
protected:
    using barycentric<T>::xi;
    using barycentric<T>::eta;
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    explicit qubic_triangle() = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          5 4
                                                          |  \
                                                          6 9 3
                                                          |    \
                                                          2-7-8-0
    */
    static constexpr std::array<std::array<T, 2>, 10> nodes = {    1.,    0.,
                                                                   0.,    1.,
                                                                   0.,    0.,
                                                                2./3., 1./3.,
                                                                1./3., 2./3.,
                                                                   0., 2./3.,
                                                                   0., 1./3.,
                                                                1./3.,    0.,
                                                                2./3.,    0.,
                                                                1./3., 1./3. };

    // Базисные функции в барицентрических координатах имеют вид: N_i = 0.5 L_i (3 L_i - 1)(3 L_i - 2), i = 1..3,
    //                                                            N_3 = 4.5 L_1 L_2 (3 L_1 - 1),
    //                                                            N_4 = 4.5 L_1 L_2 (3 L_2 - 1),
    //                                                            N_5 = 4.5 L_2 L_3 (3 L_2 - 1),
    //                                                            N_6 = 4.5 L_2 L_3 (3 L_3 - 1),
    //                                                            N_7 = 4.5 L_3 L_1 (3 L_3 - 1),
    //                                                            N_8 = 4.5 L_3 L_1 (3 L_1 - 1),
    //                                                            N_9 =  27 L_1 L_2 L_3
    static constexpr auto basis = std::make_tuple(
        0.5 * L1 * (3.*L1 - 1.) * (3.*L1 - 2.),
        0.5 * L2 * (3.*L2 - 1.) * (3.*L2 - 2.),
        0.5 * L3 * (3.*L3 - 1.) * (3.*L3 - 2.),
        4.5 * L1 * L2 * (3.*L1 - 1.),
        4.5 * L1 * L2 * (3.*L2 - 1.),
        4.5 * L2 * L3 * (3.*L2 - 1.),
        4.5 * L2 * L3 * (3.*L3 - 1.),
        4.5 * L3 * L1 * (3.*L3 - 1.),
        4.5 * L3 * L1 * (3.*L1 - 1.),
        27. * L1 * L2 * L3
    );

    static inline const std::array<std::function<T(const std::array<T, 2>&)>, 10>
        N    = symdiff::to_function<T, 2>(basis),
        Nxi  = symdiff::to_function<T, 2>(symdiff::derivative<xi>(basis)),
        Neta = symdiff::to_function<T, 2>(symdiff::derivative<eta>(basis));
};

}

#endif