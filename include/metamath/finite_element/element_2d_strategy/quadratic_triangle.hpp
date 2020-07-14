#ifndef FINITE_ELEMENT_QUADRATIC_TRIANGLE_ELEMENT_HPP
#define FINITE_ELEMENT_QUADRATIC_TRIANGLE_ELEMENT_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class T>
class quadratic_triangle : public barycentric<T> {
protected:
    using barycentric<T>::xi;
    using barycentric<T>::eta;
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    explicit quadratic_triangle() = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          4  3
                                                          |    \
                                                          2--5--0
    */
    static constexpr std::array<std::array<T, 2>, 6> nodes = { 1.0, 0.0,
                                                               0.0, 1.0,
                                                               0.0, 0.0,
                                                               0.5, 0.5,
                                                               0.0, 0.5,
                                                               0.5, 0.0 };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i (2 L_i - 1), i = 1..3,
    //                                                            N_3 = 4 L_1 L_2,
    //                                                            N_4 = 4 L_2 L_3,
    //                                                            N_5 = 4 L_3 L_1.
    static constexpr auto basis = std::make_tuple(
        L1 * (2. * L1 - 1),
        L2 * (2. * L2 - 1),
        L3 * (2. * L3 - 1),
        4. *  L1 * L2,
        4. *  L2 * L3,
        4. *  L3 * L1
    );

    static inline const std::array<std::function<T(const std::array<T, 2>&)>, 6>
        N    = symdiff::to_function<T, 2>(basis),
        Nxi  = symdiff::to_function<T, 2>(symdiff::derivative<xi>(basis)),
        Neta = symdiff::to_function<T, 2>(symdiff::derivative<eta>(basis));
};

}

#endif