#pragma once

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class T, size_t Order>
class triangle;

template<class T>
class triangle<T, 0> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 1> nodes = { T{1} / T{3}, T{1} / T{3} };
    static inline constexpr auto basis = std::make_tuple(metamath::symbolic::integral_constant<1>{});

    explicit triangle() = default;
    ~triangle() override = default;
};

template<class T>
class triangle<T, 1> : public barycentric<T> {
protected:
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    /*
        1\
        | \
        2--0
    */
    static inline constexpr std::array<std::array<T, 2>, 3>
        nodes = { T{1}, T{0},
                  T{0}, T{1},
                  T{0}, T{0} };

    static inline constexpr auto basis = std::make_tuple(L1, L2, L3);

    explicit triangle() = default;
    ~triangle() override = default;
};

template<class T>
class triangle<T, 2> : public barycentric<T> {
    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<2> _2{};
    static inline constexpr metamath::symbolic::integral_constant<4> _4{};

protected:
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    /*
        1\
        | \
        4  3
        |    \
        2--5--0
    */
    static inline constexpr std::array<std::array<T, 2>, 6>
        nodes = { T{1.0}, T{0.0},
                  T{0.0}, T{1.0},
                  T{0.0}, T{0.0},
                  T{0.5}, T{0.5},
                  T{0.0}, T{0.5},
                  T{0.5}, T{0.0} };

    static inline constexpr auto basis = std::make_tuple(
        L1 * (_2 * L1 - _1),
        L2 * (_2 * L2 - _1),
        L3 * (_2 * L3 - _1),
        _4 *  L1 * L2,
        _4 *  L2 * L3,
        _4 *  L3 * L1
    );

    explicit triangle() = default;
    ~triangle() override = default;
};

template<class T>
class triangle<T, 3> : public barycentric<T> {
    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<2> _2{};
    static inline constexpr metamath::symbolic::integral_constant<3> _3{};
    static inline constexpr metamath::symbolic::integral_constant<9> _9{};
    static inline constexpr metamath::symbolic::integral_constant<27> _27{};

protected:
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    /*
        1\
        5 4
        |  \
        6 9 3
        |    \
        2-7-8-0
    */
    static inline constexpr std::array<std::array<T, 2>, 10>
        nodes = {        T{1},        T{0},
                         T{0},        T{1},
                         T{0},        T{0},
                  T{2} / T{3}, T{1} / T{3},
                  T{1} / T{3}, T{2} / T{3},
                         T{0}, T{2} / T{3},
                         T{0}, T{1} / T{3},
                  T{1} / T{3},        T{0},
                  T{2} / T{3},        T{0},
                  T{1} / T{3}, T{1} / T{3} };

    static inline constexpr auto basis = std::make_tuple(
        L1 * (_3 * L1 - _1) * (_3 * L1 - _2) / _2,
        L2 * (_3 * L2 - _1) * (_3 * L2 - _2) / _2,
        L3 * (_3 * L3 - _1) * (_3 * L3 - _2) / _2,
        _9  * L1 * L2       * (_3 * L1 - _1) / _2,
        _9  * L1 * L2       * (_3 * L2 - _1) / _2,
        _9  * L2 * L3       * (_3 * L2 - _1) / _2,
        _9  * L2 * L3       * (_3 * L3 - _1) / _2,
        _9  * L3 * L1       * (_3 * L3 - _1) / _2,
        _9  * L3 * L1       * (_3 * L1 - _1) / _2,
        _27 * L1 * L2       *       L3
    );

    explicit triangle() = default;
    ~triangle() override = default;
};

}