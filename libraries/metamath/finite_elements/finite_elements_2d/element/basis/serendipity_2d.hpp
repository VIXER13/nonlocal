#ifndef FINITE_ELEMENT_2D_BASIS_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_BASIS_SERENDIPITY_HPP

#include "symbolic_base.hpp"

#include "geometry_2d.hpp"
#include "geometric_primitives/rectangle.hpp"

namespace metamath::finite_element {

template<class T, size_t Order>
class serendipity;

template<class T>
class serendipity_base : public geometry_2d<T, rectangle_element_geometry> {
protected:
    static inline constexpr symbolic::variable<2> p{};
    T _p;

    explicit serendipity_base(const T p) : _p{p} {}
    ~serendipity_base() override = default;

public:
    T get_parameter() const noexcept { return _p; }
    void set_parameter(const T p) noexcept { _p = p; }
};

template<class T>
class serendipity<T, 0> : public geometry_2d<T, rectangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 1> nodes = { T{0}, T{0} };
    static inline constexpr auto basis = std::make_tuple(metamath::symbolic::integral_constant<1>{});

    explicit serendipity() = default;
    ~serendipity() override = default;
};

template<class T>
class serendipity<T, 1> : public geometry_2d<T, rectangle_element_geometry> {
    // Билинейный элемент

    using _base = geometry_2d<T, rectangle_element_geometry>;

    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<4> _4{};

protected:
    using _base::x;
    using _base::y;

    /*
        3---2
        |   |
        0---1
    */
    static inline constexpr std::array<std::array<T, 2>, 4>
        nodes = { T{-1}, T{-1},
                  T{ 1}, T{-1},
                  T{ 1}, T{ 1},
                  T{-1}, T{ 1} };

    static inline constexpr auto basis = std::make_tuple(
        (_1 - x) * (_1 - y) / _4,
        (_1 + x) * (_1 - y) / _4,
        (_1 + x) * (_1 + y) / _4,
        (_1 - x) * (_1 + y) / _4
    );

    explicit serendipity() = default;
    ~serendipity() override = default;
};

template<class T>
class serendipity<T, 2> : public serendipity_base<T> {
    // Астионенко И.А., Гучек П.И., Литвиненко Е.И., Хомченко А.Н. Применениче альтренативных серендиповых моделей при решении задачи о кручении призматических стержней.

    using _base = serendipity_base<T>;

    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<3> _3{};
    static inline constexpr metamath::symbolic::integral_constant<5> _5{};
    static inline constexpr metamath::symbolic::integral_constant<9> _9{};
    static inline constexpr metamath::symbolic::integral_constant<16> _16{};

protected:
    using _base::x;
    using _base::y;
    using _base::p;

    /*
        6---5---4
        |       |
        7       3
        |       |
        0---1---2
    */
    static inline constexpr std::array<std::array<T, 2>, 8>
        nodes = { T{-1}, T{-1},
                  T{ 0}, T{-1},
                  T{ 1}, T{-1},
                  T{ 1}, T{ 0},
                  T{ 1}, T{ 1},
                  T{ 0}, T{ 1},
                  T{-1}, T{ 1},
                  T{-1}, T{ 0} };

    static inline constexpr auto basis = std::make_tuple(
         (_1-x  ) * (_1-y) * ((_9*p-_1)*(_1+x+y) + (_9*p+_3)*x*y) / _16,
        -(_1-x*x) * (_1-y) * ((_9*p-_5)          + (_9*p+_3)*y  ) / _16,
         (_1+x  ) * (_1-y) * ((_9*p-_1)*(_1-x+y) - (_9*p+_3)*x*y) / _16,
        -(_1-y*y) * (_1+x) * ((_9*p-_5)          - (_9*p+_3)*x  ) / _16,
         (_1+x  ) * (_1+y) * ((_9*p-_1)*(_1-x-y) + (_9*p+_3)*x*y) / _16,
        -(_1-x*x) * (_1+y) * ((_9*p-_5)          - (_9*p+_3)*y  ) / _16,
         (_1-x  ) * (_1+y) * ((_9*p-_1)*(_1+x-y) - (_9*p+_3)*x*y) / _16,
        -(_1-y*y) * (_1-x) * ((_9*p-_5)          + (_9*p+_3)*x  ) / _16
    );

    explicit serendipity() : _base{T{2} / T{9}} {}
    ~serendipity() override = default;
};

template<class T>
class serendipity<T, 3> : public serendipity_base<T> {
    // Астионенко И.А., Литвиненко Е.И., Хомченко А.Н. Конструирование многопараметрических полиномов на бикубическом элементе серендипова семейства.

    using _base = serendipity_base<T>;

    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<2> _2{};
    static inline constexpr metamath::symbolic::integral_constant<9> _9{};
    static inline constexpr metamath::symbolic::integral_constant<18> _18{};
    static inline constexpr metamath::symbolic::integral_constant<32> _32{};
    static inline constexpr metamath::symbolic::integral_constant<54> _54{};
    static inline constexpr metamath::symbolic::integral_constant<64> _64{};

protected:
    using _base::x;
    using _base::y;
    using _base::p;

    /*
        9---8---7---6
        |           |
        10          5
        |           |
        11          4
        |           |
        0---1---2---3
    */
    static inline constexpr std::array<std::array<T, 2>, 12>
        nodes = {      -T{1},      -T{1},
                  -T{1}/T{3},      -T{1},
                   T{1}/T{3},      -T{1},
                        T{1},      -T{1},
                        T{1}, -T{1}/T{3},
                        T{1},  T{1}/T{3},
                        T{1},       T{1},
                   T{1}/T{3},       T{1},
                  -T{1}/T{3},       T{1},
                       -T{1},       T{1},
                       -T{1},  T{1}/T{3},
                       -T{1}, -T{1}/T{3}
        };

    static inline constexpr auto basis = std::make_tuple(
         (_1-x  ) * (_1-y  ) * (_9 *(x*x+y*y + (_2 *p+_1)*(x*y+x+y)) + _18*p - _1) / _32,
        -(_1-x*x) * (_1-y  ) * (_54*x        + (_18*p+_9)*y          + _18*p - _9) / _64,
         (_1-x*x) * (_1-y  ) * (_54*x        - (_18*p+_9)*y          - _18*p + _9) / _64,
         (_1+x  ) * (_1-y  ) * (_9 *(x*x+y*y - (_2 *p+_1)*(x*y+x-y)) + _18*p - _1) / _32,
        -(_1+x  ) * (_1-y*y) * (_54*y        - (_18*p+_9)*x          + _18*p - _9) / _64,
         (_1+x  ) * (_1-y*y) * (_54*y        + (_18*p+_9)*x          - _18*p + _9) / _64,
         (_1+x  ) * (_1+y  ) * (_9 *(x*x+y*y + (_2 *p+_1)*(x*y-x-y)) + _18*p - _1) / _32,
         (_1-x*x) * (_1+y  ) * (_54*x        + (_18*p+_9)*y          - _18*p + _9) / _64,
        -(_1-x*x) * (_1+y  ) * (_54*x        - (_18*p+_9)*y          + _18*p - _9) / _64,
         (_1-x  ) * (_1+y  ) * (_9 *(x*x+y*y - (_2 *p+_1)*(x*y-x+y)) + _18*p - _1) / _32,
         (_1-x  ) * (_1-y*y) * (_54*y        - (_18*p+_9)*x          - _18*p + _9) / _64,
        -(_1-x  ) * (_1-y*y) * (_54*y        + (_18*p+_9)*x          + _18*p - _9) / _64
    );

    explicit serendipity() : _base{T{1} / T{8}} {}
    ~serendipity() override = default;
};

template<class T>
class serendipity<T, 4> : public geometry_2d<T, rectangle_element_geometry> {
    // Астионенко И.А., Козуб Н.А., Литвиненко Е.И., Хомченко А.Н. Управляемые серендиповы поверхности, сохраняющие межэлементную непрерывность.

    using _base = geometry_2d<T, rectangle_element_geometry>;

    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<2> _2{};
    static inline constexpr metamath::symbolic::integral_constant<3> _3{};
    static inline constexpr metamath::symbolic::integral_constant<61> _61{};
    static inline constexpr metamath::symbolic::integral_constant<141> _141{};
    static inline constexpr metamath::symbolic::integral_constant<203> _203{};
    static inline constexpr metamath::symbolic::integral_constant<311> _311{};
    static inline constexpr metamath::symbolic::integral_constant<500> _500{};
    static inline constexpr metamath::symbolic::integral_constant<561> _561{};
    static inline constexpr metamath::symbolic::integral_constant<800> _800{};
    static inline constexpr metamath::symbolic::integral_constant<1141> _1141{};
    static inline constexpr metamath::symbolic::integral_constant<1600> _1600{};
    static inline constexpr metamath::symbolic::integral_constant<2000> _2000{};
    static inline constexpr metamath::symbolic::integral_constant<3000> _3000{};
    static inline constexpr metamath::symbolic::integral_constant<4000> _4000{};

protected:
    using _base::x;
    using _base::y;

    /*
        12--11--10--9--8
        |              |
        13             7
        |              |
        14             6
        |              |
        15             5
        |              |
        0---1---2---3--4
    */
    static inline constexpr std::array<std::array<T, 2>, 16>
        nodes = { T{-1.0}, T{-1.0},
                  T{-0.5}, T{-1.0},
                  T{ 0.0}, T{-1.0},
                  T{ 0.5}, T{-1.0},
                  T{ 1.0}, T{-1.0},
                  T{ 1.0}, T{-0.5},
                  T{ 1.0}, T{ 0.0},
                  T{ 1.0}, T{ 0.5},
                  T{ 1.0}, T{ 1.0},
                  T{ 0.5}, T{ 1.0},
                  T{ 0.0}, T{ 1.0},
                  T{-0.5}, T{ 1.0},
                  T{-1.0}, T{ 1.0},
                  T{-1.0}, T{ 0.5},
                  T{-1.0}, T{ 0.0},
                  T{-1.0}, T{-0.5}
        };

    static inline constexpr auto basis = std::make_tuple(
         (_1-x  ) * (_1-y  ) * (_1 + _2*(x + y)) * (_561 + _61*(x+y) - _500*(x*x + y*y) + _311*x*y) / _3000,
        -(_1-x*x) * (_1-y  ) * (_203  + _203*y + _1600*x/_3) * (_1-_2*x) / _800,
         (_1-x*x) * (_1-y  ) * (_1141 + _141*y - _4000*x*x) / _2000,
        -(_1-x*x) * (_1-y  ) * (_203  + _203*y - _1600*x/_3) * (_1+_2*x) / _800,
         (_1+x  ) * (_1-y  ) * (_1 - _2*(x - y)) * (_561 - _61*(x-y) - _500*(x*x + y*y) - _311*x*y) / _3000,
        -(_1+x  ) * (_1-y*y) * (_203  - _203*x + _1600*y/_3) * (_1-_2*y) / _800,
         (_1+x  ) * (_1-y*y) * (_1141 - _141*x - _4000*y*y) / _2000,
        -(_1+x  ) * (_1-y*y) * (_203  - _203*x - _1600*y/_3) * (_1+_2*y) / _800,
         (_1+x  ) * (_1+y  ) * (_1 - _2*(x + y)) * (_561 - _61*(x+y) - _500*(x*x + y*y) + _311*x*y) / _3000,
        -(_1-x*x) * (_1+y  ) * (_203  - _203*y - _1600*x/_3) * (_1+_2*x) / _800,
         (_1-x*x) * (_1+y  ) * (_1141 - _141*y - _4000*x*x) / _2000,
        -(_1-x*x) * (_1+y  ) * (_203  - _203*y + _1600*x/_3) * (_1-_2*x) / _800,
         (_1-x  ) * (_1+y  ) * (_1 + _2*(x - y)) * (_561 + _61*(x-y) - _500*(x*x + y*y) - _311*x*y) / _3000,
        -(_1-x  ) * (_1-y*y) * (_203  + _203*x - _1600*y/_3) * (_1+_2*y) / _800,
         (_1-x  ) * (_1-y*y) * (_1141 + _141*x - _4000*y*y) / _2000,
        -(_1-x  ) * (_1-y*y) * (_203  + _203*x + _1600*y/_3) * (_1-_2*y) / _800
    );

    explicit serendipity() = default;
    ~serendipity() override = default;
};

template<class T>
class serendipity<T, 5> : public geometry_2d<T, rectangle_element_geometry> {
    // Хомченко А.Н., Астионенко И.А. Серендиповы поверхности высших порядков: особенности формообразования.

    using _base = geometry_2d<T, rectangle_element_geometry>;

    static inline constexpr metamath::symbolic::integral_constant<1> _1{};
    static inline constexpr metamath::symbolic::integral_constant<3> _3{};
    static inline constexpr metamath::symbolic::integral_constant<5> _5{};
    static inline constexpr metamath::symbolic::integral_constant<9> _9{};
    static inline constexpr metamath::symbolic::integral_constant<25> _25{};
    static inline constexpr metamath::symbolic::integral_constant<125> _125{};
    static inline constexpr metamath::symbolic::integral_constant<384> _384{};
    static inline constexpr metamath::symbolic::integral_constant<768> _768{};
    static inline constexpr metamath::symbolic::integral_constant<1536> _1536{};

protected:
    using _base::x;
    using _base::y;

    /*
        15--14--13--12--11-10
        |                   |
        16                  9
        |                   |
        17                  8
        |                   |
        18                  7
        |                   |
        19                  6
        |                   |
        0---1---2---3---4---5
    */
    static inline constexpr std::array<std::array<T, 2>, 20>
        nodes = { T{-1.0}, T{-1.0},
                  T{-0.6}, T{-1.0},
                  T{-0.2}, T{-1.0},
                  T{ 0.2}, T{-1.0},
                  T{ 0.6}, T{-1.0},
                  T{ 1.0}, T{-1.0},
                  T{ 1.0}, T{-0.6},
                  T{ 1.0}, T{-0.2},
                  T{ 1.0}, T{ 0.2},
                  T{ 1.0}, T{ 0.6},
                  T{ 1.0}, T{ 1.0},
                  T{ 0.6}, T{ 1.0},
                  T{ 0.2}, T{ 1.0},
                  T{-0.2}, T{ 1.0},
                  T{-0.6}, T{ 1.0},
                  T{-1.0}, T{ 1.0},
                  T{-1.0}, T{ 0.6},
                  T{-1.0}, T{ 0.2},
                  T{-1.0}, T{-0.2},
                  T{-1.0}, T{-0.6}
        };

    static inline constexpr auto basis = std::make_tuple(
              (_1-x  ) * (_1-y) * (_384 - _125 * ((_1-x*x) * (_3+_5*x*x) + (_1-y*y)*(_3+_5*y*y))) / _1536,
        _25 * (_1-x*x) * (_1-y) * (-_1+_25*x*x) * (_3-_5*x) / _1536,
        _25 * (_1-x*x) * (_1-y) * ( _9-_25*x*x) * (_1-_5*x) / _768,
        _25 * (_1-x*x) * (_1-y) * ( _9-_25*x*x) * (_1+_5*x) / _768,
        _25 * (_1-x*x) * (_1-y) * (-_1+_25*x*x) * (_3+_5*x) / _1536,
              (_1+x  ) * (_1-y) * (_384 - _125 * ((_1-x*x) * (_3+_5*x*x) + (_1-y*y)*(_3+_5*y*y))) / _1536,
        _25 * (_1-y*y) * (_1+x) * (-_1+_25*y*y) * (_3-_5*y) / _1536,
        _25 * (_1-y*y) * (_1+x) * ( _9-_25*y*y) * (_1-_5*y) / _768,
        _25 * (_1-y*y) * (_1+x) * ( _9-_25*y*y) * (_1+_5*y) / _768,
        _25 * (_1-y*y) * (_1+x) * (-_1+_25*y*y) * (_3+_5*y) / _1536,
              (_1+x  ) * (_1+y) * (_384 - _125 * ((_1-x*x) * (_3+_5*x*x) + (_1-y*y)*(_3+_5*y*y))) / _1536,
        _25 * (_1-x*x) * (_1+y) * (-_1+_25*x*x) * (_3+_5*x) / _1536,
        _25 * (_1-x*x) * (_1+y) * ( _9-_25*x*x) * (_1+_5*x) / _768,
        _25 * (_1-x*x) * (_1+y) * ( _9-_25*x*x) * (_1-_5*x) / _768,
        _25 * (_1-x*x) * (_1+y) * (-_1+_25*x*x) * (_3-_5*x) / _1536,
              (_1-x  ) * (_1+y) * (_384 - _125 * ((_1-x*x) * (_3+_5*x*x) + (_1-y*y)*(_3+_5*y*y))) / _1536,
        _25 * (_1-y*y) * (_1-x) * (-_1+_25*y*y) * (_3+_5*y) / _1536,
        _25 * (_1-y*y) * (_1-x) * ( _9-_25*y*y) * (_1+_5*y) / _768,
        _25 * (_1-y*y) * (_1-x) * ( _9-_25*y*y) * (_1-_5*y) / _768,
        _25 * (_1-y*y) * (_1-x) * (-_1+_25*y*y) * (_3-_5*y) / _1536
    );

    explicit serendipity() = default;
    ~serendipity() override = default;
};

}

#endif