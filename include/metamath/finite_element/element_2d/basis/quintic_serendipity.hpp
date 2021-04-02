#ifndef FINITE_ELEMENT_2D_BASIS_QUINTIC_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_BASIS_QUINTIC_SERENDIPITY_HPP

// Базисные функции элемента взяты из статьи:
// Хомченко А.Н., Астионенко И.А. Серендиповы поверхности высших порядков: особенности формообразования.

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class quintic_serendipity : public geometry_2d<T, rectangle_element_geometry> {
protected:
    using geometry_2d<T, rectangle_element_geometry>::xi;
    using geometry_2d<T, rectangle_element_geometry>::eta;

    explicit quintic_serendipity() = default;
    ~quintic_serendipity() override = default;

    // Нумерация узлов на пентичном серендиповом элементе:  15--14--13--12--11-10
    //                                                      |                   |
    //                                                      16                  9
    //                                                      |                   |
    //                                                      17                  8
    //                                                      |                   |
    //                                                      18                  7
    //                                                      |                   |
    //                                                      19                  6
    //                                                      |                   |
    //                                                      0---1---2---3---4---5
    static constexpr std::array<std::array<T, 2>, 20> nodes = { T{-1.0}, T{-1.0},
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
                                                                T{-1.0}, T{-0.6} };

    static constexpr auto basis = std::make_tuple(
        T{ 1}/T{1536} * (T{1}-xi     ) * (T{1}-eta) * (T{384} - T{125} * (T{1}-xi*xi) * (T{3}+T{5}*xi*xi) - T{125}*(T{1}-eta*eta)*(T{3}+T{5}*eta*eta)),
        T{25}/T{1536} * (T{1}-xi*xi  ) * (T{1}-eta) * (T{-1}+T{25}*xi*xi) * (T{3}-T{5}*xi),
        T{25}/T{ 768} * (T{1}-xi*xi  ) * (T{1}-eta) * (T{ 9}-T{25}*xi*xi) * (T{1}-T{5}*xi),
        T{25}/T{ 768} * (T{1}-xi*xi  ) * (T{1}-eta) * (T{ 9}-T{25}*xi*xi) * (T{1}+T{5}*xi),
        T{25}/T{1536} * (T{1}-xi*xi  ) * (T{1}-eta) * (T{-1}+T{25}*xi*xi) * (T{3}+T{5}*xi),
        T{ 1}/T{1536} * (T{1}+xi     ) * (T{1}-eta) * (T{384} - T{125} * (T{1}-xi*xi) * (T{3}+T{5}*xi*xi) - T{125}*(T{1}-eta*eta)*(T{3}+T{5}*eta*eta)),
        T{25}/T{1536} * (T{1}-eta*eta) * (T{1}+xi ) * (T{-1}+T{25}*eta*eta) * (T{3}-T{5}*eta),
        T{25}/T{ 768} * (T{1}-eta*eta) * (T{1}+xi ) * (T{ 9}-T{25}*eta*eta) * (T{1}-T{5}*eta),
        T{25}/T{ 768} * (T{1}-eta*eta) * (T{1}+xi ) * (T{ 9}-T{25}*eta*eta) * (T{1}+T{5}*eta),
        T{25}/T{1536} * (T{1}-eta*eta) * (T{1}+xi ) * (T{-1}+T{25}*eta*eta) * (T{3}+T{5}*eta),
        T{ 1}/T{1536} * (T{1}+xi     ) * (T{1}+eta) * (T{384} - T{125} * (T{1}-xi*xi) * (T{3}+T{5}*xi*xi) - T{125}*(T{1}-eta*eta)*(T{3}+T{5}*eta*eta)),
        T{25}/T{1536} * (T{1}-xi*xi  ) * (T{1}+eta) * (T{-1}+T{25}*xi*xi) * (T{3}+T{5}*xi),
        T{25}/T{ 768} * (T{1}-xi*xi  ) * (T{1}+eta) * (T{ 9}-T{25}*xi*xi) * (T{1}+T{5}*xi),
        T{25}/T{ 768} * (T{1}-xi*xi  ) * (T{1}+eta) * (T{ 9}-T{25}*xi*xi) * (T{1}-T{5}*xi),
        T{25}/T{1536} * (T{1}-xi*xi  ) * (T{1}+eta) * (T{-1}+T{25}*xi*xi) * (T{3}-T{5}*xi),
        T{ 1}/T{1536} * (T{1}-xi     ) * (T{1}+eta) * (T{384} - T{125} * (T{1}-xi*xi) * (T{3}+T{5}*xi*xi) - T{125}*(T{1}-eta*eta)*(T{3}+T{5}*eta*eta)),
        T{25}/T{1536} * (T{1}-eta*eta) * (T{1}-xi ) * (T{-1}+T{25}*eta*eta) * (T{3}+T{5}*eta),
        T{25}/T{ 768} * (T{1}-eta*eta) * (T{1}-xi ) * (T{ 9}-T{25}*eta*eta) * (T{1}+T{5}*eta),
        T{25}/T{ 768} * (T{1}-eta*eta) * (T{1}-xi ) * (T{ 9}-T{25}*eta*eta) * (T{1}-T{5}*eta),
        T{25}/T{1536} * (T{1}-eta*eta) * (T{1}-xi ) * (T{-1}+T{25}*eta*eta) * (T{3}-T{5}*eta)
    );
};

}

#endif