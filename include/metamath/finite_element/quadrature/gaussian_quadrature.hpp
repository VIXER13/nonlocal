#ifndef GAUSSIAN_QUADRATURE_HPP
#define GAUSSIAN_QUADRATURE_HPP

#include <cmath>
#include "geometry_1d.hpp"

namespace metamath::finite_element {

// Наследование квадратур от класса геометрии подразумевает возможность использования нестандартных квадратур,
// а так же многомерных квадратур, которые не получаются путём декартова произведения одномерных квадратур.

template<class T>
class gauss1 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss1() noexcept = default;
    ~gauss1() override = default;

    static constexpr std::array<std::array<T, 1>, 1> nodes = { T{0} };
    static constexpr std::array<T, 1> weights = { T{2} };
};

template<class T>
class gauss2 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss2() noexcept = default;
    ~gauss2() override = default;

    static const inline std::array<std::array<T, 1>, 2> nodes = { T{-1} / std::sqrt(T{3}), T{1} / std::sqrt(T{3}) };
    static constexpr std::array<T, 2> weights = { T{1}, T{1} };
};

template<class T>
class gauss3 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss3() noexcept = default;
    ~gauss3() override = default;

    static const inline std::array<std::array<T, 1>, 3> nodes = { -std::sqrt(T{3}/T{5}), T{0}, std::sqrt(T{3}/T{5}) };
    static constexpr std::array<T, 3> weights = { T{5}/T{9}, T{8}/T{9}, T{5}/T{9} };
};

template<class T>
class gauss4 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss4() noexcept = default;
    ~gauss4() override = default;

    static const inline std::array<std::array<T, 1>, 4>
            nodes = { -std::sqrt(T{3}/T{7} + T{2}/T{7} * std::sqrt(T{6}/T{5})),
                      -std::sqrt(T{3}/T{7} - T{2}/T{7} * std::sqrt(T{6}/T{5})),
                       std::sqrt(T{3}/T{7} - T{2}/T{7} * std::sqrt(T{6}/T{5})),
                       std::sqrt(T{3}/T{7} + T{2}/T{7} * std::sqrt(T{6}/T{5})) };
    static const inline std::array<T, 4>
            weights = { (T{18} - std::sqrt(T{30})) / T{36},
                        (T{18} + std::sqrt(T{30})) / T{36},
                        (T{18} + std::sqrt(T{30})) / T{36},
                        (T{18} - std::sqrt(T{30})) / T{36} };
};

template<class T>
class gauss5 : public geometry_1d<T, standart_segment_geometry> {
protected:
    explicit gauss5() noexcept = default;
    ~gauss5() override = default;

    static const inline std::array<std::array<T, 1>, 5>
            nodes = { T{-1}/T{3} * std::sqrt(T{5} + T{2} * std::sqrt(T{10}/T{7})),
                      T{-1}/T{3} * std::sqrt(T{5} - T{2} * std::sqrt(T{10}/T{7})),
                      T{ 0},
                      T{ 1}/T{3} * std::sqrt(T{5} - T{2} * std::sqrt(T{10}/T{7})),
                      T{ 1}/T{3} * std::sqrt(T{5} + T{2} * std::sqrt(T{10}/T{7})) };
    static const inline std::array<T, 5>
            weights = { (T{322} - T{13} * std::sqrt(T{70})) / T{900},
                        (T{322} + T{13} * std::sqrt(T{70})) / T{900},
                        T{128} / T{225},
                        (T{322} + T{13} * std::sqrt(T{70})) / T{900},
                        (T{322} - T{13} * std::sqrt(T{70})) / T{900} };
};

}

#endif