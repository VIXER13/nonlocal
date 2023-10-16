#ifndef FINITE_ELEMENT_GAUSSIAN_QUADRATURE_HPP
#define FINITE_ELEMENT_GAUSSIAN_QUADRATURE_HPP

#include "geometry_1d.hpp"
#include "geometry_primitives.hpp"

#include <cmath>

namespace metamath::finite_element {

template<class T, size_t N>
class gauss;

template<class T>
class gauss<T, 1> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 1> nodes = { T{0} };
    static inline constexpr std::array<T, 1> weights = { T{2} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 2> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 2> nodes = { T{-1} / std::sqrt(T{3}), T{1} / std::sqrt(T{3}) };
    static inline constexpr std::array<T, 2> weights = { T{1}, T{1} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 3> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 3> nodes = { -std::sqrt(T{3}/T{5}), T{0}, std::sqrt(T{3}/T{5}) };
    static inline constexpr std::array<T, 3> weights = { T{5}/T{9}, T{8}/T{9}, T{5}/T{9} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 4> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 4>
        nodes = { -std::sqrt(T{3}/T{7} + T{2}/T{7} * std::sqrt(T{6}/T{5})),
                  -std::sqrt(T{3}/T{7} - T{2}/T{7} * std::sqrt(T{6}/T{5})),
                   std::sqrt(T{3}/T{7} - T{2}/T{7} * std::sqrt(T{6}/T{5})),
                   std::sqrt(T{3}/T{7} + T{2}/T{7} * std::sqrt(T{6}/T{5})) };
    static inline constexpr std::array<T, 4>
        weights = { (T{18} - std::sqrt(T{30})) / T{36},
                    (T{18} + std::sqrt(T{30})) / T{36},
                    (T{18} + std::sqrt(T{30})) / T{36},
                    (T{18} - std::sqrt(T{30})) / T{36} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 5> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 5>
        nodes = { T{-1}/T{3} * std::sqrt(T{5} + T{2} * std::sqrt(T{10}/T{7})),
                  T{-1}/T{3} * std::sqrt(T{5} - T{2} * std::sqrt(T{10}/T{7})),
                  T{ 0},
                  T{ 1}/T{3} * std::sqrt(T{5} - T{2} * std::sqrt(T{10}/T{7})),
                  T{ 1}/T{3} * std::sqrt(T{5} + T{2} * std::sqrt(T{10}/T{7})) };
    static inline constexpr std::array<T, 5>
        weights = { (T{322} - T{13} * std::sqrt(T{70})) / T{900},
                    (T{322} + T{13} * std::sqrt(T{70})) / T{900},
                    T{128} / T{225},
                    (T{322} + T{13} * std::sqrt(T{70})) / T{900},
                    (T{322} - T{13} * std::sqrt(T{70})) / T{900} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

}

#endif