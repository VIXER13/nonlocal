#pragma once

#include <metamath/functions/sqrt.hpp>
#include <metamath/finite_elements/finite_elements_1d/geometry/geometry_primitives.hpp>
#include <metamath/finite_elements/finite_elements_1d/geometry/geometry_1d.hpp>

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
    static inline constexpr std::array<T, 2> nodes = { T{-1} / functions::sqrt(T{3}), T{1} / functions::sqrt(T{3}) };
    static inline constexpr std::array<T, 2> weights = { T{1}, T{1} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 3> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 3> nodes = { -functions::sqrt(T{3}/T{5}), T{0}, functions::sqrt(T{3}/T{5}) };
    static inline constexpr std::array<T, 3> weights = { T{5}/T{9}, T{8}/T{9}, T{5}/T{9} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 4> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 4>
        nodes = { -functions::sqrt(T{3}/T{7} + T{2}/T{7} * functions::sqrt(T{6}/T{5})),
                  -functions::sqrt(T{3}/T{7} - T{2}/T{7} * functions::sqrt(T{6}/T{5})),
                   functions::sqrt(T{3}/T{7} - T{2}/T{7} * functions::sqrt(T{6}/T{5})),
                   functions::sqrt(T{3}/T{7} + T{2}/T{7} * functions::sqrt(T{6}/T{5})) };
    static inline constexpr std::array<T, 4>
        weights = { (T{18} - functions::sqrt(T{30})) / T{36},
                    (T{18} + functions::sqrt(T{30})) / T{36},
                    (T{18} + functions::sqrt(T{30})) / T{36},
                    (T{18} - functions::sqrt(T{30})) / T{36} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

template<class T>
class gauss<T, 5> : public geometry_1d<T, standart_segment_geometry> {
protected:
    static inline constexpr std::array<T, 5>
        nodes = { T{-1}/T{3} * functions::sqrt(T{5} + T{2} * functions::sqrt(T{10}/T{7})),
                  T{-1}/T{3} * functions::sqrt(T{5} - T{2} * functions::sqrt(T{10}/T{7})),
                  T{ 0},
                  T{ 1}/T{3} * functions::sqrt(T{5} - T{2} * functions::sqrt(T{10}/T{7})),
                  T{ 1}/T{3} * functions::sqrt(T{5} + T{2} * functions::sqrt(T{10}/T{7})) };
    static inline constexpr std::array<T, 5>
        weights = { (T{322} - T{13} * functions::sqrt(T{70})) / T{900},
                    (T{322} + T{13} * functions::sqrt(T{70})) / T{900},
                    T{128} / T{225},
                    (T{322} + T{13} * functions::sqrt(T{70})) / T{900},
                    (T{322} - T{13} * functions::sqrt(T{70})) / T{900} };

    constexpr explicit gauss() noexcept = default;
    ~gauss() noexcept override = default;
};

}
