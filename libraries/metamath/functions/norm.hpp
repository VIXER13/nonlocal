#ifndef METAMATH_NORM_HPP
#define METAMATH_NORM_HPP

#include "power.hpp"
#include "operators.hpp"

#include <array>
#include <cmath>
#include <numeric>

namespace metamath::functions {

template<size_t Exp, class T, size_t D>
T powered_norm(const std::array<T, D>& x) {
    return std::accumulate(x.begin(), x.end(), T{0}, [](const T sum, const T x) {
        if constexpr (Exp & 1)
            return sum + power<Exp>(std::abs(x));
        return sum + power<Exp>(x);
    });
}

template<class T, size_t D>
T powered_norm(const std::array<T, D>& x, const T exp) {
    return std::accumulate(x.begin(), x.end(), T{0}, [exp](const T sum, const T x) {
        return sum + std::pow(std::abs(x), exp);
    });
}

template<size_t Exp, class T, size_t D>
T norm(const std::array<T, D>& x) {
    if constexpr (Exp == 1)
        return powered_norm<Exp>(x);
    if constexpr (Exp == 2)
        return std::sqrt(powered_norm<Exp>(x));
    if constexpr (Exp == 3)
        return std::cbrt(powered_norm<Exp>(x));
    return std::pow(powered_norm<Exp>(x), T{1} / Exp);
}

template<class T, size_t D>
T norm(const std::array<T, D>& x, const T exp) {
    return std::pow(powered_norm(x, exp), T{1} / exp);
}

template<size_t N, class T, size_t D>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y) {
    return powered_norm<N>(x - y);
}

template<class T, size_t D>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y, const T exp) {
    return powered_norm(x - y, exp);
}

template<size_t N, class T, size_t D>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return powered_norm<N>(z);
}

template<class T, size_t D>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r, const T exp) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return powered_norm(z, exp);
}

template<size_t N = 2, class T, size_t D>
T distance(const std::array<T, D>& x, const std::array<T, D>& y) {
    return norm<N>(x - y);
}

template<class T, size_t D>
T distance(const std::array<T, D>& x, const std::array<T, D>& y, const T exp) {
    return norm(x - y, exp);
}

template<size_t N = 2, class T, size_t D>
T distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return norm<N>(z);
}

template<class T, size_t D>
T distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r, const T exp) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return norm(z, exp);
}

}

#endif