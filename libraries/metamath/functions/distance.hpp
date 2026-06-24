#pragma once

#include "norm.hpp"

namespace metamath::functions {

template<size_t N = 2, class T, size_t D>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y) {
    using namespace metamath::operators;
    return powered_norm<N>(x - y);
}

template<class T, size_t D, types::arithmetic Exp>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y, const Exp exp) {
    using namespace metamath::operators;
    return powered_norm(x - y, exp);
}

template<size_t N = 2, class T, size_t D>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return powered_norm<N>(z);
}

template<class T, size_t D, types::arithmetic Exp>
T powered_distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r, const Exp exp) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return powered_norm(z, exp);
}

template<size_t N = 2, class T, size_t D>
T distance(const std::array<T, D>& x, const std::array<T, D>& y) {
    using namespace metamath::operators;
    return norm<N>(x - y);
}

template<class T, size_t D, types::arithmetic Exp>
T distance(const std::array<T, D>& x, const std::array<T, D>& y, const Exp exp) {
    using namespace metamath::operators;
    return norm(x - y, exp);
}

template<size_t N = 2, class T, size_t D>
T distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return norm<N>(z);
}

template<class T, size_t D, types::arithmetic Exp>
T distance(const std::array<T, D>& x, const std::array<T, D>& y, const std::array<T, D>& r, const Exp exp) {
    std::array<T, D> z;
    for(const size_t i : std::ranges::iota_view{0u, D})
        z[i] = (x[i] - y[i]) / r[i];
    return norm(z, exp);
}

}