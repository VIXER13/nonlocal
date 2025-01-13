#pragma once

#include "power.hpp"
#include <array>
#include <cmath>

namespace nonlocal::influence {

class _superellipsoid_region {
    explicit _superellipsoid_region() noexcept = default;

    template<uintmax_t N, class T>
    static T norm_pow(const std::array<T, 3>& x, const std::array<T, 3>& y, const std::array<T, 3>& r) noexcept {
        using metamath::function::power_u;
        if constexpr (N % 2)
            return power_u<N>(std::abs(x[0] - y[0]) / r[0]) +
                   power_u<N>(std::abs(x[1] - y[1]) / r[1]) +
                   power_u<N>(std::abs(x[2] - y[2]) / r[2]);
        else
            return power_u<N>((x[0] - y[0]) / r[0]) +
                   power_u<N>((x[1] - y[1]) / r[1]) +
                   power_u<N>((x[2] - y[2]) / r[2]);
    }

public:
    template<class T, uintmax_t N>
    friend class constant_3d;

    template<class T, uintmax_t P, uintmax_t Q, uintmax_t N>
    friend class polynomial_3d;

    template<class T>
    friend class normal_distribution_3d;
};

template<class T, uintmax_t N = 2>
class constant_3d final {
    static_assert(N > 0, "Parameter N must be greater than 0.");

    std::array<T, 3> _r = {};
    T _norm = T{0};

public:
    explicit constant_3d(const T& r) noexcept { set_radius(r); }
    explicit constant_3d(const std::array<T, 3>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 3>& r) noexcept {
        _r = r;
        using metamath::function::power;
        _norm = T{3} * N * N * std::tgamma(T{3} / N) / (r[0] * r[1] * r[2] * power<3>(std::tgamma(T{1} / N)));
    }

    const std::array<T, 3>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        return _superellipsoid_region::norm_pow<N>(x, y, _r) < 1 ? _norm : 0;
    }
};

template<class T, uintmax_t P, uintmax_t Q, uintmax_t N>
class polynomial_3d final {
    static_assert(P > 0, "Parameter P must be greater than 0.");
    static_assert(Q > 0, "Parameter Q must be greater than 0.");
    static_assert(N > 0, "Parameter N must be greater than 0.");

    std::array<T, 3> _r = {};
    T _norm = 0;

public:
    explicit polynomial_3d(const T& r) noexcept { set_radius(r); }
    explicit polynomial_3d(const std::array<T, 2>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 2>& r) noexcept {
        _r = r;
        using metamath::function::power;
        _norm = 8 * r[0] * r[1] * r[2] * std::beta(T{3} / P, Q + T{1}) * power<3>(std::tgamma(T{1} / N)) / (P * T{N} * N * std::tgamma(T{3} / N));
    }

    const std::array<T, 3>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        const T h = _superellipsoid_region::norm_pow<N>(x, y, _r);
        using metamath::function::power_u;
        if constexpr (P % N) // TODO: Optimize std::pow
            return h < 1 ? _norm * power_u<Q>(1 - std::pow(h, T{P} / N)) : 0;
        else
            return h < 1 ? _norm * power_u<Q>(1 - power_u<P / N>(h)) : 0;
    }
};

template<class T>
class normal_distribution_3d final {
    std::array<T, 3> _r = {}, _disp_mul = {};
    T _norm = 0;

public:
    explicit normal_distribution_3d(const T& r) noexcept { set_radius(r); }
    explicit normal_distribution_3d(const std::array<T, 3>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r, r}); }
    void set_radius(const std::array<T, 3>& r) noexcept {
        _r = r;
        _disp_mul = {-T{0.5} / (_r[0] * _r[0]),
                     -T{0.5} / (_r[1] * _r[1]),
                     -T{0.5} / (_r[2] * _r[2])};
        _norm = T{0.5} / (std::sqrt(2 * T{M_PI}) * T{M_PI} * r[0] * r[1] * r[2]);
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 3>& x, const std::array<T, 3>& y) const noexcept {
        using metamath::function::power;
        return _norm * std::exp(_disp_mul[0] * power<2>(x[0] - y[0]) +
                                _disp_mul[1] * power<2>(x[1] - y[1]) +
                                _disp_mul[2] * power<2>(x[2] - y[2]));
    }
};

}