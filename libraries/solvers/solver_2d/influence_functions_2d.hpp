#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_2D_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_2D_HPP

#include "metamath.hpp"

#include <array>
#include <cinttypes>
#include <cmath>

namespace nonlocal::influence {

class _superellipse_region {
    explicit _superellipse_region() noexcept = default;

    template<uintmax_t N, class T>
    static T norm_pow(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        using metamath::functions::power;
        if constexpr (N % 2)
            return power<N>(std::abs(x[0] - y[0]) / r[0]) +
                   power<N>(std::abs(x[1] - y[1]) / r[1]);
        else
            return power<N>((x[0] - y[0]) / r[0]) +
                   power<N>((x[1] - y[1]) / r[1]);
    }

public:
    template<class T, uintmax_t N>
    friend class constant_2d;

    template<class T, uintmax_t P, uintmax_t Q, uintmax_t N>
    friend class polynomial_2d;

    template<class T>
    friend class normal_distribution_2d;
};

template<class T, uintmax_t N = 2>
class constant_2d final {
    static_assert(N > 0, "Parameter N must be greater than 0.");

    std::array<T, 2> _r = {};
    T _norm = T{0};

public:
    explicit constant_2d(const T& r) noexcept { set_radius(r); }
    explicit constant_2d(const std::array<T, 2>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 2>& r) noexcept {
        _r = r;
        using metamath::functions::power;
        _norm = std::tgamma(1 + T{2} / N) / (4 * r[0] * r[1] * power<2>(std::tgamma(1 + T{1} / N)));
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        return _superellipse_region::norm_pow<N>(x, y, _r) < 1 ? _norm : 0;
    }
};

template<class T, uintmax_t P, uintmax_t Q, uintmax_t N = 2>
class polynomial_2d final {
    static_assert(P > 0, "Parameter P must be greater than 0.");
    static_assert(Q > 0, "Parameter Q must be greater than 0.");
    static_assert(N > 0, "Parameter N must be greater than 0.");

    std::array<T, 2> _r = {};
    T _norm = 0;

public:
    explicit polynomial_2d(const T& r) noexcept { set_radius(r); }
    explicit polynomial_2d(const std::array<T, 2>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 2>& r) noexcept {
        _r = r;
        _norm = P * T{N} / (4 * r[0] * r[1] * std::beta(T{1} / N, T{1} / N) * std::beta(T{2} / P, Q + T{1}));
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        const T h = _superellipse_region::norm_pow<N>(x, y, _r);
        using metamath::functions::power;
        if constexpr (P % N) // TODO: Optimize std::pow
            return h < 1 ? _norm * power<Q>(1 - std::pow(h, T{P} / N)) : 0;
        else
            return h < 1 ? _norm * power<Q>(1 - power<P / N>(h)) : 0;
    }
};

template<class T>
class normal_distribution_2d final {
    std::array<T, 2> _r = {}, _disp_mul = {};
    T _norm = 0;

public:
    explicit normal_distribution_2d(const T& r) noexcept { set_radius(r); }
    explicit normal_distribution_2d(const std::array<T, 2>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 2>& r) noexcept {
        _r = r;
        _disp_mul = {-T{0.5} / (_r[0] * _r[0]),
                     -T{0.5} / (_r[1] * _r[1])};
        _norm = T{0.5} / (T{M_PI} * r[0] * r[1]);
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        using metamath::functions::power;
        return _norm * std::exp(_disp_mul[0] * power<2>(x[0] - y[0]) +
                                _disp_mul[1] * power<2>(x[1] - y[1]));
    }
};

}

#endif