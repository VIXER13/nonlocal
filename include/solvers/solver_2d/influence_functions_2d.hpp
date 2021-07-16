#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_HPP

// В данном модуле описаны различные возмоные функции виляния

#include <cmath>
#include <array>
#include "power.hpp"

namespace nonlocal::influence {

class _superellipse_region {
    explicit _superellipse_region() noexcept = default;

    template<intmax_t n, class T>
    static T norm_pow(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        using metamath::function::power;
        if constexpr (n % 2)
            return power<n>(std::abs(x[0] - y[0]) / r[0]) + power<n>(std::abs(x[1] - y[1]) / r[1]);
        else
            return power<n>((x[0] - y[0]) / r[0]) + power<n>((x[1] - y[1]) / r[1]);
    }

public:
    template<class T, intmax_t n>
    friend class constant_2d;

    template<class T, intmax_t p, intmax_t q, intmax_t n>
    friend class polynomial_2d;

    template<class T>
    friend class normal_distribution_2d;
};

template<class T, intmax_t n = 2>
class constant_2d final {
    static_assert(n > 0, "Parameter n must be greater than 0.");

    std::array<T, 2> _r = {};
    T _norm = 0;

public:
    explicit constant_2d(const T& r) noexcept { set_radius(r); }
    explicit constant_2d(const std::array<T, 2>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 2>& r) noexcept {
        _r = r;
        using metamath::function::power;
        _norm = std::tgamma(1 + T{2} / n) / (4 * r[0] * r[1] * power<2>(std::tgamma(1 + T{1} / n)));
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        return _superellipse_region::norm_pow<n>(x, y, _r) < 1 ? _norm : T{0};
    }
};

template<class T, intmax_t p, intmax_t q, intmax_t n = 2>
class polynomial_2d final {
    static_assert(p > 0, "Parameter p must be greater than 0.");
    static_assert(q > 0, "Parameter q must be greater than 0.");
    static_assert(n > 0, "Parameter n must be greater than 0.");

    std::array<T, 2> _r = {};
    T _norm = 0;

public:
    explicit polynomial_2d(const T& r) noexcept { set_radius(r); }
    explicit polynomial_2d(const std::array<T, 2>& r) noexcept { set_radius(r); }

    void set_radius(const T& r) noexcept { set_radius(std::array{r, r}); }
    void set_radius(const std::array<T, 2>& r) noexcept {
        _r = r;
        using metamath::function::power;
        _norm = p * T{n} / (4 * r[0] * r[1] * std::beta(T{1} / n, T{1} / n) * std::beta(T{2} / p, T{q + 1}));
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        const T h = _superellipse_region::norm_pow<n>(x, y, _r);
        using metamath::function::power;
        if constexpr (p % n) // TODO: Optimize std::pow
            return h < 1 ? _norm * power<q>(T{1} - power<p / n>(h) * std::pow(h, T{p % n} / n)) : T{0};
        else
            return h < 1 ? _norm * power<q>(T{1} - power<p / n>(h)) : T{0};
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
        _disp_mul = {T{0.5} / (_r[0] * _r[0]), T{0.5} / (_r[1] * _r[1])};
        _norm = T{0.5} / (M_PI * r[0] * r[1]);
    }

    const std::array<T, 2>& radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        using metamath::function::power;
        return _norm * std::exp(-_disp_mul[0] * power<2>(x[0] - y[0]) - _disp_mul[1] * power<2>(x[1] - y[1]));
    }
};

}

#endif