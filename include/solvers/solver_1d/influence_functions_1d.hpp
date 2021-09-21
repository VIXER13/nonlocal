#ifndef NONLOCAL_INFLUENCE_FUNCTIONS_1D_HPP
#define NONLOCAL_INFLUENCE_FUNCTIONS_1D_HPP

#include "power.hpp"
#include <cmath>

namespace nonlocal::influence {

template<class T>
class constant_1d final {
    T _r, _norm;

public:
    explicit constant_1d(const T r) noexcept { set_radius(r); }

    void set_radius(const T r) noexcept {
        _r = r;
        _norm = T{0.5} / (_r);
    }

    T radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const T x, const T y) const noexcept {
        return std::abs(x - y) < _r ? _norm : T{0};
    }
};

template<class T, intmax_t p, intmax_t q>
class polynomial_1d final {
    static_assert(p > 0, "Parameter p must be greater than 0.");
    static_assert(q > 0, "Parameter q must be greater than 0.");

    T _r, _norm;

public:
    explicit polynomial_1d(const T r) noexcept { set_radius(r); }

    void set_radius(const T r) noexcept {
        _r = r;
        _norm = T{p} / (2 * _r * std::beta(T{1} / T{p}, T{q + 1}));
    }

    T radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const T x, const T y) const noexcept {
        const T h = std::abs(x - y);
        using metamath::function::power;
        return h < _r ? _norm * power<q>(T{1} - power<p>(h / _r)) : T{0};
    }
};

template<class T>
class normal_distribution_1d final {
    T _r, _norm, _disp_mul;

public:
    explicit normal_distribution_1d(const T r) noexcept { set_radius(r); }

    void set_radius(const T r) noexcept {
        _r = r;
        _norm = T{1} / (_r * std::sqrt(T{2} * M_PI));
        _disp_mul = -T{0.5} / (_r * _r);
    }

    T radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const T x, const T y) const noexcept {
        using metamath::function::power;
        return _norm * std::exp(_disp_mul * power<2>(x - y));
    }
};

}

#endif