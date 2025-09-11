#pragma once

#include <metamath/metamath.hpp>

#include <cmath>

namespace nonlocal::solver_1d::influence {

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

template<class T, uintmax_t P, uintmax_t Q>
class polynomial_1d final {
    static_assert(P > 0, "Parameter P must be greater than 0.");
    static_assert(Q > 0, "Parameter Q must be greater than 0.");

    T _r, _norm;

public:
    explicit polynomial_1d(const T r) noexcept { set_radius(r); }

    void set_radius(const T r) noexcept {
        _r = r;
        _norm = P / (2 * _r * std::beta(T{1} / P, Q + T{1}));
    }

    T radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const T x, const T y) const noexcept {
        const T h = std::abs(x - y);
        using metamath::functions::power;
        return h < _r ? _norm * power<Q>(1 - power<P>(h / _r)) : 0;
    }
};

template<class T>
class normal_distribution_1d final {
    T _r, _norm, _disp_mul;

public:
    explicit normal_distribution_1d(const T r) noexcept { set_radius(r); }

    void set_radius(const T r) noexcept {
        _r = r;
        _norm = 1 / (_r * std::sqrt(2 * std::numbers::pi_v<T>));
        _disp_mul = -T{0.5} / (_r * _r);
    }

    T radius() const noexcept { return _r; }
    T norm() const noexcept { return _norm; }

    T operator()(const T x, const T y) const noexcept {
        using metamath::functions::power;
        return _norm * std::exp(_disp_mul * power<2>(x - y));
    }
};

}