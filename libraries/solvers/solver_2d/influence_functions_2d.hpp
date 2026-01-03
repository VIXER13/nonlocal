#pragma once

#include <mesh/mesh_2d/search_function.hpp>

namespace nonlocal::solver_2d::influence {

template<std::floating_point T>
class constant_2d final {
    const mesh::influence<T> _influence;
    const T _area = T{1};

    static T calc_area(const mesh::influence<T>& influence) {
        using metamath::functions::power;
        const auto& [distance, radius] = influence;
        const T n = mesh::get_value(mesh::get_value(distance));
        return std::tgamma(1 + 2 / n) / (4 * radius[0] * radius[1] * power<2>(std::tgamma(1 + 1 / n)));
    }

public:
    explicit constant_2d(const mesh::influence<T>& influence)
        : _influence{influence}
        , _area{calc_area(influence)} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        return distance(x, y, radius) <= T{1} ? _area : T{0};
    }
};

template<std::floating_point T>
class polynomial_2d final {
    const mesh::influence<T> _influence;
    const mesh::size_t_or<T> _p;
    const mesh::size_t_or<T> _q;
    const T _norm = T{1};

    static T calc_norm(const mesh::influence<T>& influence, const T p, const T q) {
        const auto& [distance, radius] = influence;
        const auto value = mesh::get_value(distance);
        if (std::holds_alternative<size_t>(value) && std::get<size_t>(value) == metamath::constants::Infinity<size_t>)
            return p / (8 * radius[0] * radius[1] * std::beta(2 / p, q + 1));
        const T n = mesh::get_value(value);
        return n * p / (4 * radius[0] * radius[1] * std::beta(1 / n, 1 / n) * std::beta(2 / p, q + 1));
    }

public:
    explicit polynomial_2d(const mesh::influence<T>& influence, const mesh::size_t_or<T>& p, const mesh::size_t_or<T>& q)
        : _influence{influence}, _p{p}, _q{q}, _norm{calc_norm(influence, mesh::get_value(p), mesh::get_value(q))} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        using metamath::functions::power;
        const auto& [distance, radius] = _influence;
        const T n = mesh::get_value(mesh::get_value(distance));
        const T dist = distance(x, y, radius);
        return dist < T{1} ? _norm * std::pow(T{1} - std::pow(dist, mesh::get_value(_p) / n), mesh::get_value(_q)) : T{0};
    }
};

template<std::floating_point T>
class exponential_2d final {
    using power_function = T(*)(const T, const T);

    const mesh::influence<T> _influence;
    power_function _power_func;
    const T _degree = T{2};
    const T _q = T{-0.5};
    const T _norm = T{1};

    static T identity(const T value, const T) noexcept {
        return value;
    }

    static T sqr(const T value, const T) noexcept {
        return metamath::functions::power<2>(value);
    }

    static T sqrt(const T value, const T) noexcept {
        return std::sqrt(value);
    }

    static T cbrt(const T value, const T) noexcept {
        return std::cbrt(value);
    }

    static T calc_norm(const mesh::influence<T>& influence, const T p, const T q) {
        const auto& [distance, radius] = influence;
        const auto value = mesh::get_value(distance);
        if (std::holds_alternative<size_t>(value) && std::get<size_t>(value) == metamath::constants::Infinity<size_t>)
            return p * std::pow(q, 2 / p) / (8 * radius[0] * radius[1] * std::tgamma(2 / p));
        const T n = mesh::get_value(value);
        return (n * p * std::pow(T{4}, 1 / n) * std::pow(q, 2 / p)) / (8 * radius[0] * radius[1] * std::tgamma(2 / p) * std::beta(T{0.5}, 1 / n));
    }

    static T init_degree(const mesh::influence<T>& influence, const T p) {
        const auto& [distance, radius] = influence;
        return p / mesh::get_value(mesh::get_value(distance));
    }

    static power_function init_power_function(const mesh::size_t_or<T>& n_variant, const mesh::size_t_or<T>& p_variant) {
        if (std::holds_alternative<size_t>(p_variant) && std::holds_alternative<size_t>(n_variant)) {
            const size_t p = std::get<size_t>(p_variant);
            const size_t n = std::get<size_t>(n_variant);
            if (p == n)
                return identity;
            if (2 * p == n)
                return sqrt;
            if (3 * p == n)
                return cbrt;
            if (p == 2 * n)
                return sqr;
        }
        return std::pow<T, T>;
    }

public:
    explicit exponential_2d(const mesh::influence<T>& influence, const mesh::size_t_or<T>& p, const T q = T{0.5})
        : _influence{influence}
        , _power_func(init_power_function(mesh::get_value(influence.distance), p))
        , _degree{init_degree(influence, mesh::get_value(p))}
        , _q{-q}
        , _norm{calc_norm(influence, mesh::get_value(p), q)} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        const T dist = distance(x, y, radius);
        return _norm * std::exp(_q * _power_func(dist, _degree));
    }
};

}