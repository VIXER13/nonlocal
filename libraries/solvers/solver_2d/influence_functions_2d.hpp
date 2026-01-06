#pragma once

#include <mesh/mesh_2d/search_function.hpp>

namespace nonlocal::solver_2d::influence {

class _influence_function_2d final {
    template<std::floating_point T>
    using power_function = T(*)(const T, const T);

    template<std::floating_point T>
    static T identity(const T value, const T) noexcept {
        return value;
    }

    template<size_t N, std::floating_point T>
    static T power(const T value, const T) noexcept {
        return metamath::functions::power<N>(value);
    }

    template<std::floating_point T>
    static T sqrt(const T value, const T) noexcept {
        return std::sqrt(value);
    }

    template<std::floating_point T>
    static T cbrt(const T value, const T) noexcept {
        return std::cbrt(value);
    }

    template<std::floating_point T>
    static power_function<T> init_power_function(const mesh::size_t_or<T>& n_variant) {
        if (std::holds_alternative<size_t>(n_variant)) {
            const size_t n = std::get<size_t>(n_variant);
            if (n == 1zu)
                return identity<T>;
            if (n == 2zu)
                return power<2, T>;
        }
        return std::pow<T, T>;
    }

    template<std::floating_point T>
    static power_function<T> init_power_function(const mesh::size_t_or<T>& n_variant, const mesh::size_t_or<T>& p_variant) {
        if (std::holds_alternative<size_t>(p_variant) && std::holds_alternative<size_t>(n_variant)) {
            const size_t p = std::get<size_t>(p_variant);
            const size_t n = std::get<size_t>(n_variant);
            if (p == n)
                return identity<T>;
            if (2 * p == n)
                return sqrt<T>;
            if (3 * p == n)
                return cbrt<T>;
            if (p == 2 * n)
                return power<2, T>;
        }
        return std::pow<T, T>;
    }

    explicit constexpr _influence_function_2d() noexcept = default;

public:
    template<std::floating_point T>
    friend class polynomial_2d;

    template<std::floating_point T>
    friend class exponential_2d;
};

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
    using _impl = _influence_function_2d;

    const mesh::influence<T> _influence;
    _impl::power_function<T> _power_func_inner;
    _impl::power_function<T> _power_func_outer;
    T _degree_inner = T{2};
    T _degree_outer = T{2};
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
        : _influence{influence}
        , _power_func_inner{_impl::init_power_function(mesh::get_value(influence.distance), p)}
        , _power_func_outer{_impl::init_power_function(q)}
        , _degree_inner{mesh::get_value(p) / mesh::get_value(mesh::get_value(influence.distance))}
        , _degree_outer{mesh::get_value(q)}
        , _norm{calc_norm(influence, mesh::get_value(p), mesh::get_value(q))} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        const T dist = distance(x, y, radius);
        return dist < T{1} ? _norm * _power_func_outer(T{1} - _power_func_inner(dist, _degree_inner), _degree_outer) : T{0};
    }
};

template<std::floating_point T>
class exponential_2d final {
    using _impl = _influence_function_2d;

    const mesh::influence<T> _influence;
    _impl::power_function<T> _power_func;
    const T _degree = T{2};
    const T _q = T{-0.5};
    const T _norm = T{1};

    static T calc_norm(const mesh::influence<T>& influence, const T p, const T q) {
        const auto& [distance, radius] = influence;
        const auto value = mesh::get_value(distance);
        if (std::holds_alternative<size_t>(value) && std::get<size_t>(value) == metamath::constants::Infinity<size_t>)
            return p * std::pow(q, 2 / p) / (8 * radius[0] * radius[1] * std::tgamma(2 / p));
        const T n = mesh::get_value(value);
        return (n * p * std::pow(T{4}, 1 / n) * std::pow(q, 2 / p)) / (8 * radius[0] * radius[1] * std::tgamma(2 / p) * std::beta(T{0.5}, 1 / n));
    }

public:
    explicit exponential_2d(const mesh::influence<T>& influence, const mesh::size_t_or<T>& p, const T q = T{0.5})
        : _influence{influence}
        , _power_func(_impl::init_power_function(mesh::get_value(influence.distance), p))
        , _degree{mesh::get_value(p) / mesh::get_value(mesh::get_value(influence.distance))}
        , _q{-q}
        , _norm{calc_norm(influence, mesh::get_value(p), q)} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        const T dist = distance(x, y, radius);
        return _norm * std::exp(_q * _power_func(dist, _degree));
    }
};

}