#pragma once

#include <mesh/mesh_2d/search_function.hpp>

namespace nonlocal::solver_2d::influence {

class _influence_function_2d final {
    template<std::floating_point T>
    static metamath::types::size_t_or<T> division(const metamath::types::size_t_or<T>& lhs, const metamath::types::size_t_or<T>& rhs) {
        if (std::holds_alternative<size_t>(rhs)) {
            if (std::get<size_t>(rhs) == metamath::constants::Infinity<size_t>) // used infinity norm
                return lhs;
            if (std::holds_alternative<size_t>(lhs) && std::get<size_t>(lhs) % std::get<size_t>(rhs) == 0)
                return std::get<size_t>(lhs) / std::get<size_t>(rhs);
        }
        return metamath::types::get_value(lhs) / metamath::types::get_value(rhs);
    }

    template<std::floating_point T>
    static T power(const T value, const metamath::types::size_t_or<T>& exp) noexcept {
        if (std::holds_alternative<T>(exp))
            return metamath::functions::power(value, std::get<T>(exp));
        const size_t n = std::get<size_t>(exp);
        if (n == 1zu)
            return value;
        return metamath::functions::power(value, n);
    }

    explicit constexpr _influence_function_2d() noexcept = default;

public:
    template<std::floating_point T>
    friend class polynomial;

    template<std::floating_point T>
    friend class exponential;
};

template<std::floating_point T>
class constant final {
    mesh::influence<T> _influence;
    T _area = T{1};

    static T calc_area(const mesh::influence<T>& influence) {
        using metamath::functions::power;
        const auto& [distance, radius] = influence;
        const T n = metamath::types::get_value(distance.exponent());
        return std::tgamma(1 + 2 / n) / (4 * radius[0] * radius[1] * power<2>(std::tgamma(1 + 1 / n)));
    }

public:
    explicit constant(const mesh::influence<T>& influence)
        : _influence{influence}
        , _area{calc_area(influence)} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        return distance(x, y, radius) <= T{1} ? _area : T{0};
    }
};

template<std::floating_point T>
class polynomial final {
    using _impl = _influence_function_2d;

    mesh::influence<T> _influence;
    metamath::types::size_t_or<T> _p;
    metamath::types::size_t_or<T> _q;
    T _norm = T{1};

    static T calc_norm(const mesh::influence<T>& influence, const T p, const T q) {
        const auto& [distance, radius] = influence;
        const auto value = distance.exponent();
        if (std::holds_alternative<size_t>(value) && std::get<size_t>(value) == metamath::constants::Infinity<size_t>)
            return p / (8 * radius[0] * radius[1] * std::beta(2 / p, q + 1));
        const T n = metamath::types::get_value(value);
        return n * p / (4 * radius[0] * radius[1] * std::beta(1 / n, 1 / n) * std::beta(2 / p, q + 1));
    }

public:
    explicit polynomial(const mesh::influence<T>& influence, const metamath::types::size_t_or<T>& p, const metamath::types::size_t_or<T>& q)
        : _influence{influence}
        , _p{_impl::division(p, influence.distance.exponent())}
        , _q{q}
        , _norm{calc_norm(influence, metamath::types::get_value(p), metamath::types::get_value(q))} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        const T dist = distance(x, y, radius);
        return dist < T{1} ? _norm * _impl::power(T{1} - _impl::power(dist, _p), _q) : T{0};
    }
};

template<std::floating_point T>
class exponential final {
    using _impl = _influence_function_2d;

    mesh::influence<T> _influence;
    metamath::types::size_t_or<T> _p;
    T _q;
    T _norm;

    static T calc_norm(const mesh::influence<T>& influence, const T p, const T q) {
        const auto& [distance, radius] = influence;
        const auto value = distance.exponent();
        if (std::holds_alternative<size_t>(value) && std::get<size_t>(value) == metamath::constants::Infinity<size_t>)
            return p * std::pow(q, 2 / p) / (8 * radius[0] * radius[1] * std::tgamma(2 / p));
        const T n = metamath::types::get_value(value);
        return (n * p * std::pow(T{4}, 1 / n) * std::pow(q, 2 / p)) / (8 * radius[0] * radius[1] * std::tgamma(2 / p) * std::beta(T{0.5}, 1 / n));
    }

public:
    // The default parameters define the normal distribution function
    explicit exponential(const mesh::influence<T>& influence, const metamath::types::size_t_or<T>& p = 2zu, const T q = T{0.5})
        : _influence{influence}
        , _p{_impl::division(p, influence.distance.exponent())}
        , _q{-q}
        , _norm{calc_norm(influence, metamath::types::get_value(p), q)} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const auto& [distance, radius] = _influence;
        const T dist = distance(x, y, radius);
        return _norm * std::exp(_q * _impl::power(dist, _p));
    }
};

// influence function for default calculations scenarious
template<std::floating_point T>
class fast_polynomial final {
    std::array<T, 2> _radius;
    T _norm;

public:
    explicit fast_polynomial(const std::array<T, 2>& radius)
        : _radius{radius}
        , _norm{T{2} / (std::numbers::pi_v<T> * radius[0] * radius[1])} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        const T dist = metamath::functions::powered_distance<2>(x, y, _radius);
        return dist < T{1} ? _norm * (T{1} - dist) : T{0};
    }
};

}