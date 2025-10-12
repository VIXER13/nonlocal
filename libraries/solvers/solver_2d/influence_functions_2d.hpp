#pragma once

#include <math_expression/math_expression.hpp>
#include <metamath/metamath.hpp>

namespace nonlocal::solver_2d::influence {

template<std::floating_point T>
class influence_2d_base {
    std::array<T, 2> _radius = {};
    T _norm = T{0};

protected:
    influence_2d_base() noexcept = default;
    virtual ~influence_2d_base() noexcept = default;

    void set_parameters(const std::array<T, 2>& radius, const T norm) noexcept {
        _radius = radius;
        _norm = norm;
    }

    const std::array<T, 2>& radius() const noexcept { return _radius; }
    T norm() const noexcept { return _norm; }
};

template<std::floating_point T, size_t N = 2>
class constant_2d final : public influence_2d_base<T> {
    using _base = influence_2d_base<T>;

public:
    explicit constant_2d(const T radius) noexcept { set_radius(radius); }
    explicit constant_2d(const std::array<T, 2>& radius) noexcept { set_radius(radius); }
    ~constant_2d() noexcept override = default;

    void set_radius(const T radius) noexcept { set_radius(std::array{radius, radius}); }
    void set_radius(const std::array<T, 2>& radius) noexcept {
        using metamath::functions::power;
        if constexpr (N == std::numeric_limits<size_t>::max())
            _base::set_parameters(radius, T{1} / (radius[0] * radius[1]));
        else
            _base::set_parameters(radius, std::tgamma(1 + T{2} / N) / (4 * radius[0] * radius[1] * power<2>(std::tgamma(1 + T{1} / N))));
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        if constexpr (N == std::numeric_limits<size_t>::max())
            return metamath::functions::distance<N>(x, y, _base::radius()) < T{1} ? _base::norm() : T{0};
        return metamath::functions::powered_distance<N>(x, y, _base::radius()) < T{1} ? _base::norm() : T{0};
    }
};

template<std::floating_point T>
class constant_2d<T, 0> final : public influence_2d_base<T> {
    using _base = influence_2d_base<T>;

    T _n = T{2};

public:
    explicit constant_2d(const T radius, const T n) noexcept { set_parameters(radius, n); }
    explicit constant_2d(const std::array<T, 2>& radius, const T n) noexcept { set_parameters(radius, n); }
    ~constant_2d() noexcept override = default;

    void set_parameters(const T radius, const T n) noexcept { set_parameters(std::array{radius, radius}, n); }
    void set_parameters(const std::array<T, 2>& radius, const T n) noexcept {
        using metamath::functions::power;
        _n = n;
        _base::set_parameters(radius, std::tgamma(1 + T{2} / n) / (4 * radius[0] * radius[1] * power<2>(std::tgamma(1 + T{1} / n))));
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        return metamath::functions::powered_distance(x, y, _base::radius(), _n) < T{1} ? _base::norm() : T{0};
    }
};

template<std::floating_point T, size_t P, size_t Q, size_t N = 2>
class polynomial_2d final : public influence_2d_base<T> {
    static_assert(P > 0, "Parameter P must be greater than 0.");
    static_assert(Q > 0, "Parameter Q must be greater than 0.");
    static_assert(N > 0, "Parameter N must be greater than 0.");

    using _base = influence_2d_base<T>;

public:
    explicit polynomial_2d(const T& radius) noexcept { set_radius(radius); }
    explicit polynomial_2d(const std::array<T, 2>& radius) noexcept { set_radius(radius); }
    ~polynomial_2d() noexcept override = default;

    void set_radius(const T& radius) noexcept { set_radius(std::array{radius, radius}); }
    void set_radius(const std::array<T, 2>& radius) noexcept {
        const T coeff = N == std::numeric_limits<size_t>::max() ? 0.5 : N / std::beta(T{1} / N, T{1} / N);
        _base::set_parameters(radius, coeff * P / (4 * radius[0] * radius[1] * std::beta(T{2} / P, Q + T{1})));
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        using namespace metamath::functions;
        if constexpr (P % N) {
            const T dist = distance<N>(x, y, _base::radius());
            return dist < T{1} ? _base::norm() * power<Q>(T{1} - power<P>(dist)) : T{0};
        }
        const T dist = metamath::functions::powered_distance<N>(x, y, _base::radius());
        return dist < T{1} ? _base::norm() * power<Q>(1 - power<P / N>(dist)) : 0;
    }
};

template<std::floating_point T>
class polynomial_2d<T, 0, 0, 0> final : public influence_2d_base<T> {
    using _base = influence_2d_base<T>;

    T _p = 2;
    T _q = 1;
    T _n = 2;

public:
    explicit polynomial_2d(const T& radius, const T p, const T q, const T n = T{2}) noexcept {
        set_parameters(radius, p, q, n);
    }
    explicit polynomial_2d(const std::array<T, 2>& radius, const T p, const T q, const T n = T{2}) noexcept {
        set_parameters(radius, p, q, n);
    }
    ~polynomial_2d() noexcept override = default;

    void set_parameters(const T& radius, const T p, const T q, const T n = T{2}) noexcept {
        set_radius(std::array{radius, radius}, p, q, n);
    }
    void set_parameters(const std::array<T, 2>& radius, const T p, const T q, const T n = T{2}) noexcept {
        _base::set_parameters(radius, n * p / (4 * radius[0] * radius[1] * std::beta(1 / n, 1 / n) * std::beta(2 / p, q + 1)));
        _p = p;
        _q = q;
        _n = n;
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        using metamath::functions::power;
        const T distance = metamath::functions::powered_distance(x, y, _base::radius(), _n);
        return distance < T{1} ? _base::norm() * std::pow(T{1} - std::pow(distance, _p / _n), _q) : T{0};
    }
};

template<std::floating_point T, size_t P = 2, size_t N = 2>
class exponential_2d final : public influence_2d_base<T> {
    static_assert(P > 0, "Parameter P must be greater than 0.");
    static_assert(N > 0, "Parameter N must be greater than 0.");

    using _base = influence_2d_base<T>;

    T _q = -0.5;

public:
    explicit exponential_2d(const T& radius, const T q = 0.5) noexcept { set_parameters(radius, q); }
    explicit exponential_2d(const std::array<T, 2>& radius, const T q = 0.5) noexcept { set_parameters(radius, q); }
    ~exponential_2d() noexcept override = default;

    void set_parameters(const T& radius, const T q = 0.5) noexcept { set_parameters(std::array{radius, radius}, q); }
    void set_parameters(const std::array<T, 2>& radius, const T q = 0.5) noexcept {
        const T coeff = N == std::numeric_limits<size_t>::max() ? T{1} : N * std::pow(T{4}, T{1} / N) / std::beta(T{0.5}, T{1} / N);
        _base::set_parameters(radius, (coeff * P * std::pow(q, T{2} / P)) / (8 * radius[0] * radius[1] * std::tgamma(T{2} / P)));
        _q = -q;
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        using namespace metamath::functions;
        if constexpr (P % N) {
            const T dist = distance<N>(x, y, _base::radius());
            return _base::norm() * std::exp(_q * power<P>(dist));
        }
        const T dist = powered_distance<N>(x, y, _base::radius());
        return _base::norm() * std::exp(_q * power<P / N>(dist));
    }
};

template<std::floating_point T>
class exponential_2d<T, 0, 0> final : public influence_2d_base<T> {
    using _base = influence_2d_base<T>;

    T _p = 2;
    T _q = -0.5;
    T _n = 2;

public:
    explicit exponential_2d(const T& radius, const T p, const T q, const T n = 2) noexcept {
        set_parameters(radius, p, q, n);
    }
    explicit exponential_2d(const std::array<T, 2>& radius, const T p, const T q, const T n = 2) noexcept {
        set_parameters(radius, p, q, n);
    }
    ~exponential_2d() noexcept override = default;

    void set_parameters(const T& radius, const T p, const T q, const T n = T{2}) noexcept {
        set_radius(std::array{radius, radius}, p, q, n);
    }
    void set_parameters(const std::array<T, 2>& radius, const T p, const T q, const T n = T{2}) noexcept {
        const T coeff = n * std::pow(T{4}, 1 / n) / std::beta(T{0.5}, 1 / n);
        _base::set_parameters(radius, (coeff * p * std::pow(q, 2 / p)) / (8 * radius[0] * radius[1] * std::tgamma(2 / p)));
        _p = p;
        _q = -q;
        _n = n;
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const noexcept {
        using namespace metamath::functions;
        const T distance = metamath::functions::powered_distance(x, y, _base::radius(), _n);
        return _base::norm() * std::exp(_q * std::pow(distance, _p / _n));
    }
};

template<class T>
class custom_2d {
    formula::math_expression<T> _expression;

public: 
    explicit custom_2d(const std::string& expression)
        : _expression{expression} {
        if (_expression.variables_count() != 4)
            throw std::domain_error{"The variables number in the influence function shall be 4."};
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        return _expression({x[0], x[1], y[0], y[1]});
    }
};

template<class T>
class custom_limited_area_2d final : public custom_2d<T> {
    std::array<T, 2> _radius = {};

public:
    explicit custom_limited_area_2d(const std::string& expression, const std::array<T, 2> radius)
        : custom_2d<T>{expression}
        , _radius{radius} {}
    
    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        return metamath::functions::distance<2u>(x, y, _radius) < T{1} ? custom_2d<T>::operator()(x, y) : T{0};
    }
};

template<class T>
class polynomial_with_angle_2d final {
    std::array<T, 2> _radius = {};
    const T _norm = 0;

public:
    explicit polynomial_with_angle_2d(const std::array<T, 2>& radius) noexcept 
        : _radius{radius[0] * radius[0], radius[1] * radius[1]}
        , _norm{T{2} / (std::numbers::pi_v<T> * radius[0] * radius[1])} {}

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y) const {
        using metamath::functions::power;
        const T x2y2 = metamath::functions::powered_norm<2>(x);
        const T distance = ( (x2y2 - x[0] * y[0] - x[1] * y[1]) * _radius[1] + 
                            power<2>(x[1] * y[0] - x[0] * y[1]) * _radius[0]) /
                           (_radius[0] * _radius[1] * x2y2);
        return distance < T{1} ? _norm * (T{1} - distance) : T{0};
    }
};

template<class T>
using constant_2d_dynamic = constant_2d<T, 0>;
template<class T>
using polynomial_2d_dynamic = polynomial_2d<T, 0, 0, 0>;
template<class T>
using exponential_2d_dynamic = exponential_2d<T, 0, 0>;
template<class T>
using normal_distribution_2d = exponential_2d<T, 2, 2>;

}