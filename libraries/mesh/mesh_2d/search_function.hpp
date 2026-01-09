#pragma once

#include <metamath/types/visitor.hpp>

#include <array>
#include <functional>
#include <unordered_map>
#include <variant>
#include <vector>

namespace nonlocal::mesh {

template<class T>
using size_t_or = std::variant<T, size_t>;

template<class T>
T get_value(const size_t_or<T>& value) {
    if (std::holds_alternative<T>(value))
        return std::get<T>(value);
    return T(std::get<size_t>(value));
}

template<class T>
struct powered_distance final {
    // x - first point, y - second point, r - normalizator (nonlocal radius)
    using distance_t = T(powered_distance<T>::*)(const std::array<T, 2>&, const std::array<T, 2>&, const std::array<T, 2>&);

    size_t_or<T> n = 2zu;
    distance_t distance = &powered_distance<T>::calculate<2zu>;

    explicit powered_distance() noexcept = default;
    explicit powered_distance(const mesh::size_t_or<T> parameter)
        : n{parameter}
        , distance{init_distance(parameter)} {}

    distance_t init_distance(const mesh::size_t_or<T> parameter) {
        return std::visit(metamath::types::visitor{
            [this](const T value) -> distance_t {
                return &powered_distance<T>::calculate_with_exp<T>; 
            },
            [this](const size_t value) -> distance_t {
                if (value == 1zu)
                    return &powered_distance<T>::calculate<1zu>;
                if (value == 2zu)
                    return &powered_distance<T>::calculate<2zu>;
                if (value == metamath::constants::Infinity<size_t>)
                    return &powered_distance<T>::calculate<metamath::constants::Infinity<size_t>>;
                return &powered_distance<T>::calculate_with_exp<size_t>;
            }
        }, parameter);
    }

    template<size_t N>
    T calculate(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) {
        return metamath::functions::powered_distance<N, T, 2zu>(x, y, r);
    }

    template<class U>
    T calculate_with_exp(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) {
        return metamath::functions::powered_distance<T, 2zu>(x, y, r, std::get<U>(n));
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) {
        return (this->*distance)(x, y, r);
    }
};

template<std::floating_point T>
struct powered_distance_with_rotation final {
    static constexpr size_t_or<T> n = 2zu;

    static T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) {
        using metamath::functions::power;
        using metamath::functions::powered_norm;
        return (power<2>(r[0] * (x[0] * y[1] - x[1] * y[0])) +
                power<2>(r[1] * (x[0] * (x[0] - y[0]) + x[1] * (x[1] - y[1])))) /
               (power<2>(r[0] * r[1]) * powered_norm<2>(x));
    }
};

template<class T>
struct distance_function final {
    std::variant<powered_distance<T>, powered_distance_with_rotation<T>> function = powered_distance<T>{};

    const size_t_or<T>& exponent() const {
        if (std::holds_alternative<powered_distance<T>>(function))
            return std::get<powered_distance<T>>(function).n;
        return std::get<powered_distance_with_rotation<T>>(function).n;
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) const {
        if (std::holds_alternative<powered_distance<T>>(function))
            return const_cast<powered_distance<T>&>(std::get<powered_distance<T>>(function))(x, y, r);
        return std::get<powered_distance_with_rotation<T>>(function)(x, y, r);
    }
};

template<class T>
struct influence final {
    distance_function<T> distance;
    std::array<T, 2> radius;
};

template<class T>
using influences = std::unordered_map<std::string, influence<T>>;

template<class T, class I>
using neighbours_t = std::pair<influences<T>, std::vector<std::vector<I>>>;

}