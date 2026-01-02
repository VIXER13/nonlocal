#pragma once

#include <metamath/types/visitor.hpp>

#include <array>
#include <functional>
#include <unordered_map>
#include <variant>
#include <vector>

namespace nonlocal::mesh {

template<class T>
using lp_norm_parameter = std::variant<T, size_t>;

template<class T>
using distance_function = std::function<T(const std::array<T, 2>&, const std::array<T, 2>&, const std::array<T, 2>&)>;

template<class T>
T extract(const lp_norm_parameter<T>& value) {
    return std::visit(metamath::types::visitor{
        [](const T value) { return value; },
        [](const size_t value) { return T(value); }
    }, value);
}

template<class T>
struct powered_distance final {
    // x - first point, y - second point, r - normalizator (nonlocal radius)
    //using distance_t = T(*)(const std::array<T, 2>&, const std::array<T, 2>&, const std::array<T, 2>&);
    typedef T(powered_distance<T>::*distance_t)(const std::array<T, 2>&, const std::array<T, 2>&, const std::array<T, 2>&);

    const lp_norm_parameter<T> n = 2zu;
    const distance_t distance;

    explicit powered_distance(const mesh::lp_norm_parameter<T> parameter)
        : n{parameter}
        , distance{init_distance(parameter)} {}

    distance_t init_distance(const mesh::lp_norm_parameter<T> parameter) {
        return std::visit(metamath::types::visitor{
            [this](const T value) -> distance_t {
                return &powered_distance<T>::floating_exp; 
            },
            [this](const size_t value) -> distance_t {
                if (value == 2zu)
                    return &powered_distance<T>::ellipse;
                return &powered_distance<T>::integer_exp;
            }
        }, parameter);
    }

    T ellipse(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        return metamath::functions::powered_distance<2zu, T, 2zu>(x, y, r);
    }

    T integer_exp(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        return metamath::functions::powered_distance<T, 2zu>(x, y, r, std::get<size_t>(n));
    }

    T floating_exp(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        return metamath::functions::powered_distance<T, 2zu>(x, y, r, std::get<T>(n));
    }

    T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        return (this->*distance)(x, y, r);
    }
};

template<std::floating_point T>
struct powered_distance_with_rotation final {
    static T operator()(const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& r) noexcept {
        using metamath::functions::power;
        using metamath::functions::powered_norm;
        return (power<2>(r[0] * (x[0] * y[1] - x[1] * y[0])) +
                power<2>(r[1] * (x[0] * (x[0] - y[0]) + x[1] * (x[1] - y[1])))) /
               (power<2>(r[0] * r[1]) * powered_norm<2>(x));
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
