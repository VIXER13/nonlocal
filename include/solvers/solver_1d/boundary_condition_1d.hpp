#ifndef BOUNDARY_CONDITION_1D_HPP
#define BOUNDARY_CONDITION_1D_HPP

#include "../solvers_constants.hpp"
#include <algorithm>
#include <array>
#include <functional>
#include <utility>

namespace nonlocal {

template<class T>
using stationary_boundary_1d_t = std::pair<boundary_condition_t, T>;

template<class T>
using nonstatinary_boundary_1d_t = std::pair<boundary_condition_t, std::function<T(T)>>;

template<class Boundary_Condition>
constexpr boundary_condition_t boundary_type(const Boundary_Condition& condition) noexcept {
    return std::get<boundary_condition_t>(condition);
}

template<class Boundary_Condition>
constexpr std::array<boundary_condition_t, 2> boundary_type(const std::array<Boundary_Condition, 2>& conditions) noexcept {
    std::array<boundary_condition_t, 2> boundary_types{};
    std::transform(conditions.cbegin(), conditions.cend(), boundary_types.begin(),
        [](const Boundary_Condition& condition) { return boundary_type(condition); });
    return boundary_types;
}

template<class T>
constexpr T boundary_value(const stationary_boundary_1d_t<T>& condition) noexcept {
    return std::get<T>(condition);
}

template<class T>
T boundary_value(const nonstatinary_boundary_1d_t<T>& condition, const T t) {
    return std::get<std::function<T(T)>>(condition)(t);
}

template<class T>
stationary_boundary_1d_t<T> nonstationary_boundary_to_stationary(const nonstatinary_boundary_1d_t<T>& condition, const T t) {
    return { condition.first, condition.second(t) };
}

template<class T>
std::array<stationary_boundary_1d_t<T>, 2> nonstationary_boundary_to_stationary(const std::array<nonstatinary_boundary_1d_t<T>, 2>& conditions, const T t) {
    std::array<stationary_boundary_1d_t<T>, 2> stationary_conditions{};
    std::transform(conditions.cbegin(), conditions.cend(), stationary_conditions.begin(),
        [t](const nonstatinary_boundary_1d_t<T>& condition) { return  nonstationary_boundary_to_stationary(condition, t); });
    return stationary_conditions;
}

}

#endif