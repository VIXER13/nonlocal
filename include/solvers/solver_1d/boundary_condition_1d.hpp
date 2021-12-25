#ifndef BOUNDARY_CONDITION_1D_HPP
#define BOUNDARY_CONDITION_1D_HPP

#include "../solvers_constants.hpp"
#include <algorithm>
#include <array>
#include <functional>
#include <utility>

namespace nonlocal {

template<class B, class T>
using stationary_boundary_1d_t = std::pair<B, T>;

template<class B, class T>
using nonstatinary_boundary_1d_t = std::pair<B, std::function<T(T)>>;

template<class B, class T, template<class, class> class Boundary_Condition>
constexpr B boundary_type(const Boundary_Condition<B, T>& condition) noexcept {
    return std::get<B>(condition);
}

template<class B, class T, template<class, class> class Boundary_Condition>
constexpr std::array<B, 2> boundary_type(const std::array<Boundary_Condition<B, T>, 2>& conditions) noexcept {
    std::array<B, 2> boundary_types{};
    std::transform(conditions.cbegin(), conditions.cend(), boundary_types.begin(),
        [](const Boundary_Condition<B, T>& condition) { return boundary_type(condition); });
    return boundary_types;
}

template<class B>
constexpr std::array<boundary_condition_t, 2> make_general_condition(const std::array<B, 2> conditions) noexcept {
    return {boundary_condition_t(conditions.front()), boundary_condition_t(conditions.back())};
}

template<class B, class T>
constexpr T boundary_value(const stationary_boundary_1d_t<B, T>& condition) noexcept {
    return std::get<T>(condition);
}

template<class B, class T>
T boundary_value(const nonstatinary_boundary_1d_t<B, T>& condition, const T t) {
    return std::get<std::function<T(T)>>(condition)(t);
}

template<class B, class T>
stationary_boundary_1d_t<B, T> nonstationary_boundary_to_stationary(const nonstatinary_boundary_1d_t<B, T>& condition, const T t) {
    return { condition.first, condition.second(t) };
}

template<class B, class T>
std::array<stationary_boundary_1d_t<B, T>, 2> nonstationary_boundary_to_stationary(const std::array<nonstatinary_boundary_1d_t<B, T>, 2>& conditions, const T t) {
    std::array<stationary_boundary_1d_t<B, T>, 2> stationary_conditions{};
    std::transform(conditions.cbegin(), conditions.cend(), stationary_conditions.begin(),
        [t](const nonstatinary_boundary_1d_t<B, T>& condition) { return  nonstationary_boundary_to_stationary(condition, t); });
    return stationary_conditions;
}

}

#endif