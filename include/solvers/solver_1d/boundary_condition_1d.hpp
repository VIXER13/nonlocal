#ifndef BOUNDARY_CONDITION_1D_HPP
#define BOUNDARY_CONDITION_1D_HPP

#include "../solvers_constants.hpp"
#include <algorithm>
#include <array>
#include <functional>
#include <utility>

namespace nonlocal {

template<class B, class T>
struct stationary_boundary_1d_t final {
    B type = B(boundary_condition_t::SECOND_KIND);
    T val = T{0};
};

template<class B, class T>
struct nonstatinary_boundary_1d_t final {
    B type = B(boundary_condition_t::SECOND_KIND);
    std::function<T(T)> func = [](const T) constexpr noexcept { return 0; };
};

template<class B, class T, template<class, class> class Boundary_Condition>
constexpr std::array<B, 2> boundary_type(const std::array<Boundary_Condition<B, T>, 2>& conditions) noexcept {
    std::array<B, 2> boundary_types{};
    std::transform(conditions.cbegin(), conditions.cend(), boundary_types.begin(),
        [](const Boundary_Condition<B, T>& condition) constexpr noexcept { return condition.type; });
    return boundary_types;
}

template<class B, class T>
stationary_boundary_1d_t<B, T> to_stationary(const nonstatinary_boundary_1d_t<B, T>& condition, const T t) {
    return {condition.type, condition.func(t)};
}

template<class B, class T>
std::array<stationary_boundary_1d_t<B, T>, 2> to_stationary(const std::array<nonstatinary_boundary_1d_t<B, T>, 2>& conditions, const T t) {
    std::array<stationary_boundary_1d_t<B, T>, 2> stationary_conditions{};
    std::transform(conditions.cbegin(), conditions.cend(), stationary_conditions.begin(),
        [t](const nonstatinary_boundary_1d_t<B, T>& condition) { return to_stationary(condition, t); });
    return stationary_conditions;
}

}

#endif