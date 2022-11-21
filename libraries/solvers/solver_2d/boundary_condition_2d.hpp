#ifndef NONLOCAL_BOUNDARY_CONDITION_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_2D_HPP

#include "../solvers_constants.hpp"
#include <functional>
#include <ranges>
#include <string>
#include <unordered_map>

namespace nonlocal {

template<class B, class T, size_t DoF, class... Args>
struct boundary_condition_2d_t final {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    struct node final {
        B type = B(boundary_condition_t::SECOND_KIND);
        std::function<T(Args...)> func = [](Args...) constexpr noexcept { return T{0}; };
    };

    std::conditional_t<DoF == 1, node, std::array<node, DoF>> condition;

    node& operator[]([[maybe_unused]] const size_t i) noexcept {
        if constexpr (DoF == 1) return condition;
        else                    return condition[i];
    }

    const node& operator[]([[maybe_unused]] const size_t i) const noexcept {
        if constexpr (DoF == 1) return condition;
        else                    return condition[i];
    }
};

template<class B, class T, size_t DoF>
using stationary_boundary_2d_t = boundary_condition_2d_t<B, T, DoF, const std::array<T, 2>&>;

template<class B, class T, size_t DoF>
using nonstationary_boundary_2d_t = boundary_condition_2d_t<B, T, DoF, const T, const std::array<T, 2>&>;

template<class B, class T, size_t DoF, class... Args, template<class, class, size_t, class...> class Boundary_Condition>
std::array<B, DoF> boundary_type(const Boundary_Condition<B, T, DoF, Args...>& condition) noexcept {
    std::array<B, DoF> bound_type;
    for(const size_t b : std::views::iota(size_t{0}, DoF))
        bound_type[b] = condition[b].type;
    return bound_type;
}

template<class B, class T, size_t DoF, class... Args, template<class, class, size_t, class...> class Boundary_Condition>
std::unordered_map<std::string, std::array<B, DoF>>
boundary_type(const std::unordered_map<std::string, Boundary_Condition<B, T, DoF, Args...>>& conditions) {
    std::unordered_map<std::string, std::array<B, DoF>> bound_types;
    for(const auto& [bound_name, condition] : conditions)
        bound_types[bound_name] = boundary_type(condition);
    return bound_types;
}

template<class B, class T, size_t DoF>
stationary_boundary_2d_t<B, T, DoF> to_stationary(const nonstationary_boundary_2d_t<B, T, DoF>& condition, const T t) {
    stationary_boundary_2d_t<B, T, DoF> stationary;
    for(const size_t i : std::views::iota(size_t{0}, DoF))
        stationary[i] = {
            condition[i].type,
            [t, func = condition[i].func](const std::array<T, 2>& x) { return func(t, x); }
        };
    return stationary;
}

template<class B, class T, size_t DoF>
std::unordered_map<std::string, stationary_boundary_2d_t<B, T, DoF>>
to_stationary(const std::unordered_map<std::string, nonstationary_boundary_2d_t<B, T, DoF>>& conditions, const T t) {
    std::unordered_map<std::string, stationary_boundary_2d_t<B, T, DoF>> stationary;
    for(const auto& [bound_name, condition] : conditions)
        stationary[bound_name] = to_stationary(condition, t);
    return stationary;
}

}

#endif