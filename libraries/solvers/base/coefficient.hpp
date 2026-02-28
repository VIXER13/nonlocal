#pragma once

#include <metamath/types/visitor.hpp>

#include <algorithm>
#include <array>
#include <concepts>
#include <functional>
#include <ranges>
#include <variant>

namespace nonlocal {

template<std::floating_point T, size_t Dimension>
using point = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;

template<std::floating_point T, size_t Dimension>
using spatial_dependency = std::function<T(const point<T, Dimension>&)>;

template<std::floating_point T, size_t Dimension>
using solution_dependency = std::function<T(const point<T, Dimension>&, const T)>;

template<std::floating_point T, size_t Dimension>
using coefficient_t = std::variant<T, spatial_dependency<T, Dimension>, solution_dependency<T, Dimension>>;

template<std::floating_point T, size_t Dimension>
bool is_constant(const coefficient_t<T, Dimension>& coefficient) noexcept {
    return std::holds_alternative<T>(coefficient);
}

template<std::floating_point T, size_t Dimension, size_t N>
bool is_constant(const std::array<coefficient_t<T, Dimension>, N>& coefficient) {
    static constexpr auto checker = [](const coefficient_t<T, Dimension>& coefficient) { return is_constant<T, Dimension>(coefficient); };
    return std::all_of(coefficient.begin(), coefficient.end(), checker);
}

template<std::floating_point T, size_t Dimension>
T evaluate(const coefficient_t<T, Dimension>& coefficient, const point<T, Dimension>& point, const T solution) {
    return std::visit(metamath::types::visitor{
        [](const T value) noexcept { return value; },
        [&point](const spatial_dependency<T, 2u>& value) { return value(point); },
        [&point, solution](const solution_dependency<T, 2u>& value) { return value(point, solution); }
    }, coefficient);
}

template<std::floating_point T, size_t Dimension, size_t N>
std::array<T, N> evaluate(const std::array<coefficient_t<T, Dimension>, N>& coefficient,
                          const point<T, Dimension>& point, const T solution) {
    std::array<T, N> result;
    for (const size_t i : std::ranges::iota_view{0zu, N})
        result[i] = evaluate<T, Dimension>(coefficient[i], point, solution);
    return result;
}

namespace utils {

template<std::floating_point T, size_t Dimension, class Operator>
coefficient_t<T, Dimension> operation(const coefficient_t<T, Dimension>& lhs, const coefficient_t<T, Dimension>& rhs, const Operator& oper) {
    return std::visit(metamath::types::visitor{
        [&rhs, &oper](const T& lhs) {
            return std::visit(metamath::types::visitor{
                [&lhs, &oper](const T& rhs) -> coefficient_t<T, Dimension> {
                    return oper(lhs, rhs);
                },
                [&lhs, &oper](const spatial_dependency<T, Dimension>& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x) { return oper(lhs, rhs(x)); };
                },
                [&lhs, &oper](const solution_dependency<T, Dimension>& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x, const T solution) { return oper(lhs, rhs(x, solution)); };
                },
            }, rhs);
        },
        [&rhs, &oper](const spatial_dependency<T, Dimension>& lhs) {
            return std::visit(metamath::types::visitor{
                [&lhs, &oper](const T& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x) { return oper(lhs(x), rhs); };
                },
                [&lhs, &oper](const spatial_dependency<T, Dimension>& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x) { return oper(lhs(x), rhs(x)); };
                },
                [&lhs, &oper](const solution_dependency<T, Dimension>& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x, const T solution) { return oper(lhs(x), rhs(x, solution)); };
                },
            }, rhs);
        },
        [&rhs, &oper](const solution_dependency<T, Dimension>& lhs) {
            return std::visit(metamath::types::visitor{
                [&lhs, &oper](const T& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x, const T solution) { return oper(lhs(x, solution), rhs); };
                },
                [&lhs, &oper](const spatial_dependency<T, Dimension>& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x, const T solution) { return oper(lhs(x, solution), rhs(x)); };
                },
                [&lhs, &oper](const solution_dependency<T, Dimension>& rhs) -> coefficient_t<T, Dimension> {
                    return [lhs, rhs, oper](const std::array<T, Dimension>& x, const T solution) { return oper(lhs(x, solution), rhs(x, solution)); };
                },
            }, rhs);
        }
    }, lhs);
}

template<std::floating_point T, size_t Dimension>
coefficient_t<T, Dimension> operator+(const coefficient_t<T, Dimension>& lhs, const coefficient_t<T, Dimension>& rhs) {
    return operation(lhs, rhs, std::plus{});
}

template<std::floating_point T, size_t Dimension>
coefficient_t<T, Dimension> operator-(const coefficient_t<T, Dimension>& lhs, const coefficient_t<T, Dimension>& rhs) {
    return operation(lhs, rhs, std::minus{});
}

template<std::floating_point T, size_t Dimension>
coefficient_t<T, Dimension> operator*(const coefficient_t<T, Dimension>& lhs, const coefficient_t<T, Dimension>& rhs) {
    return operation(lhs, rhs, std::multiplies{});
}

template<std::floating_point T, size_t Dimension>
coefficient_t<T, Dimension> operator/(const coefficient_t<T, Dimension>& lhs, const coefficient_t<T, Dimension>& rhs) {
    return operation(lhs, rhs, std::divides{});
}

}

}