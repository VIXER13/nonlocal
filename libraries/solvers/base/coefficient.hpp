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
struct point final : public std::array<T, Dimension> {
    point() = default;
    point(const std::array<T, Dimension>& x) : std::array<T, Dimension>{x} {}
};

template<std::floating_point T>
struct point<T, 1zu> final {
    T value = T{0};
    point() = default;
    point(T value) : value{value} {}
    T& operator=(T value) { value = value; }
    operator T() const noexcept { return value; }
};

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
    static constexpr auto checker = [](const coefficient_t<T, Dimension>& coefficient) { return is_constant(coefficient); };
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
        result[i] = evaluate(coefficient[i], point, solution);
    return result;
}

namespace utils {

template<std::floating_point T, size_t Dimension, class Operator>
coefficient_t<T, Dimension> operation(const coefficient_t<T, Dimension>& lhs, const coefficient_t<T, Dimension>& rhs, const Operator& oper) {
    using R = coefficient_t<T, Dimension>;
    using U = spatial_dependency<T, Dimension>;
    using S = solution_dependency<T, Dimension>;
    using P = point<T, Dimension>;
    return std::visit(metamath::types::visitor{
        [&oper](const T  lhs, const T  rhs) -> R { return oper(lhs, rhs); },
        [&oper](const T  lhs, const U& rhs) -> R { return [oper, lhs, rhs](const P& x)            { return oper(lhs,       rhs(x)   ); }; },
        [&oper](const T  lhs, const S& rhs) -> R { return [oper, lhs, rhs](const P& x, const T s) { return oper(lhs,       rhs(x, s)); }; },
        [&oper](const U& lhs, const T  rhs) -> R { return [oper, lhs, rhs](const P& x)            { return oper(lhs(x),    rhs      ); }; },
        [&oper](const U& lhs, const U& rhs) -> R { return [oper, lhs, rhs](const P& x)            { return oper(lhs(x),    rhs(x)   ); }; },
        [&oper](const U& lhs, const S& rhs) -> R { return [oper, lhs, rhs](const P& x, const T s) { return oper(lhs(x),    rhs(x, s)); }; },
        [&oper](const S& lhs, const T  rhs) -> R { return [oper, lhs, rhs](const P& x, const T s) { return oper(lhs(x, s), rhs      ); }; },
        [&oper](const S& lhs, const U& rhs) -> R { return [oper, lhs, rhs](const P& x, const T s) { return oper(lhs(x, s), rhs(x)   ); }; },
        [&oper](const S& lhs, const S& rhs) -> R { return [oper, lhs, rhs](const P& x, const T s) { return oper(lhs(x, s), rhs(x, s)); }; },
    }, lhs, rhs);
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