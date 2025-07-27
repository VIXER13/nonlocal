#pragma once

#include <solvers/base/equation_parameters.hpp>

#include <array>
#include <memory>
#include <string>
#include <unordered_map>

namespace nonlocal {

template<class T, physics_t Physics>
class boundary_condition_2d {
protected:
    template<class U>
    static constexpr auto from_value(const U& value) {
        if constexpr (std::is_arithmetic_v<U>)
            return [result = T(value)](const std::array<T, 2>&) constexpr noexcept { return result; };
        else if constexpr (std::is_same_v<U, coefficient_t<T, 2u>>) {
            if (std::holds_alternative<T>(value))
                return spatial_dependency<T, 2u>{[result = std::get<T>(value)](const std::array<T, 2>&) constexpr noexcept { return result; }};
            if (std::holds_alternative<spatial_dependency<T, 2u>>(value))
                return std::get<spatial_dependency<T, 2u>>(value);
            throw std::domain_error{"A boundary condition cannot be constructed from this object."};
        } else // Otherwise, we assume that is a functor
            return value;
    }

    constexpr explicit boundary_condition_2d() noexcept = default;

public:
    virtual ~boundary_condition_2d() noexcept = default;
    virtual T operator()(const std::array<T, 2>& x) const = 0;
};

template<class T, physics_t Physics, size_t DoF>
class boundary_conditions_2d final : public std::array<std::unique_ptr<boundary_condition_2d<T, Physics>>, DoF> {
    using _base = std::array<std::unique_ptr<boundary_condition_2d<T, Physics>>, DoF>;

public:
    using _base::array;
    using _base::operator=;
};

template<class T, physics_t Physics>
class boundary_conditions_2d<T, Physics, 1> final : public std::unique_ptr<boundary_condition_2d<T, Physics>> {
    using _base = std::unique_ptr<boundary_condition_2d<T, Physics>>;

public:
    using _base::unique_ptr;
};

template<class T, physics_t Physics, size_t DoF>
using boundaries_conditions_2d = std::unordered_map<std::string, boundary_conditions_2d<T, Physics, DoF>>;

template<class T, physics_t Physics>
class first_kind_2d : public boundary_condition_2d<T, Physics> {
protected:
    constexpr explicit first_kind_2d() noexcept = default;

public:
    ~first_kind_2d() noexcept override = default;
};

template<class T, physics_t Physics>
class second_kind_2d : public boundary_condition_2d<T, Physics> {
protected:
    constexpr explicit second_kind_2d() noexcept = default;

public:
    ~second_kind_2d() noexcept override = default;
};

}