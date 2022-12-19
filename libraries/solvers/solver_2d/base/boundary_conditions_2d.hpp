#ifndef NONLOCAL_BOUNDARY_CONDITIONS_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITIONS_2D_HPP

#include "../../solvers_constants.hpp"

#include <array>
#include <memory>
#include <string>
#include <unordered_map>

namespace nonlocal {

template<class T, physics_t Physics>
class boundary_condition_2d {
protected:
    static constexpr auto from_value(const T value) noexcept {
        return [value](const std::array<T, 2>&) constexpr noexcept { return value; };
    }

    constexpr explicit boundary_condition_2d() noexcept = default;

public:
    virtual ~boundary_condition_2d() noexcept = default;
    virtual T operator()(const std::array<T, 2>& x) const = 0;
};

template<class T, physics_t Physics, size_t DoF>
class boundary_conditions_2d final : public std::array<std::unique_ptr<boundary_condition_2d<T, Physics>>, DoF> {};

template<class T, physics_t Physics>
class boundary_conditions_2d<T, Physics, 1> final : public std::unique_ptr<boundary_condition_2d<T, Physics>> {
    using _base = std::unique_ptr<boundary_condition_2d<T, Physics>>;

public:
    boundary_conditions_2d() = default;
    boundary_conditions_2d(boundary_conditions_2d<T, Physics, 1>&& condition) = default;
    boundary_conditions_2d(_base&& condition)
        : _base{std::move(condition)} {}
    boundary_conditions_2d(boundary_condition_2d<T, Physics>* const condition)
        : _base{condition} {}

    boundary_conditions_2d<T, Physics, 1>& operator=(boundary_conditions_2d<T, Physics, 1>&& condition) = default;

    boundary_conditions_2d<T, Physics, 1>& operator=(_base&& condition) {
        *this = condition.release();
        return *this;
    }

    boundary_conditions_2d<T, Physics, 1>& operator=(boundary_condition_2d<T, Physics>* const condition) {
        *this = boundary_conditions_2d{condition};
        return *this;
    }
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

#endif