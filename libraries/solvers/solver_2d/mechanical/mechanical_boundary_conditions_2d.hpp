#pragma once

#include "boundary_conditions_2d.hpp"

#include <functional>

namespace nonlocal::mechanical {

template<class T>
class displacement_2d final : public first_kind_2d<T, physics_t::MECHANICAL> {
    using first_kind_2d<T, physics_t::MECHANICAL>::from_value;
    const std::function<T(const std::array<T, 2>&)> _displacement;

public:
    template<class U>
    explicit displacement_2d(const U& displacement)
        : _displacement{from_value(displacement)} {}
    ~displacement_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _displacement(x);
    }
};

template<class T>
class pressure_2d final : public second_kind_2d<T, physics_t::MECHANICAL> {
    using second_kind_2d<T, physics_t::MECHANICAL>::from_value;
    const std::function<T(const std::array<T, 2>&)> _pressure;

public:
    template<class U>
    explicit pressure_2d(const U& pressure)
        : _pressure{from_value(pressure)} {}
    ~pressure_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _pressure(x);
    }
};

template<class T>
using mechanical_boundary_condition_2d = boundary_condition_2d<T, physics_t::MECHANICAL>;

template<class T>
using mechanical_boundary_conditions_2d = boundary_conditions_2d<T, physics_t::MECHANICAL, 2>;

template<class T>
using mechanical_boundaries_conditions_2d = boundaries_conditions_2d<T, physics_t::MECHANICAL, 2>;

}