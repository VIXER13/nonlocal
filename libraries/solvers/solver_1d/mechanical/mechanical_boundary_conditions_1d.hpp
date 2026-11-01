#pragma once

#include <solvers/solver_1d/base/boundary_conditions_1d.hpp>

#include <concepts>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
using displacement_1d = prescribed_value_1d<T, physics_t::MECHANICAL>;

template<std::floating_point T>
using pressure_1d = prescribed_load_1d<T, physics_t::MECHANICAL>;

template<std::floating_point T>
class spring_1d : public robin_1d<T, physics_t::MECHANICAL> {
    using _base = robin_1d<T, physics_t::MECHANICAL>;

  public:
    using _base::_base;

    T stiffness() const noexcept {
        return _base::coefficient();
    }
};

template<class T>
class combined_loading_1d : public pressure_1d<T>, public spring_1d<T> {
  public:
    explicit combined_loading_1d(const T pressure, const T stiffness, const T displacement) noexcept
        : pressure_1d<T>{ pressure }, spring_1d<T>{ stiffness, displacement } {}

    T operator()() const override {
        return pressure_1d<T>::operator()() + spring_1d<T>::operator()();
    }

    boundary_data_1d<T> data() const override {
        auto result = spring_1d<T>::data();
        result.load = pressure_1d<T>::operator()() + spring_1d<T>::operator()();
        return result;
    }
};

template<std::floating_point T>
using mechanical_boundary_condition_1d = boundary_condition_1d<T, physics_t::MECHANICAL>;

template<std::floating_point T>
using mechanical_boundaries_conditions_1d = boundaries_conditions_1d<T, physics_t::MECHANICAL>;

} // namespace nonlocal::solver_1d::mechanical
