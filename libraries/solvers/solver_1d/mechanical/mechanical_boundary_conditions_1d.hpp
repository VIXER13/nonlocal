#pragma once

#include <metamath/metamath.hpp>
#include <solvers/solver_1d/base/boundary_conditions_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
class displacement_1d final : public first_kind_1d<T, physics_t::MECHANICAL>{
    T _displacement = T{0};

public:
    explicit displacement_1d(const T displacement) noexcept
        : _displacement{displacement} {}
    ~displacement_1d() noexcept override = default;

    T operator()() const override {
        return _displacement;
    }
};

template<std::floating_point T>
class pressure_1d : public virtual second_kind_1d<T, physics_t::MECHANICAL> {
    T _pressure = T{0};

public:
    explicit pressure_1d(const T pressure) noexcept
        : _pressure{pressure} {}
    ~pressure_1d() noexcept override = default;

    T operator()() const override {
        return _pressure;
    }
};

template<std::floating_point T>
class spring_1d : public virtual second_kind_1d<T, physics_t::MECHANICAL> {
    T _stiffness = T{0};
    T _displacement = T{0};

public:
    explicit spring_1d(const T stiffness, const T displacement) noexcept
        : _stiffness{stiffness}
        , _displacement{displacement} {}
    ~spring_1d() noexcept override = default;

    T operator()() const override {
        return _stiffness * _displacement;
    }

    T stiffness() const noexcept {
        return _stiffness;
    }
};

template<class T>
class combined_loading_1d : public pressure_1d<T>
                          , public spring_1d<T> {
public:
    explicit combined_loading_1d(const T pressure,
                                 const T stiffness, const T displacement) noexcept
        : pressure_1d<T>{pressure}
        , spring_1d<T>{stiffness, displacement} {}
    ~combined_loading_1d() noexcept override = default;

    T operator()() const override {
        return pressure_1d<T>::operator()() + spring_1d<T>::operator()();
    }
};

template<std::floating_point T>
using mechanical_boundary_condition_1d = boundary_condition_1d<T, physics_t::MECHANICAL>;

template<std::floating_point T>
using mechanical_boundaries_conditions_1d = boundaries_conditions_1d<T, physics_t::MECHANICAL>;

}