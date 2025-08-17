#pragma once

#include <metamath/metamath.hpp>
#include <solvers/solver_1d/base/boundary_conditions_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T>
class temperature_1d final : public first_kind_1d<T, physics_t::THERMAL>{
    T _temperature = T{0};

public:
    explicit temperature_1d(const T temperature) noexcept
        : _temperature{temperature} {}
    ~temperature_1d() noexcept override = default;

    T operator()() const override {
        return _temperature;
    }
};

template<class T>
class flux_1d : public virtual second_kind_1d<T, physics_t::THERMAL> {
    T _flux = T{0};

public:
    explicit flux_1d(const T flux) noexcept
        : _flux{flux} {}
    ~flux_1d() noexcept override = default;

    T operator()() const override {
        return _flux;
    }
};

template<class T>
class convection_1d : public virtual second_kind_1d<T, physics_t::THERMAL> {
    T _heat_transfer = T{0};
    T _ambient_temperature = T{0};

public:
    explicit convection_1d(const T heat_transfer, const T ambient_temperature) noexcept
        : _heat_transfer{heat_transfer}
        , _ambient_temperature{ambient_temperature} {}
    ~convection_1d() noexcept override = default;

    T operator()() const override {
        return _heat_transfer * _ambient_temperature;
    }

    T heat_transfer() const noexcept {
        return _heat_transfer;
    }
};

template<class T>
class radiation_1d : public virtual second_kind_1d<T, physics_t::THERMAL> {
    T _emissivity = T{1};

public:
    explicit radiation_1d(const T emissivity) noexcept
        : _emissivity{emissivity} {}
    ~radiation_1d() noexcept override = default;

    T operator()() const override {
        return T{0};
    }

    T emissivity() const noexcept {
        return _emissivity;
    }
};

template<class T>
class combined_flux_1d : public flux_1d<T>
                       , public convection_1d<T>
                       , public radiation_1d<T> {
public:
    explicit combined_flux_1d(const T flux,
                              const T heat_transfer, const T ambient_temperature,
                              const T emissivity) noexcept
        : flux_1d<T>{flux}
        , convection_1d<T>{heat_transfer, ambient_temperature}
        , radiation_1d<T>{emissivity} {}
    ~combined_flux_1d() noexcept override = default;

    T operator()() const override {
        return flux_1d<T>::operator()() + convection_1d<T>::operator()() + radiation_1d<T>::operator()();
    }
};

template<class T>
using thermal_boundary_condition_1d = boundary_condition_1d<T, physics_t::THERMAL>;

template<class T>
using thermal_boundaries_conditions_1d = boundaries_conditions_1d<T, physics_t::THERMAL>;

}