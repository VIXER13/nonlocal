#pragma once

#include "boundary_conditions_2d.hpp"

#include <functional>

namespace nonlocal::thermal {

template<class T>
class temperature_2d final : public first_kind_2d<T, physics_t::THERMAL> {
    using first_kind_2d<T, physics_t::THERMAL>::from_value;
    const std::function<T(const std::array<T, 2>&)> _temperature;

public:
    template<class U>
    explicit temperature_2d(const U& temperature)
        : _temperature{from_value(temperature)} {}
    ~temperature_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _temperature(x);
    }
};

template<class T>
class flux_2d : public virtual second_kind_2d<T, physics_t::THERMAL> {
    using second_kind_2d<T, physics_t::THERMAL>::from_value;
    const std::function<T(const std::array<T, 2>&)> _flux;

public:
    template<class U>
    explicit flux_2d(const U& flux)
        : _flux{from_value(flux)} {}
    ~flux_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _flux(x);
    }
};

template<class T>
class convection_2d : public virtual second_kind_2d<T, physics_t::THERMAL> {
    using second_kind_2d<T, physics_t::THERMAL>::from_value;
    const std::function<T(const std::array<T, 2>&)> _ambient_temperature;
    const T _heat_transfer;

public:
    template<class U>
    explicit convection_2d(const T heat_transfer, const U& ambient_temperature)
        : _ambient_temperature{from_value(ambient_temperature)}
        , _heat_transfer{heat_transfer} {}
    ~convection_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _heat_transfer * _ambient_temperature(x);
    }

    constexpr T heat_transfer() const noexcept {
        return _heat_transfer;
    }
};

template<class T>
class radiation_2d : public virtual second_kind_2d<T, physics_t::THERMAL> {
    T _emissivity = T{1};

public:
    explicit radiation_2d(const T emissivity)
        : _emissivity{emissivity} {}
    ~radiation_2d() noexcept override = default;

    T operator()(const std::array<T, 2>&) const override {
        return T{0};
    }

    T emissivity() const noexcept {
        return _emissivity;
    }
};

template<class T>
class combined_flux_2d : public flux_2d<T>
                       , public convection_2d<T>
                       , public radiation_2d<T> {
public:
    template<class Flux, class Ambient_Temperature>
    explicit combined_flux_2d(const Flux& flux,
                              const T heat_transfer, const Ambient_Temperature& ambient_temperature,
                              const T emissivity)
        : flux_2d<T>{flux}
        , convection_2d<T>{heat_transfer, ambient_temperature}
        , radiation_2d<T>{emissivity} {}
    ~combined_flux_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return flux_2d<T>::operator()(x) + convection_2d<T>::operator()(x); // the radiation term is implemented in a separate function
    }
};

template<class T>
using thermal_boundary_condition_2d = boundary_conditions_2d<T, nonlocal::physics_t::THERMAL, 1>;

template<class T>
using thermal_boundaries_conditions_2d = boundaries_conditions_2d<T, physics_t::THERMAL, 1>;

}