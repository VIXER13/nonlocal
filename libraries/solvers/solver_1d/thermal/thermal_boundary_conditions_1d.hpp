#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP

#include "boundary_conditions_1d.hpp"
#include "metamath.hpp"

namespace nonlocal::thermal {

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
    explicit flux_1d(const T flux)
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
    explicit convection_1d(const T heat_transfer, const T ambient_temperature)
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
    T _boundary_temperature = T{0};

public:
    explicit radiation_1d(const T boundary_temperature)
        : _boundary_temperature{boundary_temperature} {}
    ~radiation_1d() noexcept override = default;

    T operator()() const override {
        return -STEFAN_BOLTZMANN_CONSTANT<T> * metamath::functions::power<4>(_boundary_temperature);
    }

    virtual T radiation() const {
        return T{4} * STEFAN_BOLTZMANN_CONSTANT<T> * metamath::functions::power<3>(_boundary_temperature);
    }
};

template<class T>
class combined_flux_1d : public flux_1d<T>
                       , public convection_1d<T>
                       , public radiation_1d<T> {
    using _flux = flux_1d<T>;
    using _convection = convection_1d<T>;
    using _radiation = radiation_1d<T>;

    T _absorption = T{0};
    T _emissivity = T{0};

public:
    explicit combined_flux_1d(const T absorption, const T flux,
                              const T heat_transfer, const T ambient_temperature,
                              const T emissivity, const T boundary_temperature)
        : _flux{flux}
        , _convection{heat_transfer, ambient_temperature}
        , _radiation{boundary_temperature}
        , _absorption{absorption}
        , _emissivity{emissivity} {}
    ~combined_flux_1d() noexcept override = default;

    T operator()() const override {
        return _absorption * _flux::operator()() + _convection::operator()() + _emissivity * _radiation::operator()();
    }

    T radiation() const override {
        return _emissivity * _radiation::radiation();
    }
};

template<class T>
using thermal_boundary_condition_1d = boundary_condition_1d<T, physics_t::THERMAL>;

template<class T>
using thermal_boundaries_conditions_1d = boundaries_conditions_1d<T, physics_t::THERMAL>;

}

#endif