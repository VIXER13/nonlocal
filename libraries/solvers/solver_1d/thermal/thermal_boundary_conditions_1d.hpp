#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP

#include "boundary_conditions_1d.hpp"
#include "../../solvers_constants.hpp"
#include "metamath.hpp"

#include <array>
#include <memory>

namespace nonlocal::thermal {

template<class T>
class thermal_boundary_condition_1d : public virtual boundary_condition_1d<T> {
protected:
    explicit thermal_boundary_condition_1d() noexcept = default;

public:
    virtual ~thermal_boundary_condition_1d() noexcept = default;
};

template<class T>
class temperature_1d final : public first_kind_1d<T>
                           , public thermal_boundary_condition_1d<T> {
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
class flux_1d final : public second_kind_1d<T>
                    , public thermal_boundary_condition_1d<T> {
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
class convection_1d final : public second_kind_1d<T>
                          , public thermal_boundary_condition_1d<T> {
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
class radiation_1d final : public second_kind_1d<T>
                         , public thermal_boundary_condition_1d<T> {
    T _boundary_temperature = T{0};

public:
    explicit radiation_1d(const T boundary_temperature)
        : _boundary_temperature{boundary_temperature} {}
    ~radiation_1d() noexcept override = default;

    T operator()() const override {
        return STEFAN_BOLTZMANN_CONSTANT<T> * metamath::functions::power<4>(_boundary_temperature);
    }

    T radiation() const noexcept {
        return T{4} * STEFAN_BOLTZMANN_CONSTANT<T> * metamath::functions::power<3>(_boundary_temperature);
    }
};

template<class T>
class combined_flux_1d final : public second_kind_1d<T>
                             , public thermal_boundary_condition_1d<T> {
    T _absorption = T{0};
    flux_1d<T> _flux{0};
    convection_1d<T> _convection{0, 0};
    T _emissivity = T{0};
    radiation_1d<T> _radiation{0};

public:
    explicit combined_flux_1d(const T absorption, const T flux,
                              const T heat_transfer, const T ambient_temperature,
                              const T emissivity, const T boundary_temperature)
        : _absorption{absorption}, _flux{flux}
        , _convection{heat_transfer, ambient_temperature}
        , _emissivity{emissivity}, _radiation{boundary_temperature} {}
    ~combined_flux_1d() noexcept override = default;

    T operator()() const override {
        return _absorption * _flux() + _convection() - _emissivity * _radiation();
    }

    T heat_transfer() const noexcept {
        return _convection.heat_transfer();
    }

    T radiation() const noexcept {
        return _emissivity * _radiation.radiation();
    }
};

template<class T>
using boundary_condition_1d = std::unique_ptr<thermal_boundary_condition_1d<T>>;

template<class T>
using boundaries_conditions_1d = std::array<boundary_condition_1d<T>, 2>;

}

#endif