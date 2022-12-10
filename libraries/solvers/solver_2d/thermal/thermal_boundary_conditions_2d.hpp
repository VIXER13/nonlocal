#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_2D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_2D_HPP

#include "boundary_conditions_2d.hpp"

#include <functional>

namespace nonlocal::thermal {

template<class T>
class thermal_boundary_condition_2d : public virtual boundary_condition_2d<T> {
protected:
    explicit thermal_boundary_condition_2d() noexcept = default;

public:
    virtual ~thermal_boundary_condition_2d() noexcept = default;
};

template<class T>
class temperature_2d final : public first_kind_2d<T>
                           , public thermal_boundary_condition_2d<T> {
    using first_kind_2d<T>::function_from_value;
    const std::function<T(const std::array<T, 2>&)> _temperature;

public:
    explicit temperature_2d(const T temperature)
        : _temperature{function_from_value(temperature)} {}
    explicit temperature_2d(const std::function<T(const std::array<T, 2>&)>& temperature)
        : _temperature{temperature} {}
    ~temperature_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _temperature(x);
    }
};

template<class T>
class flux_2d final : public second_kind_2d<T>
                    , public thermal_boundary_condition_2d<T> {
    using second_kind_2d<T>::function_from_value;
    const std::function<T(const std::array<T, 2>&)> _flux;

public:
    explicit flux_2d(const T flux)
        : _flux{function_from_value(flux)} {}
    explicit flux_2d(const std::function<T(const std::array<T, 2>&)>& flux) noexcept
        : _flux{flux} {}
    ~flux_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _flux(x);
    }
};

template<class T>
class convection_2d final : public second_kind_2d<T>
                          , public thermal_boundary_condition_2d<T> {
    using second_kind_2d<T>::function_from_value;
    const std::function<T(const std::array<T, 2>&)> _ambient_temperature;
    const T _heat_transfer;

public:
    explicit convection_2d(const T heat_transfer, const T ambient_temperature)
        : _ambient_temperature{function_from_value(ambient_temperature)}
        , _heat_transfer{heat_transfer} {}
    explicit convection_2d(const T heat_transfer, const std::function<T(const std::array<T, 2>&)>& ambient_temperature)
        : _ambient_temperature{ambient_temperature}
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
class radiation_2d : public second_kind_2d<T>
                   , public thermal_boundary_condition_2d<T> {
    using second_kind_2d<T>::function_from_value;
    // parameters

public:
    explicit radiation_2d(/*init params*/) {}
    ~radiation_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        // What goes to the right part
        return 0;
    }

    // some others functions...
};

template<class T>
using boundary_condition_2d = std::unique_ptr<thermal_boundary_condition_2d<T>>;

template<class T>
using boundaries_conditions_2d = std::unordered_map<std::string, boundary_condition_2d<T>>;

}

#endif