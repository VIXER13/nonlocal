#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_2D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_2D_HPP

#include "boundary_conditions_2d.hpp"

#include <functional>

namespace nonlocal::thermal {

class thermal_boundary_condition_2d : public virtual boundary_condition_2d {
protected:
    explicit thermal_boundary_condition_2d() noexcept = default;

public:
    virtual ~thermal_boundary_condition_2d() noexcept = default;
};

template<class T>
class temperature_2d : public first_kind_2d<T>
                     , public thermal_boundary_condition_2d {
protected:
    const std::function<T(const std::array<T, 2>&)> _temperature;

public:
    explicit temperature_2d(const T temperature)
        : _temperature{[temperature](const std::array<T, 2>&) constexpr noexcept { return temperature; }} {}
    explicit temperature_2d(const std::function<T(const std::array<T, 2>&)>& temperature)
        : _temperature{temperature} {}
    ~temperature_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _temperature(x);
    }
};

template<class T>
class flux_2d : public second_kind_2d<T>
              , public thermal_boundary_condition_2d {
protected:
    const std::function<T(const std::array<T, 2>&)> _flux;

public:
    explicit flux_2d(const T flux)
        : _flux{[flux](const std::array<T, 2>&) constexpr noexcept { return flux; }} {}
    explicit flux_2d(const std::function<T(const std::array<T, 2>&)>& flux) noexcept
        : _flux{flux} {}
    ~flux_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return _flux(x);
    }
};

template<class T>
class convection_2d : public second_kind_2d<T>
                    , public thermal_boundary_condition_2d {
protected:
    const std::function<T(const std::array<T, 2>&)> _ambient_temperature;
    const T _heat_transfer;

public:
    explicit convection_2d(const T heat_transfer, const T ambient_temperature)
        : _ambient_temperature{[ambient_temperature](const std::array<T, 2>&) constexpr noexcept { return ambient_temperature }}
        , _heat_transfer{heat_transfer} {}
    explicit convection_2d(const T heat_transfer, 
                           const std::function<T(const std::array<T, 2>&)>& ambient_temperature)
        : _ambient_temperature{ambient_temperature}
        , _heat_transfer{heat_transfer} {}
    ~convection_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        return heat_transfer * ambient_temperature(x);
    }
};

template<class T>
class radiation_2d : public second_kind_2d<T>
                   , public thermal_boundary_condition_2d {
protected:
    // parameters

public:
    explicit radiation_2d(/*init params*/) {}
    ~radiation_2d() noexcept override = default;

    T operator()(const std::array<T, 2>& x) const override {
        // What goes to the right part
        return 0;
    }
};

}

#endif