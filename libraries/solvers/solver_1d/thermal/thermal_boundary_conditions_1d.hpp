#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP

#include "boundary_conditions_1d.hpp"

namespace nonlocal::thermal {

class thermal_boundary_condition_1d : public virtual boundary_condition_1d {
protected:
    explicit thermal_boundary_condition_1d() noexcept = default;

public:
    virtual ~thermal_boundary_condition_1d() noexcept = default;
};

template<class T>
class temperature_1d : public first_kind_1d<T>
                     , public thermal_boundary_condition_1d {
public:
    T temperature = T{0};

    explicit temperature_1d(const T value) noexcept
        : temperature(value) {}
    ~temperature_1d() noexcept override = default;

    T operator()() const override {
        return temperature;
    }
};

template<class T>
class flux_1d : public second_kind_1d<T>
              , public thermal_boundary_condition_1d {
public:
    T flux = T{0};

    explicit flux_1d(const T value)
        : flux{value} {}
    ~flux_1d() noexcept override = default;

    T operator()() const override {
        return flux;
    }
};

template<class T>
class convection_1d : public second_kind_1d<T>
                    , public thermal_boundary_condition_1d {
public:
    T heat_transfer = T{0};
    T ambient_temperature = T{0};

    explicit convection_1d(const T transfer, const T temperature)
        : heat_transfer{transfer}
        , ambient_temperature{temperature} {}
    ~convection_1d() noexcept override = default;

    T operator()() const override {
        return heat_transfer * ambient_temperature;
    }
};

template<class T>
class radiation_1d : public second_kind_1d<T>
                   , public thermal_boundary_condition_1d {
public:
    // parameters

    explicit radiation_1d(/*init params*/) {}
    ~radiation_1d() noexcept override = default;

    T operator()() const override {
        // What goes to the right part
        return 0;
    }
};

}

#endif