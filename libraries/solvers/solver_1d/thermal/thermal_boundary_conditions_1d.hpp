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
public:
    //Convection
    T heat_transfer = T{0};
    T ambient_temperature = T{0};

    //Radiation
    T emissivity = T{0};
    T bound_temperature = T{0};
    T time_step = T{0};

    //Falling flux
    T absorption = T{0};
    T flux = T{0};

    explicit convection_1d(const T transfer, const T temperature,
                           const T emis = T{0}, const T temp = T{0}, const T tau = T{0},
                           const T absorp = T{0}, const T fl = T{0})
        : heat_transfer{transfer}, ambient_temperature{temperature},
          emissivity{emis * STEFAN_BOLTZMANN_CONSTANT<T>}, bound_temperature{temp}, time_step{tau},
          absorption{absorp}, flux{fl} {}
    ~convection_1d() noexcept override = default;

    T operator()() const override {
        using namespace metamath::functions;
        return heat_transfer * ambient_temperature - emissivity * power<3>(bound_temperature) + absorption * flux;
    }

    T matrix_value() const {
        using namespace metamath::functions;
        return 4 * time_step * emissivity * power<4>(bound_temperature);
    }
};

template<class T>
using boundary_condition_1d = std::unique_ptr<thermal_boundary_condition_1d<T>>;

template<class T>
using boundaries_conditions_1d = std::array<boundary_condition_1d<T>, 2>;

}

#endif