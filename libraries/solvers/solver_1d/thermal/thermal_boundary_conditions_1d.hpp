#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP

#include "boundary_conditions_1d.hpp"
#include "../../solvers_constants.hpp"
#include <cmath>
//#include "../../../metamath.hpp"

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
        : temperature{value} {}
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
        return heat_transfer * ambient_temperature - emissivity * std::pow(3, bound_temperature) + absorption * flux;
    }

    T matrix_value() const {
        return 4 * time_step * emissivity * std::pow(4, bound_temperature);
    }
};

}

#endif