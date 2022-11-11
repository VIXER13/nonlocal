#ifndef NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP
#define NONLOCAL_THERMAL_BOUNDARY_CONDITION_1D_HPP

#include "boundary_conditions_1d.hpp"

namespace nonlocal::thermal {

class stationary_thermal_boundary_condition_1d : public virtual stationary_boundary_condition_1d {
protected:
    explicit stationary_thermal_boundary_condition_1d() noexcept = default;

public:
    virtual ~stationary_thermal_boundary_condition_1d() noexcept = default;
};

template<class T>
class nonstationary_thermal_boundary_condition_1d : public virtual nonstationary_boundary_condition_1d<T> {
protected:
    explicit nonstationary_thermal_boundary_condition_1d() noexcept = default;

public:
    virtual ~nonstationary_thermal_boundary_condition_1d() noexcept = default;
};

template<class T>
class stationary_temperature_1d : public stationary_first_kind_1d<T>
                                , public stationary_thermal_boundary_condition_1d {
public:
    T temperature = T{0};

    explicit stationary_temperature_1d(const T value) noexcept
        : temperature(value) {}
    ~stationary_temperature_1d() noexcept override = default;

    T operator()() const override {
        return temperature;
    }
};

template<class T, class Functor>
class nonstationary_temperature_1d : public nonstationary_first_kind_1d<T>
                                   , public nonstationary_thermal_boundary_condition_1d<T> {
public:
    Functor temperature;

    explicit nonstationary_temperature_1d(Functor&& functor)
        : temperature{std::move(functor)} {}
    ~nonstationary_temperature_1d() noexcept override = default;

    T operator()(const T time) const override {
        return temperature(time);
    }

    std::unique_ptr<stationary_boundary_condition_1d> to_stationary(const T time) const override {
        return std::make_unique<stationary_temperature_1d<T>>(temperature(time));
    }
};

template<class T>
class stationary_flux_1d : public stationary_second_kind_1d<T>
                         , public stationary_thermal_boundary_condition_1d {
public:
    T flux = T{0};

    explicit stationary_flux_1d(const T value)
        : flux{value} {}
    ~stationary_flux_1d() noexcept override = default;

    T operator()() const override {
        return flux;
    }
};

template<class T, class Functor>
class nonstationary_flux_1d : public nonstationary_second_kind_1d<T>
                            , public nonstationary_thermal_boundary_condition_1d<T> {
public:
    Functor flux;

    explicit nonstationary_flux_1d(Functor&& functor)
        : flux{std::move(functor)} {}
    ~nonstationary_flux_1d() noexcept override = default;

    T operator()(const T time) const override {
        return flux(time);
    }

    std::unique_ptr<stationary_boundary_condition_1d> to_stationary(const T time) const override {
        return std::make_unique<stationary_flux_1d>(flux(time));
    }
};

template<class T>
class stationary_convection_1d : public stationary_second_kind_1d<T>
                               , public stationary_thermal_boundary_condition_1d {
public:
    T heat_transfer = T{0};
    T ambient_temperature = T{0};

    explicit stationary_convection_1d(const T transfer, const T temperature)
        : heat_transfer{transfer}, ambient_temperature{temperature} {}
    ~stationary_convection_1d() noexcept override = default;

    T operator()() const override {
        return heat_transfer * ambient_temperature;
    }
};

template<class T, class Functor>
class nonstationary_convection_1d : public nonstationary_second_kind_1d<T>
                                  , public nonstationary_thermal_boundary_condition_1d<T> {
public:
    T heat_transfer = T{0};
    Functor ambient_temperature;

    explicit nonstationary_convection_1d(const T transfer, Functor&& functor)
        : heat_transfer{transfer}, ambient_temperature{std::move(functor)} {}
    ~nonstationary_convection_1d() noexcept override = default;

    T operator()(const T time) const override {
        return heat_transfer * ambient_temperature(time);
    }

    std::unique_ptr<stationary_boundary_condition_1d> to_stationary(const T time) const override {
        return std::make_unique<stationary_convection_1d>(heat_transfer, ambient_temperature(time));
    }
};

}

#endif