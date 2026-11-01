#pragma once

#include <metamath/metamath.hpp>
#include <solvers/solver_1d/base/boundary_conditions_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T>
using temperature_1d = prescribed_value_1d<T, physics_t::THERMAL>;

template<class T>
using flux_1d = prescribed_load_1d<T, physics_t::THERMAL>;

template<class T>
class convection_1d : public robin_1d<T, physics_t::THERMAL> {
    using _base = robin_1d<T, physics_t::THERMAL>;

  public:
    using _base::_base;

    T heat_transfer() const noexcept {
        return _base::coefficient();
    }
};

template<class T>
class radiation_1d : public virtual second_kind_1d<T, physics_t::THERMAL> {
    T _emissivity = T{ 1 };

  public:
    explicit radiation_1d(const T emissivity) noexcept : _emissivity{ emissivity } {}

    T operator()() const override {
        return T{ 0 };
    }

    T emissivity() const noexcept {
        return _emissivity;
    }

    boundary_data_1d<T> data() const override {
        using namespace metamath::constants;
        using namespace metamath::functions;
        return { .nonlinear = nonlinear_boundary_operator_1d<T>{
                     .value = [emissivity = _emissivity](const T u) constexpr noexcept {
                         return emissivity * Stefan_Boltzmann_Constant<T> * power<4>(u);
                     },
                     .derivative = [emissivity = _emissivity](const T u) constexpr noexcept {
                         return T{ 4 } * Stefan_Boltzmann_Constant<T> * emissivity * power<3>(u);
                     }
                 } };
    }
};

template<class T>
class combined_flux_1d : public flux_1d<T>, public convection_1d<T>, public radiation_1d<T> {
  public:
    explicit combined_flux_1d(const T flux, const T heat_transfer, const T ambient_temperature, const T emissivity) noexcept
        : flux_1d<T>{ flux }, convection_1d<T>{ heat_transfer, ambient_temperature }, radiation_1d<T>{ emissivity } {}

    T operator()() const override {
        return flux_1d<T>::operator()() + convection_1d<T>::operator()() + radiation_1d<T>::operator()();
    }

    boundary_data_1d<T> data() const override {
        auto result = convection_1d<T>::data();
        result.load = flux_1d<T>::operator()() + convection_1d<T>::operator()();
        result.nonlinear = radiation_1d<T>::data().nonlinear;
        return result;
    }
};

template<class T>
using thermal_boundary_condition_1d = boundary_condition_1d<T, physics_t::THERMAL>;

template<class T>
using thermal_boundaries_conditions_1d = boundaries_conditions_1d<T, physics_t::THERMAL>;

} // namespace nonlocal::solver_1d::thermal
