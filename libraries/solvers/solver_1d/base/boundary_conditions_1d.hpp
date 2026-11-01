#pragma once

#include "boundary_data_1d.hpp"

#include <constants/nonlocal_constants.hpp>

#include <array>
#include <memory>

namespace nonlocal::solver_1d {

template<class T, physics_t Physics>
class boundary_condition_1d {
  protected:
    constexpr explicit boundary_condition_1d() noexcept = default;

  public:
    virtual ~boundary_condition_1d() noexcept = default;
    virtual T operator()() const = 0;
    virtual boundary_data_1d<T> data() const = 0;
};

template<class T, physics_t Physics>
class first_kind_1d : public boundary_condition_1d<T, Physics> {
  protected:
    constexpr explicit first_kind_1d() noexcept = default;

  public:
    ~first_kind_1d() noexcept override = default;

    boundary_data_1d<T> data() const override {
        return { .prescribed_value = (*this)() };
    }
};

template<class T, physics_t Physics>
class second_kind_1d : public boundary_condition_1d<T, Physics> {
  protected:
    constexpr explicit second_kind_1d() noexcept = default;

  public:
    ~second_kind_1d() noexcept override = default;

    boundary_data_1d<T> data() const override {
        return { .load = (*this)() };
    }
};

template<class T, physics_t Physics>
class prescribed_value_1d final : public first_kind_1d<T, Physics> {
    T _value = T{ 0 };

  public:
    explicit prescribed_value_1d(const T value) noexcept : _value{ value } {}

    T operator()() const override {
        return _value;
    }
};

template<class T, physics_t Physics>
class prescribed_load_1d : public virtual second_kind_1d<T, Physics> {
    T _load = T{ 0 };

  public:
    explicit prescribed_load_1d(const T load) noexcept : _load{ load } {}

    T operator()() const override {
        return _load;
    }
};

template<class T, physics_t Physics>
class robin_1d : public virtual second_kind_1d<T, Physics> {
    T _coefficient = T{ 0 };
    T _reference_value = T{ 0 };

  public:
    explicit robin_1d(const T coefficient, const T reference_value) noexcept
        : _coefficient{ coefficient }, _reference_value{ reference_value } {}

    T operator()() const override {
        return _coefficient * _reference_value;
    }

    T coefficient() const noexcept {
        return _coefficient;
    }

    boundary_data_1d<T> data() const override {
        return { .load = _coefficient * _reference_value, .robin_coefficient = _coefficient };
    }
};

template<class T, physics_t Physics>
using boundaries_conditions_1d = std::array<std::unique_ptr<boundary_condition_1d<T, Physics>>, 2>;

template<class T, physics_t Physics>
boundaries_data_1d<T> normalize_boundaries_1d(const boundaries_conditions_1d<T, Physics>& conditions) {
    return { conditions.front()->data(), conditions.back()->data() };
}

} // namespace nonlocal::solver_1d
