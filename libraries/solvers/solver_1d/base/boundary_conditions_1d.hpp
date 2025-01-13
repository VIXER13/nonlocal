#pragma once

#include "nonlocal_constants.hpp"

#include <array>
#include <memory>

namespace nonlocal {

template<class T, physics_t Physics>
class boundary_condition_1d {
protected:
    constexpr explicit boundary_condition_1d() noexcept = default;

public:
    virtual ~boundary_condition_1d() noexcept = default;
    virtual T operator()() const = 0;
};

template<class T, physics_t Physics>
class first_kind_1d : public boundary_condition_1d<T, Physics> {
protected:
    constexpr explicit first_kind_1d() noexcept = default;

public:
    ~first_kind_1d() noexcept override = default;
};

template<class T, physics_t Physics>
class second_kind_1d : public boundary_condition_1d<T, Physics> {
protected:
    constexpr explicit second_kind_1d() noexcept = default;

public:
    ~second_kind_1d() noexcept override = default;
};

template<class T, physics_t Physics>
using boundaries_conditions_1d = std::array<std::unique_ptr<boundary_condition_1d<T, Physics>>, 2>;

}