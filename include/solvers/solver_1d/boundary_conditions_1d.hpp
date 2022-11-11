#ifndef BOUNDARY_CONDITION_1D_HPP
#define BOUNDARY_CONDITION_1D_HPP

#include "../solvers_constants.hpp"

#include <memory>

namespace nonlocal {

class boundary_condition_1d {
protected:
    constexpr explicit boundary_condition_1d() noexcept = default;

public:
    virtual ~boundary_condition_1d() noexcept = default;
};

class stationary_boundary_condition_1d : public boundary_condition_1d {
protected:
    constexpr explicit stationary_boundary_condition_1d() noexcept = default;

public:
    ~stationary_boundary_condition_1d() noexcept override = default;
};

template<class T>
struct nonstationary_boundary_condition_1d : public boundary_condition_1d {
protected:
    constexpr explicit nonstationary_boundary_condition_1d() noexcept = default;

public:
    ~nonstationary_boundary_condition_1d() noexcept override = default;
    virtual std::unique_ptr<stationary_boundary_condition_1d> to_stationary(const T time) = 0;
};

template<class T>
class stationary_first_kind_1d : public virtual stationary_boundary_condition_1d {
protected:
    constexpr explicit stationary_first_kind_1d() noexcept = default;

public:
    ~stationary_first_kind_1d() noexcept override = default;
    virtual T operator()() const = 0;
};

template<class T>
class nonstationary_first_kind_1d : public virtual nonstationary_boundary_condition_1d<T> {
protected:
    constexpr explicit nonstationary_first_kind_1d() noexcept = default;

public:
    ~nonstationary_first_kind_1d() noexcept override = default;
    virtual T operator()(const T time) const = 0;
};

template<class T>
class stationary_second_kind_1d : public virtual stationary_boundary_condition_1d {
protected:
    constexpr explicit stationary_second_kind_1d() noexcept = default;

public:
    ~stationary_second_kind_1d() noexcept override = default;
    virtual T operator()() const = 0;
};

template<class T>
class nonstationary_second_kind_1d : public virtual nonstationary_boundary_condition_1d<T> {
protected:
    constexpr explicit nonstationary_second_kind_1d() noexcept = default;

public:
    ~nonstationary_second_kind_1d() noexcept override = default;
    virtual T operator()(const T time) const = 0;
};

}

#endif