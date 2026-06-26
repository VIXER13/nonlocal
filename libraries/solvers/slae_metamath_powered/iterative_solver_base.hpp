#pragma once

#include "solver_base.hpp"
// #include "preconditioner_base.hpp"
// #include "eigen_preconditioner.hpp"

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
class iterative_solver_base : public solver_base<T, I, J> {
    using _base = solver_base<T, I, J>;
    
public:
    using typename _base::floating_point_t;

private:
    // std::unique_ptr<preconditioner_base<T, I>> _preconditioner = std::make_unique<eigen_identity_preconditioner<T, I>>();
    floating_point_t _tolerance = std::numeric_limits<floating_point_t>::epsilon();
    uintmax_t _max_iterations = 10000;

protected:
    mutable uintmax_t _iterations = 0;
    mutable floating_point_t _residual = 0;

public:
    using solver_base<T, I, J>::solver_base;
    virtual ~iterative_solver_base() noexcept = default;

    // preconditioner_base<T, I>& preconditioner() noexcept {
    //     return *_preconditioner;
    // }

    // const preconditioner_base<T, I>& preconditioner() const noexcept {
    //     return *_preconditioner;
    // }

    // void init_preconditioner(std::unique_ptr<preconditioner_base<T, I>>&& preconditioner) {
    //     _preconditioner = std::move(preconditioner);
    // }

    // template<template<class, class> class Preconditioner, class... Types>
    // preconditioner_base<T, I>& init_preconditioner(Types&... args) {
    //     _preconditioner = std::make_unique<Preconditioner<T, I>>(std::forward<Types>(args)...);
    //     return preconditioner();
    // }

    floating_point_t tolerance() const noexcept {
        return _tolerance;
    }

    uintmax_t max_iterations() const noexcept {
        return _max_iterations;
    }

    void set_tolerance(const floating_point_t tolerance) noexcept {
        _tolerance = tolerance;
    }

    void set_max_iterations(const uintmax_t max_iterations) noexcept {
        _max_iterations = max_iterations;
    }

    uintmax_t iterations() const noexcept {
        return _iterations;
    }

    floating_point_t residual() const noexcept {
        return _residual;
    }
};

}