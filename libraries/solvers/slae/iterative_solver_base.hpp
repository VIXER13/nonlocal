#ifndef NONLOCAL_SLAE_ITERATIVE_SOLVER_BASE_HPP
#define NONLOCAL_SLAE_ITERATIVE_SOLVER_BASE_HPP

#include "solver_base.hpp"

namespace nonlocal::slae {

template<class T, class I>
class iterative_solver_base : public solver_base<T, I> {
    T _tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    uintmax_t _max_iterations = 10000;

protected:
    mutable uintmax_t _iterations = 0;
    mutable T _residual = 0;

public:
    using solver_base<T, I>::solver_base;
    virtual ~iterative_solver_base() noexcept = default;

    T tolerance() const noexcept {
        return _tolerance;
    }

    uintmax_t max_iterations() const noexcept {
        return _max_iterations;
    }

    void set_tolerance(const T tolerance) noexcept {
        _tolerance = tolerance;
    }

    void set_max_interations(const uintmax_t max_iterations) noexcept {
        _max_iterations = max_iterations;
    }

    uintmax_t iterations() const noexcept {
        return _iterations;
    }

    T residual() const noexcept {
        
    }
};

}

#endif