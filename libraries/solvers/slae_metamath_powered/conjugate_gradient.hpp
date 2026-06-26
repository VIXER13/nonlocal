#pragma once

// #include "independent_symmetric_matrix_vector_product.hpp"
// #include "unrelated_symmetric_matrix_vector_product.hpp"
#include "iterative_solver_base.hpp"

#include <logger/logger.hpp>

#include <utility>

namespace nonlocal::slae {

enum class product_strategy : bool {
    Independent,
    Nonintersecting
};

template<class T, std::integral I, std::integral J>
class conjugate_gradient final : public iterative_solver_base<T, I, J> {
    using _base = iterative_solver_base<T, I, J>;

    using _base::_iterations;
    using _base::_residual;

public:
    using typename _base::entity_t;
    using typename _base::floating_point_t;
    using _base::matrix;
    using _base::tolerance;
    using _base::max_iterations;
    //using _base::preconditioner;
    //using _base::init_preconditioner;
    using _base::processes_ranges;

    explicit conjugate_gradient(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _base{matrix} {}

    std::vector<entity_t> solve(
        const std::vector<entity_t>& b,
        const std::optional<std::vector<entity_t>>& x0 = std::nullopt) const override {
        using namespace metamath::linear;
        using metamath::operators::operator-;
        using metamath::operators::operator*;
        using metamath::operators::operator+=;
        using metamath::operators::operator-=;
        using metamath::operators::operator*=;
        logger::info() << "Conjugate gradient slae solver started" << std::endl;
        std::vector<entity_t> z(matrix().cols(), entity_t{}); // It used as the right part in preparation calculation before iteration process
        parallel::reduce_vector(z, b);
        std::vector<entity_t> r(matrix().cols(), entity_t{});
        std::vector<entity_t> x = x0.template value_or(std::vector<entity_t>(matrix().cols(), entity_t{}));
        r = matrix().template self_adjoint<matrix_part::Upper>() * x;
        r = z - r;
        std::vector<entity_t> p = r; //std::vector<T> p = preconditioner().solve(r);
        floating_point_t r_squared_norm = scalar_production(r, p);
        const floating_point_t b_norm = norm(z);
        _iterations = 0;
        _residual = std::sqrt(r_squared_norm) / b_norm;
        while(_iterations < max_iterations() && _residual > tolerance()) {
            z = matrix().template self_adjoint<matrix_part::Upper>() * p;
            const floating_point_t nu = r_squared_norm / scalar_production(p, z);
            x += nu * p;
            z *= nu;
            r -= z;
            z = r; // z = preconditioner().solve(r);
            const floating_point_t r_squared_norm_prev = std::exchange(r_squared_norm, scalar_production(r, z));
            const floating_point_t mu = r_squared_norm / r_squared_norm_prev;
            p *= mu;
            p += z;
            ++_iterations;
            _residual = std::sqrt(r_squared_norm) / b_norm;
        }
        logger::info() << "iterations = " << _iterations << '\n'
                       << "residual = "   << _residual << std::endl;
        return x;
    }
};

}