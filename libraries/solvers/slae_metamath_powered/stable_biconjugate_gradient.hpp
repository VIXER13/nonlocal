#pragma once

#include "iterative_solver_base.hpp"

#include <logger/logger.hpp>

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
class stable_biconjugate_gradient final : public iterative_solver_base<T, I, J> {
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

    explicit stable_biconjugate_gradient(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _base{matrix} {}

    std::vector<entity_t> solve(const std::vector<entity_t>& b,
                                const std::optional<std::vector<entity_t>>& x0 = std::nullopt) const override {
        using namespace metamath::linear;
        using metamath::operators::operator+;
        using metamath::operators::operator-;
        using metamath::operators::operator*;
        using metamath::operators::operator+=;
        using metamath::operators::operator-=;
        using metamath::operators::operator*=;
        logger::info() << "Stable BiConjugate gradient slae solver started" << std::endl;

        const J n = matrix().cols();
        std::vector<entity_t> x = x0.template value_or(std::vector<entity_t>(matrix().cols(), entity_t{}));
        std::vector<entity_t> r = matrix() * x;
        r *= floating_point_t{-1};
        r += b;
        std::vector<entity_t> r0 = r;
        std::vector<entity_t> v(n, entity_t{});
        std::vector<entity_t> p(n, entity_t{});
        std::vector<entity_t> y(n, entity_t{});
        std::vector<entity_t> z(n, entity_t{});
        std::vector<entity_t> s(n, entity_t{});
        std::vector<entity_t> t(n, entity_t{});
        floating_point_t r0_sqnorm = powered_norm(r0);
        floating_point_t rhs_norm = norm(b);
        if(rhs_norm == 0) {
            x = std::vector<entity_t>(n, entity_t{});
            return x;
        }
        auto rho   = floating_point_t{1};
        auto alpha = floating_point_t{1};
        auto w     = floating_point_t{1};
        const auto eps2 = metamath::functions::power<2>(std::numeric_limits<floating_point_t>::epsilon());
        uintmax_t restarts = 0;

        _iterations = 0;
        _residual = norm(r) / rhs_norm;
        while (_iterations < max_iterations() && _residual > tolerance()) {
            const floating_point_t rho_old = rho;
            rho = scalar_production(r0, r);
            if (std::abs(rho) < eps2) {
                // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
                // Let's restart with a new r0:
                r  = matrix() * x;
                r *= floating_point_t{-1};
                r += b;
                r0 = r;
                rho = powered_norm(r);
                r0_sqnorm = rho;
                if(restarts++ == 0)
                    _iterations = 0;
            }

            const floating_point_t beta = (rho / rho_old) * (alpha / w);
            p = r + beta * (p - w * v);
            y = p; // y = preconditioner().solve(p);
            v = matrix() * y;

            alpha = rho / scalar_production(r0, v);
            s = r - alpha * v;
            z = s; // z = preconditioner().solve(s);
            t = matrix() * z;

            const floating_point_t t_squared_norm = powered_norm(t);
            w = t_squared_norm > floating_point_t{0} ? scalar_production(t, s) / t_squared_norm : floating_point_t{0};
            x += alpha * y + w * z;
            r = s - w * t;

            _residual = norm(r) / rhs_norm;
            ++_iterations;
        }

        logger::info() << "iterations = " << _iterations << '\n'
                       << "residual = "   << _residual << std::endl;
        return x;
    }
};

}