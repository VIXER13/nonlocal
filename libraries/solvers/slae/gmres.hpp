#pragma once

#include "iterative_solver_base.hpp"

#include <logger/logger.hpp>

#include <cmath>
#include <vector>

namespace nonlocal::slae {

// Preconditioned GMRES(m) — right-preconditioned restarted GMRES.
// For a system Ax = b it solves A*M^{-1}*(M*x) = b and recovers x = M^{-1}*(M*x).
// Restart length m is set via restart(); default 30.
template<class T, std::integral I = uint32_t, std::integral J = size_t>
class gmres final : public iterative_solver_base<T, I, J> {
    using _base = iterative_solver_base<T, I, J>;
    using _base::_iterations;
    using _base::_residual;

    uintmax_t _restart = 30;

public:
    using typename _base::entity_t;
    using typename _base::floating_point_t;
    using _base::matrix;
    using _base::tolerance;
    using _base::max_iterations;
    using _base::preconditioner;

    explicit gmres(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _base{matrix} {}

    uintmax_t restart() const noexcept { return _restart; }
    void restart(const uintmax_t m) noexcept { _restart = m; }

    std::vector<entity_t> solve(const std::vector<entity_t>& b,
                                const std::optional<std::vector<entity_t>>& x0 = std::nullopt) const override {
        using namespace metamath::linear;
        using metamath::operators::operator+;
        using metamath::operators::operator-;
        using metamath::operators::operator*;
        using metamath::operators::operator+=;
        using metamath::operators::operator-=;
        using metamath::operators::operator*=;
        logger::info() << "GMRES slae solver started" << std::endl;

        const J n = matrix().cols();
        std::vector<entity_t> x = x0.value_or(std::vector<entity_t>(n, entity_t{}));

        const floating_point_t rhs_norm = norm(b);
        if (rhs_norm == floating_point_t{0}) {
            return std::vector<entity_t>(n, entity_t{});
        }

        _iterations = 0;
        _residual = floating_point_t{1};

        const uintmax_t m = _restart;

        // Hessenberg matrix H  (m+1) x m, stored column-major: H[j] = column j of length j+2
        // V: Krylov basis vectors v_0 ... v_m
        // Givens rotation coefficients
        std::vector<std::vector<entity_t>>        V(m + 1, std::vector<entity_t>(n, entity_t{}));
        std::vector<std::vector<floating_point_t>> H(m, std::vector<floating_point_t>(m + 1, floating_point_t{0}));
        std::vector<floating_point_t> cs(m, floating_point_t{0});
        std::vector<floating_point_t> sn(m, floating_point_t{0});
        std::vector<floating_point_t> e1(m + 1, floating_point_t{0});
        std::vector<floating_point_t> y(m, floating_point_t{0});

        while (_iterations < max_iterations() && _residual > tolerance()) {
            // --- compute initial residual for this restart cycle ---
            std::vector<entity_t> r = matrix() * x;
            r *= floating_point_t{-1};
            r += b;
            const floating_point_t beta = norm(r);
            _residual = beta / rhs_norm;
            if (_residual <= tolerance())
                break;

            // v_0 = r / beta
            V[0] = r;
            V[0] *= (floating_point_t{1} / beta);

            // reset per-cycle data
            for (auto& col : H) std::fill(col.begin(), col.end(), floating_point_t{0});
            std::fill(cs.begin(), cs.end(), floating_point_t{0});
            std::fill(sn.begin(), sn.end(), floating_point_t{0});
            std::fill(e1.begin(), e1.end(), floating_point_t{0});
            e1[0] = beta;

            uintmax_t j = 0; // Arnoldi step index; also counts inner iterations
            for (; j < m && _iterations < max_iterations(); ++j) {
                // right preconditioner: w = A * M^{-1} * v_j
                const std::vector<entity_t> z = preconditioner().solve(V[j]);
                std::vector<entity_t> w = matrix() * z;

                // modified Gram-Schmidt orthogonalisation
                for (uintmax_t i = 0; i <= j; ++i) {
                    H[j][i] = static_cast<floating_point_t>(scalar_product(V[i], w));
                    w -= H[j][i] * V[i];
                }
                H[j][j + 1] = norm(w);

                if (H[j][j + 1] > floating_point_t{0})
                    V[j + 1] = w * (floating_point_t{1} / H[j][j + 1]);

                // apply previous Givens rotations to new Hessenberg column
                for (uintmax_t i = 0; i < j; ++i) {
                    const floating_point_t tmp = cs[i] * H[j][i] + sn[i] * H[j][i + 1];
                    H[j][i + 1]               = -sn[i] * H[j][i] + cs[i] * H[j][i + 1];
                    H[j][i]                   = tmp;
                }

                // compute and apply new Givens rotation for (j, j+1)
                const floating_point_t r_val = std::hypot(H[j][j], H[j][j + 1]);
                if (r_val > floating_point_t{0}) {
                    cs[j] = H[j][j]     / r_val;
                    sn[j] = H[j][j + 1] / r_val;
                } else {
                    cs[j] = floating_point_t{1};
                    sn[j] = floating_point_t{0};
                }
                H[j][j]     =  cs[j] * H[j][j] + sn[j] * H[j][j + 1];
                H[j][j + 1] = floating_point_t{0};

                e1[j + 1] = -sn[j] * e1[j];
                e1[j]     =  cs[j] * e1[j];

                _residual = std::abs(e1[j + 1]) / rhs_norm;
                ++_iterations;
                if (_residual <= tolerance())
                    break;
            }

            // back-substitution to find y such that H_upper * y = e1[0..j-1]
            const uintmax_t k = j; // number of Arnoldi steps completed this cycle
            for (uintmax_t i = k; i-- > 0;) {
                y[i] = e1[i];
                for (uintmax_t l = i + 1; l < k; ++l)
                    y[i] -= H[l][i] * y[l];
                y[i] /= H[i][i];
            }

            // update solution: x += M^{-1} * (V_k * y)
            std::vector<entity_t> dx(n, entity_t{});
            for (uintmax_t i = 0; i < k; ++i)
                dx += y[i] * V[i];
            const std::vector<entity_t> correction = preconditioner().solve(dx);
            x += correction;
        }

        logger::info() << "iterations = " << _iterations << '\n'
                       << "residual = "   << _residual   << std::endl;
        return x;
    }
};

}
