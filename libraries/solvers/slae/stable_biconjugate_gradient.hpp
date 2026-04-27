#pragma once

#include "iterative_solver_base.hpp"

namespace nonlocal::slae {

template<class T, class I>
class stable_biconjugate_gradient final : public iterative_solver_base<T, I> {
    using _base = iterative_solver_base<T, I>;
    using _base::_iterations;
    using _base::_residual;

public:
    using _base::matrix;
    using _base::tolerance;
    using _base::max_iterations;
    using _base::preconditioner;
    using _base::init_preconditioner;
    using _base::processes_ranges;

    explicit stable_biconjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);

    Eigen::Matrix<T, Eigen::Dynamic, 1> solve(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0 = std::nullopt) const override;
};

template<class T, class I>
stable_biconjugate_gradient<T, I>::stable_biconjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
    : _base{matrix} {}

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> stable_biconjugate_gradient<T, I>::solve(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
    const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0) const {
    logger::info() << "Stable BiConjugate gradient slae solver started" << std::endl;

    const I n = matrix().cols();
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = x0.template value_or(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n));
    Eigen::Matrix<T, Eigen::Dynamic, 1> r  = b - matrix() * x;
    Eigen::Matrix<T, Eigen::Dynamic, 1> r0 = r;
    Eigen::Matrix<T, Eigen::Dynamic, 1> v = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> p = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> y(n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> z(n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> s(n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> t(n);
    T r0_sqnorm = r0.squaredNorm();
    T rhs_sqnorm = b.squaredNorm();
    if(rhs_sqnorm == 0) {
        x.setZero();
        return x;
    }
    T rho   = 1;
    T alpha = 1;
    T w     = 1;
    T eps2 = std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon() * r0_sqnorm;
    I restarts = 0;

    _iterations = 0;
    _residual = std::sqrt(r.squaredNorm()) / rhs_sqnorm;
    while (_iterations < max_iterations() && _residual > tolerance()) {
        const T rho_old = rho;
        rho = r0.dot(r);
        if (std::abs(rho) < eps2) {
            // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
            // Let's restart with a new r0:
            r  = b - matrix() * x;
            r0 = r;
            rho = r.squaredNorm();
            r0_sqnorm = rho;
            if(restarts++ == 0)
                _iterations = 0;
        }

        const T beta = (rho / rho_old) * (alpha / w);
        p = r + beta * (p - w * v);
        y = preconditioner().solve(p);
        v.noalias() = matrix() * y;

        alpha = rho / r0.dot(v);
        s = r - alpha * v;
        z = preconditioner().solve(s);
        t.noalias() = matrix() * z;

        const T t_squared_norm = t.squaredNorm();
        w = t_squared_norm > T{0} ? t.dot(s) / t_squared_norm : T{0};
        x += alpha * y + w * z;
        r = s - w * t;

        _residual = std::sqrt(r.squaredNorm() / rhs_sqnorm);
        ++_iterations;
    }

    logger::info() << "iterations = " << _iterations << '\n'
                   << "residual = "   << _residual << std::endl;
    return x;
}

}