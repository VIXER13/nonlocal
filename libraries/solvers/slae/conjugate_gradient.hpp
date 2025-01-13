#pragma once

#include "independent_symmetric_matrix_vector_product.hpp"
#include "unrelated_symmetric_matrix_vector_product.hpp"

#include <vector>

namespace nonlocal::slae {

enum class product_strategy : bool {
    INDEPENDENT,
    UNRELATED
};

template<class T, class I, product_strategy Strategy = product_strategy::INDEPENDENT>
class conjugate_gradient final : public std::conditional_t<
                                            Strategy == product_strategy::INDEPENDENT,
                                            independent_symmetric_matrix_vector_product<T, I>,
                                            unrelated_symmetric_matrix_vector_product<T, I>
                                        > {
    using _base = std::conditional_t<
        Strategy == product_strategy::INDEPENDENT,
        independent_symmetric_matrix_vector_product<T, I>,
        unrelated_symmetric_matrix_vector_product<T, I>
    >;
    using _base::_iterations;
    using _base::_residual;

public:
    using _base::matrix;
    using _base::tolerance;
    using _base::max_iterations;
    using _base::preconditioner;
    using _base::init_preconditioner;
    using _base::processes_ranges;

    explicit conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);

    Eigen::Matrix<T, Eigen::Dynamic, 1> solve(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0 = std::nullopt) const override;
};

template<class T, class I, product_strategy Strategy>
conjugate_gradient<T, I, Strategy>::conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
    : _base{matrix} {}

template<class T, class I, product_strategy Strategy>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient<T, I, Strategy>::solve(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
    const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0) const {
    logger::get().log() << "Conjugate gradient slae solver started" << std::endl;
    Eigen::Matrix<T, Eigen::Dynamic, 1> z = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(matrix().cols()); // It used as the right part in preparation calculation before iteration process
    parallel::reduce_vector(z, b);
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = x0.template value_or(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(z.size()));
    Eigen::Matrix<T, Eigen::Dynamic, 1> r = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(x.size());
    _base::matrix_vector_product(r, x);
    r = z - r;
    Eigen::Matrix<T, Eigen::Dynamic, 1> p = preconditioner().solve(r);
    T r_squared_norm = r.dot(p);
    const T b_norm = z.norm();
    _iterations = 0;
    _residual = std::sqrt(r_squared_norm) / b_norm;
    while(_iterations < max_iterations() && _residual > tolerance()) {
        _base::matrix_vector_product(z, p);
        const T nu = r_squared_norm / p.dot(z);
        x += nu * p;
        r -= nu * z;
        z = preconditioner().solve(r);
        const T r_squared_norm_prev = std::exchange(r_squared_norm, r.dot(z));
        const T mu = r_squared_norm / r_squared_norm_prev;
        p = z + mu * p;
        ++_iterations;
        _residual = std::sqrt(r_squared_norm) / b_norm;
    }
    logger::get().log() << "iterations = " << _iterations << '\n'
                        << "residual = "   << _residual << std::endl;
    return x;
}

}