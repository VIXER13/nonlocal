#ifndef NONLOCAL_SLAE_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_SLAE_CONJUGATE_GRADIENT_HPP

#include "independent_symmetric_matrix_vector_product.hpp"
#include "unrelated_symmetric_matrix_vector_product.hpp"

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
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b_full = b;
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = x0.template value_or(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size()));
    Eigen::Matrix<T, Eigen::Dynamic, 1> r = b_full - matrix().template selfadjointView<Eigen::Upper>() * x;
    Eigen::Matrix<T, Eigen::Dynamic, 1> p = r;
    Eigen::Matrix<T, Eigen::Dynamic, 1> z = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size());
    T r_squared_norm = r.squaredNorm();
    const T b_norm = b_full.norm();
    _iterations = 0;
    _residual = std::sqrt(r_squared_norm) / b_norm;
    while(_iterations < max_iterations() && _residual > tolerance()) {
        _base::matrix_vector_product(z, p);
        const T nu = r_squared_norm / p.dot(z);
        x += nu * p;
        r -= nu * z;
        const T r_squared_norm_prev = std::exchange(r_squared_norm, r.squaredNorm()),
                mu = r_squared_norm / r_squared_norm_prev;
        p = r + mu * p;
        ++_iterations;
        _residual = std::sqrt(r_squared_norm) / b_norm;
    }
    return x;
}

}

#endif