#pragma once

#include <Eigen/Sparse>

namespace nonlocal::slae {

template<class T, class I>
struct preconditioner_base {
    virtual ~preconditioner_base() noexcept = default;
    virtual preconditioner_base<T, I>& analyze_pattern(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) = 0;
    virtual preconditioner_base<T, I>& factorize(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) = 0;
    virtual preconditioner_base<T, I>& compute(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) = 0;
    virtual Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const = 0;
    virtual Eigen::ComputationInfo computation_info() const = 0;
};

}