#pragma once

#include "preconditioner_base.hpp"

namespace nonlocal::slae {

// Wrapper for Eigen library preconditioners
template<class T, class I, class Preconditioner = Eigen::IdentityPreconditioner>
class eigen_preconditioner final : public preconditioner_base<T, I> {

    Preconditioner _preconditioner;

public:
    eigen_preconditioner() = default;
    explicit eigen_preconditioner(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
        : _preconditioner{matrix} {}

    Preconditioner& preconditioner() noexcept {
        return _preconditioner;
    }

    const Preconditioner& preconditioner() const noexcept {
        return _preconditioner;
    }

    preconditioner_base<T, I>& analyze_pattern(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) override {
        _preconditioner.analyzePattern(matrix);
        return *this;
    }

    preconditioner_base<T, I>& factorize(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) override {
        _preconditioner.factorize(matrix);
        return *this;
    }

    preconditioner_base<T, I>& compute(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) override {
        _preconditioner.compute(matrix);
        return *this;
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const override {
        return _preconditioner.solve(x);
    }

    Eigen::ComputationInfo computation_info() const override {
        // const_cast is used because there is an inconsistency in constancy in some of the Eigen library preconditioners.
        return const_cast<Preconditioner&>(_preconditioner).info();
    }
};

template<class T, class I>
using eigen_identity_preconditioner = eigen_preconditioner<T, I>;

template<class T, class I>
using eigen_diagonal_preconditioner = eigen_preconditioner<T, I, Eigen::DiagonalPreconditioner<T>>;


template<class T, class I>
using eigen_ILLT_preconditioner = eigen_preconditioner<T, I, Eigen::IncompleteCholesky<T, Eigen::Upper, Eigen::NaturalOrdering<I>>>;

}