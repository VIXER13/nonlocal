#pragma once

#include <Eigen/Sparse>

#include <concepts>
#include <stdexcept>

namespace nonlocal::solver_1d {

struct linear_system_properties_1d final {
    bool is_symmetric = true;
    bool may_be_indefinite = false;
};

template<std::floating_point T, std::integral I>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve_linear_system_1d(
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& right_part, const linear_system_properties_1d properties) {
    if (properties.is_symmetric && !properties.may_be_indefinite) {
        const Eigen::SimplicialLLT<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>>
            solver{ matrix };
        if (solver.info() == Eigen::Success) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> solution = solver.solve(right_part);
            if (solver.info() == Eigen::Success && solution.allFinite())
                return solution;
        }
    }

    // Symmetric assemblers store only the upper triangle; LU needs the full matrix.
    Eigen::SparseMatrix<T, Eigen::ColMajor, I> full_matrix;
    if (properties.is_symmetric)
        full_matrix = matrix.template selfadjointView<Eigen::Upper>();
    else
        full_matrix = matrix;
    const Eigen::SparseLU<Eigen::SparseMatrix<T, Eigen::ColMajor, I>, Eigen::NaturalOrdering<I>> solver{ full_matrix };
    if (solver.info() != Eigen::Success)
        throw std::domain_error{ "The linear system cannot be factorized." };
    Eigen::Matrix<T, Eigen::Dynamic, 1> solution = solver.solve(right_part);
    if (solver.info() != Eigen::Success || !solution.allFinite())
        throw std::runtime_error{ "The linear system could not be solved." };
    return solution;
}

} // namespace nonlocal::solver_1d
