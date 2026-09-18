#pragma once 

#include "preconditioner_base.hpp"

#include <solvers/base/degree_of_freedom.hpp>
#include <logger/logger.hpp>

#include <Eigen/Eigen>

namespace nonlocal::slae {

template<class T, class Preconditioner>
class eigen_preconditioner_base : public preconditioner_base<T> {
    using _base = preconditioner_base<T>;

protected:
    using typename _base::entity_t;
    using typename _base::floating_point_t;
    using EI = int64_t;

    Preconditioner _preconditioner;
    mutable Eigen::Matrix<floating_point_t, Eigen::Dynamic, 1> _right_part;
    mutable Eigen::Matrix<floating_point_t, Eigen::Dynamic, 1> _solution;

    template<std::integral I, std::integral J>
    auto matrix_to_eigen(const metamath::linear::sparse_matrix<T, I, J>& matrix, const bool is_symmetric) const {
        Eigen::SparseMatrix<floating_point_t, Eigen::RowMajor, EI> eigen_matrix(DoF<T> * matrix.rows(), DoF<T> * matrix.cols());
        if constexpr (DoF<T> == 1) {
            for(const size_t s : std::ranges::iota_view{0zu, matrix.rows()})
                eigen_matrix.outerIndexPtr()[s + 1] = matrix.portrait.shifts[s + 1];
            eigen_matrix.data().resize(eigen_matrix.nonZeros());
            for(const size_t s : std::ranges::iota_view{0zu, matrix.portrait.shifts.back()}) {
                eigen_matrix.innerIndexPtr()[s] = matrix.portrait.indices[s];
                eigen_matrix.valuePtr()[s] = matrix.values[s];
            }
        } else {
            for(const size_t s : std::ranges::iota_view{0zu, matrix.rows()})
                for(const size_t d : std::ranges::iota_view{0zu, DoF<T>})
                    eigen_matrix.outerIndexPtr()[DoF<T> * s + d + 1] = DoF<T> * (matrix.portrait.shifts[s + 1] - matrix.portrait.shifts[s]);
            for(const size_t s : std::ranges::iota_view{0zu, DoF<T> * matrix.rows()})
                eigen_matrix.outerIndexPtr()[s + 1] += eigen_matrix.outerIndexPtr()[s];
            eigen_matrix.data().resize(eigen_matrix.nonZeros());
            size_t index = 0;
            for(const size_t row : std::ranges::iota_view{0zu, matrix.rows()})
                for(const size_t d_row : std::ranges::iota_view{0zu, DoF<T>})
                    for(const size_t col : matrix.portrait.indices_range(row)) {
                        const auto& block = matrix(row, col);
                        for(const size_t d_col : std::ranges::iota_view{0zu, DoF<T>}) {
                            eigen_matrix.innerIndexPtr()[index] = DoF<T> * col + d_col;
                            eigen_matrix.valuePtr()[index] = block[d_row][d_col];
                            ++index;
                        }
                    }
        }
        return eigen_matrix;
    }

    void right_part_to_eigen(const std::vector<entity_t>& rhs) const {
        if constexpr (DoF<T> == 1)
            for(const size_t i : std::ranges::iota_view{0u, rhs.size()})
                _right_part[i] = rhs[i];
        else
            for(const size_t i : std::ranges::iota_view{0u, rhs.size()})
                for(const size_t d : std::ranges::iota_view{0u, DoF<T>})
                    _right_part[DoF<T> * i + d] = rhs[i][d];
    }

    std::vector<entity_t> solution_to_std() const {
        std::vector<entity_t> result(_solution.size() / DoF<T>);
        if constexpr (DoF<T> == 1)
            for(const size_t i : std::ranges::iota_view{0u, result.size()})
                result[i] = _solution[i];
        else
            for(const size_t i : std::ranges::iota_view{0u, result.size()})
                for(const size_t d : std::ranges::iota_view{0u, DoF<T>})
                    result[i][d] = _solution[DoF<T> * i + d];
        return result;
    }

    template<std::integral I, std::integral J>
    explicit eigen_preconditioner_base(const metamath::linear::sparse_matrix<T, I, J>& matrix)
        : _right_part{EI(DoF<T> * matrix.cols())}
        , _solution{EI(DoF<T> * matrix.cols())} {
            if (matrix.rows() != matrix.cols())
                throw std::invalid_argument{"ILUT preconditioner requires a square matrix."};
            metamath::linear::validate_sparse_matrix(matrix);
        }

public:
    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        if (rhs.size() != _right_part.size() / DoF<T>)
            throw std::invalid_argument{"The size of the right-hand side vector does not match the number of columns in ildlt preconditioner matrix."};
        if (_preconditioner.info() != Eigen::Success)
            return rhs;
        right_part_to_eigen(rhs);
        _solution = _preconditioner.solve(_right_part);
        return solution_to_std();
    }
};

template<class T, std::integral I, std::integral J>
class ilut_eigen_preconditioner : public eigen_preconditioner_base<
        T, Eigen::IncompleteLUT<typename preconditioner_base<T>::floating_point_t, int64_t>
    > {
    using _base = eigen_preconditioner_base<
        T, Eigen::IncompleteLUT<typename preconditioner_base<T>::floating_point_t, int64_t>
    >;
    using _base::_preconditioner;
    using _base::matrix_to_eigen;

public:
    explicit ilut_eigen_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _base{matrix} {
        static constexpr bool Is_Symmetric = false;
        _preconditioner.compute(matrix_to_eigen(matrix, Is_Symmetric));
        if (_preconditioner.info() != Eigen::Success)
            logger::warning() << "ILUT preconditioner failed to compute. The solver will proceed without preconditioning.";
    }
};

template<class T, std::integral I, std::integral J>
class ildlt_eigen_preconditioner : public eigen_preconditioner_base<
        T, Eigen::IncompleteCholesky<typename preconditioner_base<T>::floating_point_t, Eigen::Upper, Eigen::NaturalOrdering<int64_t>>
    > {
    using _base = eigen_preconditioner_base<
        T, Eigen::IncompleteCholesky<typename preconditioner_base<T>::floating_point_t, Eigen::Upper, Eigen::NaturalOrdering<int64_t>>
    >;
    using _base::_preconditioner;
    using _base::matrix_to_eigen;

public:
    explicit ildlt_eigen_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _base{matrix} {
        static constexpr bool Is_Symmetric = true;
        _preconditioner.compute(matrix_to_eigen(matrix, Is_Symmetric));
        if (_preconditioner.info() != Eigen::Success)
            logger::warning() << "ILLT preconditioner failed to compute. The solver will proceed without preconditioning.";
    }
};

}