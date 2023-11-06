#ifndef NONLOCAL_SLAE_UNRELATED_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP
#define NONLOCAL_SLAE_UNRELATED_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP

#include "iterative_solver_base.hpp"

#include "unrelated_rows.hpp"

namespace nonlocal::slae {

template<class T, class I>
class unrelated_symmetric_matrix_vector_product : public iterative_solver_base<T, I> {
    using _base = iterative_solver_base<T, I>;

#ifdef MPI_BUILD
    mutable Eigen::Matrix<T, Eigen::Dynamic, 1> _sendbuf;
#endif
    parallel_utils::unrelated_rows<I> _unrelated;

protected:
    void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& result,
                               const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const;

public:
    using _base::matrix;
    using _base::process_rows;
    using _base::threads_count;

    explicit unrelated_symmetric_matrix_vector_product(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);
};

template<class T, class I>
unrelated_symmetric_matrix_vector_product<T, I>::unrelated_symmetric_matrix_vector_product(
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
    : iterative_solver_base<T, I>{matrix}
#ifdef MPI_BUILD
    , _sendbuf{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(matrix.cols())}
#endif
    , _unrelated{matrix} {}

template<class T, class I>
void unrelated_symmetric_matrix_vector_product<T, I>::matrix_vector_product(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& result,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1>& product =
#ifdef MPI_BUILD
    _sendbuf;
#else
    result;
#endif

    product.setZero();
    const size_t processor_shift = _base::process_rows().front();
    for(const size_t shift_index : std::ranges::iota_view<size_t, size_t>{0u, _unrelated.shifts.size() - 1}) {
#pragma omp parallel for num_threads(threads_count())
        for(size_t shift = _unrelated.shifts[shift_index]; shift < _unrelated.shifts[shift_index + 1]; ++shift) {
            const I row = _unrelated.rows[shift];
            const I glob_row = row + processor_shift;
            const I ind = matrix().outerIndexPtr()[row];
            product[glob_row] += matrix().valuePtr()[ind] * vector[matrix().innerIndexPtr()[ind]];
            for(const I i : std::ranges::iota_view{ind + 1, matrix().outerIndexPtr()[row + 1]}) {
                product[glob_row] += matrix().valuePtr()[i] * vector[matrix().innerIndexPtr()[i]];
                product[matrix().innerIndexPtr()[i]] += matrix().valuePtr()[i] * vector[glob_row];
            }
        }
    }

#ifdef MPI_BUILD
    parallel_utils::reduce_vector(result, product);
#endif
}

}

#endif