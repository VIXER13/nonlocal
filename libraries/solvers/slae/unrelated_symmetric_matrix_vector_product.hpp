#ifndef NONLOCAL_SLAE_UNRELATED_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP
#define NONLOCAL_SLAE_UNRELATED_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP

#include "iterative_solver_base.hpp"

#include "unrelated_rows.hpp"

namespace nonlocal::slae {

template<class T, class I>
class unrelated_symmetric_matrix_vector_product : public iterative_solver_base<T, I> {
    using _base = iterative_solver_base<T, I>;

    parallel_utils::unrelated_rows<I> _unrelated;

protected:
    void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& product,
                               const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const;

public:
    using _base::matrix;
    using _base::threads_count;

    explicit unrelated_symmetric_matrix_vector_product(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);
};

template<class T, class I>
unrelated_symmetric_matrix_vector_product<T, I>::unrelated_symmetric_matrix_vector_product(
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
    : iterative_solver_base<T, I>{matrix}
    , _unrelated{matrix} {}

template<class T, class I>
void unrelated_symmetric_matrix_vector_product<T, I>::matrix_vector_product(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& product,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const {
    product.setZero();
    for(const size_t shift_index : std::ranges::iota_view<size_t, size_t>{0u, _unrelated.shifts.size() - 1}) {
        const I rows_count = _unrelated.shifts[shift_index + 1] - _unrelated.shifts[shift_index];
#pragma omp parallel for num_threads(std::min(I(threads_count()), rows_count))
        for(size_t shift = _unrelated.shifts[shift_index]; shift < _unrelated.shifts[shift_index + 1]; ++shift) {
            const I row = _unrelated.rows[shift];
            const I ind = matrix().outerIndexPtr()[row];
            product[row] += matrix().valuePtr()[ind] * vector[matrix().innerIndexPtr()[ind]];
            for(const I i : std::ranges::iota_view{ind + 1, matrix().outerIndexPtr()[row + 1]}) {
                product[row] += matrix().valuePtr()[i] * vector[matrix().innerIndexPtr()[i]];
                product[matrix().innerIndexPtr()[i]] += matrix().valuePtr()[i] * vector[row];
            }
        }
    }
}

}

#endif