#ifndef NONLOCAL_SLAE_INDEPENDENT_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP
#define NONLOCAL_SLAE_INDEPENDENT_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP

#include "iterative_solver_base.hpp"

#include "OMP_utils.hpp"
#include "uniform_ranges.hpp"

namespace nonlocal::slae {

template<class T, class I>
class independent_symmetric_matrix_vector_product : public iterative_solver_base<T, I> {
    using _base = iterative_solver_base<T, I>;

    parallel_utils::OMP_ranges _thread_rows;
    mutable Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _threaded_product;

protected:
    void reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& result) const;
    void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& product,
                               const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const;

public:
    using _base::matrix;
    using _base::threads_count;

    explicit independent_symmetric_matrix_vector_product(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);

    void set_threads_count(const size_t threads_count) override;
};

template<class T, class I>
independent_symmetric_matrix_vector_product<T, I>::independent_symmetric_matrix_vector_product(
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix)
    : _base{matrix} {
        set_threads_count(threads_count());
}

template<class T, class I>
void independent_symmetric_matrix_vector_product<T, I>::reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& result) const {
    const size_t shift = _base::process_rows().front();
    for(const size_t i : std::ranges::iota_view{1u, _thread_rows.size()}) {
        const size_t row = shift + _thread_rows.get(i).front();
        _threaded_product.block(row, 0, result.size() - row, 1) += 
        _threaded_product.block(row, i, result.size() - row, 1);
    }
    parallel_utils::reduce_vector(result, _threaded_product.col(0));
}

template<class T, class I>
void independent_symmetric_matrix_vector_product<T, I>::matrix_vector_product(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& product,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const {
    const size_t shift = _base::process_rows().front();
#pragma omp parallel num_threads(threads_count())
{
#ifdef _OPENMP
    const I thread = omp_get_thread_num();
#else
    const I thread = 0;
#endif
    _threaded_product.col(thread).setZero();
    for(const I row : _thread_rows.get(thread)) {
        const I ind = matrix().outerIndexPtr()[row];
        const I glob_row = row + shift;
        _threaded_product(glob_row, thread) += matrix().valuePtr()[ind] * vector[matrix().innerIndexPtr()[ind]];
        for(const I i : std::ranges::iota_view{ind + 1, matrix().outerIndexPtr()[row + 1]}) {
            _threaded_product(glob_row, thread) += matrix().valuePtr()[i] * vector[matrix().innerIndexPtr()[i]];
            _threaded_product(matrix().innerIndexPtr()[i], thread) += matrix().valuePtr()[i] * vector[glob_row];
        }
    }
}
    reduction(product);
}

template<class T, class I>
void independent_symmetric_matrix_vector_product<T, I>::set_threads_count(const size_t threads_count) {
    _base::set_threads_count(threads_count > matrix().rows() ? matrix().rows() : threads_count);
    _thread_rows = parallel_utils::OMP_ranges{ parallel_utils::uniform_ranges(matrix(), _base::threads_count()) };
    _threaded_product.resize(matrix().cols(), _base::threads_count());
}

}

#endif