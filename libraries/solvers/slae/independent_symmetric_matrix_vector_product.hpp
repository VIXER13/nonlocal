#ifndef NONLOCAL_SLAE_INDEPENDENT_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP
#define NONLOCAL_SLAE_INDEPENDENT_SYMMETRIC_MATRIX_VECTOR_PRODUCT_HPP

#include "iterative_solver_base.hpp"

#include "uniform_ranges.hpp"

namespace nonlocal::slae {

template<class T, class I>
class independent_symmetric_matrix_vector_product : public iterative_solver_base<T, I> {
    using _base = iterative_solver_base<T, I>;

    std::vector<std::ranges::iota_view<size_t, size_t>> _threads_ranges;
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
    for(const size_t i : std::ranges::iota_view{1u, size_t(_threaded_product.cols())})
        _threaded_product.block(_threads_ranges[i].front(), 0, result.size() - _threads_ranges[i].front(), 1) += 
        _threaded_product.block(_threads_ranges[i].front(), i, result.size() - _threads_ranges[i].front(), 1);
    result = _threaded_product.col(0);
}

template<class T, class I>
void independent_symmetric_matrix_vector_product<T, I>::matrix_vector_product(
    Eigen::Matrix<T, Eigen::Dynamic, 1>& product,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) const {
#pragma omp parallel num_threads(threads_count())
{
#ifdef _OPENMP
    const I thread = omp_get_thread_num();
#else
    const I thread = 0;
#endif
    _threaded_product.col(thread).setZero();
    for(const I row : _threads_ranges[thread]) {
        const I ind = matrix().outerIndexPtr()[row];
        _threaded_product(row, thread) += matrix().valuePtr()[ind] * vector[matrix().innerIndexPtr()[ind]];
        for(const I i : std::ranges::iota_view{ind + 1, matrix().outerIndexPtr()[row + 1]}) {
            _threaded_product(row, thread) += matrix().valuePtr()[i] * vector[matrix().innerIndexPtr()[i]];
            _threaded_product(matrix().innerIndexPtr()[i], thread) += matrix().valuePtr()[i] * vector[row];
        }
    }
}
    reduction(product);
}

template<class T, class I>
void independent_symmetric_matrix_vector_product<T, I>::set_threads_count(const size_t threads_count) {
    _base::set_threads_count(threads_count > matrix().rows() ? matrix().rows() : threads_count);
    _threads_ranges = parallel_utils::uniform_ranges(matrix(), _base::threads_count());
    _threaded_product.resize(matrix().rows(), _base::threads_count());
}

}

#endif