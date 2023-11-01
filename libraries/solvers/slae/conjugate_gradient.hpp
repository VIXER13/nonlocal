#ifndef NONLOCAL_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_CONJUGATE_GRADIENT_HPP

#include "OMP_utils.hpp"
#include "uniform_ranges.hpp"
#include "unrelated_rows.hpp"

#include <Eigen/Sparse>
#include <iostream>
#include <optional>

namespace nonlocal::slae {

template<class T>
struct conjugate_gradient_parameters final {
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    uintmax_t max_iterations = 10000;
    int threads_count = parallel_utils::threads_count();
};

template<class T, class I, class Preconditioner = Eigen::IdentityPreconditioner>
class conjugate_gradient final {
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& _A;
    conjugate_gradient_parameters<T> _parameters = {};
    //parallel_utils::unrelated_rows<I> _unrelated;
    Preconditioner _preconditioner;
    std::vector<std::ranges::iota_view<size_t, size_t>> _threads_ranges;
    mutable Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _threaded_z;
    mutable uintmax_t _iteration = 0;
    mutable T _residual = 0;

    void reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& z) const;
    void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z,
                               const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const;

public:
    explicit conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                const conjugate_gradient_parameters<T>& parameters = {});

    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix() const noexcept;
    const Preconditioner& preconditioner() const noexcept;
    Preconditioner& preconditioner() noexcept;

    T tolerance() const noexcept;
    uintmax_t max_iterations() const noexcept;
    int threads_count() const noexcept;

    T residual() const noexcept;
    uintmax_t iterations() const noexcept;

    void set_tolerance(const T tolerance) noexcept;
    void set_max_iterations(const uintmax_t max_iterations) noexcept;
    void set_threads_count(const int threads_count);

    Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                                              const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0 = std::nullopt) const;
};

template<class T, class I, class Preconditioner>
conjugate_gradient<T, I, Preconditioner>::conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                             const conjugate_gradient_parameters<T>& parameters)
    : _A{A}
    , _parameters{parameters}
    //, _unrelated{A} 
    {
        set_threads_count(_parameters.threads_count);
    }

template<class T, class I, class Preconditioner>
const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& conjugate_gradient<T, I, Preconditioner>::matrix() const noexcept {
    return _A;
}

template<class T, class I, class Preconditioner>
const Preconditioner& conjugate_gradient<T, I, Preconditioner>::preconditioner() const noexcept {
    return _preconditioner;
}

template<class T, class I, class Preconditioner>
Preconditioner& conjugate_gradient<T, I, Preconditioner>::preconditioner() noexcept {
    return _preconditioner;
}

template<class T, class I, class Preconditioner>
T conjugate_gradient<T, I, Preconditioner>::tolerance() const noexcept {
    return _parameters.tolerance;
}

template<class T, class I, class Preconditioner>
uintmax_t conjugate_gradient<T, I, Preconditioner>::max_iterations() const noexcept {
    return _parameters.max_iterations;
}

template<class T, class I, class Preconditioner>
int conjugate_gradient<T, I, Preconditioner>::threads_count() const noexcept {
    return _parameters.threads_count;
}

template<class T, class I, class Preconditioner>
T conjugate_gradient<T, I, Preconditioner>::residual() const noexcept {
    return _residual;
}

template<class T, class I, class Preconditioner>
uintmax_t conjugate_gradient<T, I, Preconditioner>::iterations() const noexcept {
    return _iteration;
}

template<class T, class I, class Preconditioner>
void conjugate_gradient<T, I, Preconditioner>::set_tolerance(const T tolerance) noexcept {
    _parameters.tolerance = tolerance;
}

template<class T, class I, class Preconditioner>
void conjugate_gradient<T, I, Preconditioner>::set_max_iterations(const uintmax_t max_iterations) noexcept {
    _parameters.max_iterations = max_iterations;
}

template<class T, class I, class Preconditioner>
void conjugate_gradient<T, I, Preconditioner>::set_threads_count(const int threads_count) {
    _parameters.threads_count = threads_count;
    if (_parameters.threads_count > _A.rows())
        _parameters.threads_count = _A.rows();
    _threaded_z.resize(_A.rows(), _parameters.threads_count);
    _threads_ranges = parallel_utils::uniform_ranges(_A, _parameters.threads_count);
}

template<class T, class I, class Preconditioner>
void conjugate_gradient<T, I, Preconditioner>::reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& z) const {
    for(const size_t i : std::ranges::iota_view{1u, size_t(_threaded_z.cols())})
        _threaded_z.block(_threads_ranges[i].front(), 0, z.size() - _threads_ranges[i].front(), 1) += 
        _threaded_z.block(_threads_ranges[i].front(), i, z.size() - _threads_ranges[i].front(), 1);
    z = _threaded_z.col(0);
}

template<class T, class I, class Preconditioner>
void conjugate_gradient<T, I, Preconditioner>::matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z,
                                                                     const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const {
#pragma omp parallel default(none) shared(p) num_threads(_parameters.threads_count)
{
#ifdef _OPENMP
    const I thread = omp_get_thread_num();
#else
    const I thread = 0;
#endif
    _threaded_z.col(thread).setZero();
    for(const I row : _threads_ranges[thread]) {
        const I ind = _A.outerIndexPtr()[row];
        _threaded_z(row, thread) += _A.valuePtr()[ind] * p[_A.innerIndexPtr()[ind]];
        for(const I i : std::ranges::iota_view{ind + 1, _A.outerIndexPtr()[row + 1]}) {
            _threaded_z(row, thread) += _A.valuePtr()[i] * p[_A.innerIndexPtr()[i]];
            _threaded_z(_A.innerIndexPtr()[i], thread) += _A.valuePtr()[i] * p[row];
        }
    }
}
    reduction(z);
}

// template<class T, class I, class Preconditioner>
// void conjugate_gradient<T, I, Preconditioner>::matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z,
//                                                                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const {
//     z.setZero();
//     for(const size_t shift_index : std::ranges::iota_view<size_t, size_t>{0u, _unrelated.shifts.size() - 1}) {
//         const I rows_count = _unrelated.shifts[shift_index + 1] - _unrelated.shifts[shift_index];
//         const int threads_count = std::min(I(_parameters.threads_count), rows_count);
// #pragma omp parallel for num_threads(threads_count)
//         for(size_t shift = _unrelated.shifts[shift_index]; shift < _unrelated.shifts[shift_index + 1]; ++shift) {
//             const I row = _unrelated.rows[shift];
//             const I ind = _A.outerIndexPtr()[row];
//             z[row] += _A.valuePtr()[ind] * p[_A.innerIndexPtr()[ind]];
//             for(const I i : std::ranges::iota_view{ind + 1, _A.outerIndexPtr()[row + 1]}) {
//                 z[row] += _A.valuePtr()[i] * p[_A.innerIndexPtr()[i]];
//                 z[_A.innerIndexPtr()[i]] += _A.valuePtr()[i] * p[row];
//             }
//         }
//     }
// }

template<class T, class I, class Preconditioner>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient<T, I, Preconditioner>::solve(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0) const {
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b_full = b;
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = x0.template value_or(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size()));
    Eigen::Matrix<T, Eigen::Dynamic, 1> r = b_full - _A.template selfadjointView<Eigen::Upper>() * x;
    Eigen::Matrix<T, Eigen::Dynamic, 1> p = r;
    Eigen::Matrix<T, Eigen::Dynamic, 1> z = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size());
    T r_squared_norm = r.squaredNorm();
    const T b_norm = b_full.norm();
    _iteration = 0;
    _residual = std::sqrt(r_squared_norm) / b_norm;
    while(_iteration < _parameters.max_iterations && _residual > _parameters.tolerance) {
        matrix_vector_product(z, p);
        const T nu = r_squared_norm / p.dot(z);
        x += nu * p;
        r -= nu * z;
        const T r_squared_norm_prev = std::exchange(r_squared_norm, r.squaredNorm()),
                mu = r_squared_norm / r_squared_norm_prev;
        p = r + mu * p;
        ++_iteration;
        _residual = std::sqrt(r_squared_norm) / b_norm;
    }
    return x;
}

}

#endif