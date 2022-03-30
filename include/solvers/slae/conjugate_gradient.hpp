#ifndef NONLOCAL_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_CONJUGATE_GRADIENT_HPP

#include <eigen3/Eigen/Sparse>
#include <omp.h>
#include <optional>
#include <iostream>

namespace nonlocal::slae {

template<class T>
struct conjugate_gradient_parameters final {
    T tolerance = 1e-10;
    uintmax_t max_iterations = 10000;
    size_t threads_count = omp_get_max_threads();
};

template<class T, class I>
class conjugate_gradient final {
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& _A;
    conjugate_gradient_parameters<T> _parameters = {};
    std::vector<std::array<size_t, 2>> _threads_ranges;
    mutable Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _threaded_z;
    mutable uintmax_t _iteration = 0;
    mutable T _residual = 0;

    static std::vector<std::array<size_t, 2>> distribute_rows(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                              const size_t threads_count);

    void reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& z) const;
    void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z,
                               const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const;

public:
    explicit conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                const conjugate_gradient_parameters<T>& parameters = {});

    T tolerance() const noexcept;
    uintmax_t max_iterations() const noexcept;

    T residual() const noexcept;
    uintmax_t iterations() const noexcept;

    void set_tolerance(const T tolerance) noexcept;
    void set_max_iterations(const uintmax_t max_iterations) noexcept;

    Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                                              const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0 = std::nullopt) const;
};

template<class T, class I>
conjugate_gradient<T, I>::conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                             const conjugate_gradient_parameters<T>& parameters)
    : _A{A}
    , _parameters{parameters}
    , _threads_ranges{distribute_rows(A, _parameters.threads_count)}
    , _threaded_z(A.cols(), _parameters.threads_count) {}

template<class T, class I>
std::vector<std::array<size_t, 2>> conjugate_gradient<T, I>::distribute_rows(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                                             const size_t threads_count) {
    if (threads_count < 0)
        throw std::logic_error{"Threads count must be greater than 0."};
    const size_t mean_count = A.nonZeros() / threads_count;
    std::vector<std::array<size_t, 2>> threads_ranges(threads_count);
    threads_ranges.front().front() = 0;
    for(size_t row = 0, thread_sum = 0, curr_thread = 0; row < size_t(A.rows()); ++row) {
        thread_sum += A.outerIndexPtr()[row + 1] - A.outerIndexPtr()[row];
        if (thread_sum > mean_count) {
            thread_sum = 0;
            threads_ranges[curr_thread].back() = row;
            ++curr_thread;
            threads_ranges[curr_thread].front() = row;
        }
    }
    threads_ranges.back().back() = A.rows();
    return threads_ranges;
}

template<class T, class I>
T conjugate_gradient<T, I>::tolerance() const noexcept {
    return _parameters.tolerance;
}

template<class T, class I>
uintmax_t conjugate_gradient<T, I>::max_iterations() const noexcept {
    return _parameters.max_iterations;
}

template<class T, class I>
T conjugate_gradient<T, I>::residual() const noexcept {
    return _residual;
}

template<class T, class I>
uintmax_t conjugate_gradient<T, I>::iterations() const noexcept {
    return _iteration;
}

template<class T, class I>
void conjugate_gradient<T, I>::set_tolerance(const T tolerance) noexcept {
    _parameters.tolerance = tolerance;
}

template<class T, class I>
void conjugate_gradient<T, I>::set_max_iterations(const uintmax_t max_iterations) noexcept {
    _parameters.max_iterations = max_iterations;
}

template<class T, class I>
void conjugate_gradient<T, I>::reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& z) const {
    for(size_t i = 1; i < _threaded_z.cols(); ++i)
        _threaded_z.col(0) += _threaded_z.col(i);
    z = _threaded_z.col(0);
}

template<class T, class I>
void conjugate_gradient<T, I>::matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z,
                                                     const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const {
#pragma omp parallel default(none) shared(p)
{
    const I thread = omp_get_thread_num();
    _threaded_z.col(thread).setZero();
    for(I row = _threads_ranges[thread].front(); row < I(_threads_ranges[thread].back()); ++row) {
        const I ind = _A.outerIndexPtr()[row];
        _threaded_z(row, thread) += _A.valuePtr()[ind] * p[_A.innerIndexPtr()[ind]];
        for(I i = ind+1; i < _A.outerIndexPtr()[row+1]; ++i) {
            _threaded_z(row, thread) += _A.valuePtr()[i] * p[_A.innerIndexPtr()[i]];
            _threaded_z(_A.innerIndexPtr()[i], thread) += _A.valuePtr()[i] * p[row];
        }
    }
}
    reduction(z);
}

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient<T, I>::solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                                                                    const std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>& x0) const {
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b_full = b;
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = x0.template value_or(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size()));
    Eigen::Matrix<T, Eigen::Dynamic, 1> r = b_full - _A.template selfadjointView<Eigen::Upper>() * x;
    Eigen::Matrix<T, Eigen::Dynamic, 1> p = r;
    Eigen::Matrix<T, Eigen::Dynamic, 1> z = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size());
    T r_squaredNorm = r.squaredNorm();
    const T b_norm = b_full.norm();
    _iteration = 0;
    _residual = std::sqrt(r_squaredNorm) / b_norm;
    std::cout << _iteration << " " << _residual << std::endl;
    while(_iteration < _parameters.max_iterations && _residual > _parameters.tolerance) {
        matrix_vector_product(z, p);
        const T nu = r_squaredNorm / p.dot(z);
        x += nu * p;
        r -= nu * z;
        const T r_squaredNorm_prev = std::exchange(r_squaredNorm, r.squaredNorm()),
                mu = r_squaredNorm / r_squaredNorm_prev;
        p = r + mu * p;
        ++_iteration;
        _residual = std::sqrt(r_squaredNorm) / b_norm;
        //std::cout << _iteration << " " << _residual << std::endl;
    }
    std::cout << _iteration << " " << _residual << std::endl;
    return x;
}

}

#endif