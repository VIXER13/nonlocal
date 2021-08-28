#ifndef NONLOCAL_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_CONJUGATE_GRADIENT_HPP

#include "MPI_utils.hpp"
#include <eigen3/Eigen/Sparse>
#include <cmath>
#undef I // for new version GCC, when use I macros

namespace nonlocal::slae {

template<class T>
struct conjugate_gradient_parameters final {
    T tolerance = 2.220446e-16;
    uintmax_t max_iterations = 10000;
};

template<class T, class I>
class conjugate_gradient final {
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& _A;
    mutable Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _threaded_z;
    const MPI_utils::MPI_ranges _ranges;
    const size_t _shift;
    mutable conjugate_gradient_parameters<T> _curr_parameters{};
    conjugate_gradient_parameters<T> _parameters;

    Eigen::Matrix<T, Eigen::Dynamic, 1> calc_b_full(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) const;

    void reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& z) const;

    void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const;

public:
    explicit conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                const MPI_utils::MPI_ranges& ranges,
                                const conjugate_gradient_parameters<T> parameters = {});

    T tolerance() const;
    uintmax_t max_iterations() const;

    T residual() const;
    uintmax_t iterations() const;

    void set_tolerance(const T tolerance);
    void set_max_iterations(const uintmax_t max_iterations);

    Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) const;
};

template<class T, class I>
conjugate_gradient<T, I>::conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                             const MPI_utils::MPI_ranges& ranges,
                                             const conjugate_gradient_parameters<T> parameters)
    : _A{A}
    , _threaded_z(A.cols(), omp_get_max_threads())
    , _ranges{ranges}
    , _shift{ranges.range().front()}
    , _parameters{parameters} {}

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient<T, I>::calc_b_full(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1> b_full(_ranges.ranges().back().back());
    for(size_t i = _ranges.range().front(), j = 0; j < b.size(); ++i, ++j)
        b_full[i] = b[j];
    return MPI_utils::all_to_all<T>(b_full, _ranges);
}

template<class T, class I>
void conjugate_gradient<T, I>::reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& z) const {
    for(size_t i = 1; i < _threaded_z.cols(); ++i)
        _threaded_z.block(_shift, 0, z.size() - _shift, 1) += _threaded_z.block(_shift, i, z.size() - _shift, 1);
#if MPI_USE
    MPI_Allreduce(_threaded_z.col(0).data(), z.data(), z.size(), std::is_same_v<T, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    z = _threaded_z.col(0);
#endif
}

template<class T, class I>
void conjugate_gradient<T, I>::matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& z, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) const {
#pragma omp parallel default(none) shared(p)
{
    const int thread = omp_get_thread_num();
    _threaded_z.col(thread).setZero();
#pragma omp for schedule(dynamic)
    for(I row = 0; row < _A.rows(); ++row) {
        I glob_row = row + _shift;
        for(I i = _A.outerIndexPtr()[row]; i < _A.outerIndexPtr()[row+1]; ++i)
            _threaded_z(glob_row, thread) += _A.valuePtr()[i] * p[_A.innerIndexPtr()[i]];
        for(I i = _A.outerIndexPtr()[row]+1; i < _A.outerIndexPtr()[row+1]; ++i)
            _threaded_z(_A.innerIndexPtr()[i], thread) += _A.valuePtr()[i] * p[glob_row];
    }
}
    reduction(z);
}

template<class T, class I>
T conjugate_gradient<T, I>::tolerance() const {
    return _parameters.tolerance;
}

template<class T, class I>
uintmax_t conjugate_gradient<T, I>::max_iterations() const {
    return _parameters.max_iterations;
}

template<class T, class I>
T conjugate_gradient<T, I>::residual() const {
    return _curr_parameters.tolerance;
}

template<class T, class I>
uintmax_t conjugate_gradient<T, I>::iterations() const {
    return _curr_parameters.max_iterations;
}

template<class T, class I>
void conjugate_gradient<T, I>::set_tolerance(const T tolerance) {
    _parameters.tolerance = tolerance;
}

template<class T, class I>
void conjugate_gradient<T, I>::set_max_iterations(const uintmax_t max_iterations) {
    _parameters.max_iterations = max_iterations;
}

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient<T, I>::solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) const {
#if MPI_USE
    const Eigen::Matrix<T, Eigen::Dynamic, 1> b_full = calc_b_full(b);
#else
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b_full = b;
#endif

    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size()),
                                        r = b_full, p = b_full, z = r;
    const T b_norm = b_full.norm();
    T r_squaredNorm = r.squaredNorm();
    _curr_parameters = {.tolerance = b_norm, .max_iterations = 0};
    while(_curr_parameters.max_iterations < _parameters.max_iterations && _curr_parameters.tolerance > _parameters.tolerance) {
        matrix_vector_product(z, p);
        const T nu = r_squaredNorm / p.dot(z);
        x += nu * p;
        r -= nu * z;
        const T r_squaredNorm_prev = std::exchange(r_squaredNorm, r.squaredNorm()),
        mu = r_squaredNorm / r_squaredNorm_prev;
        p = r + mu * p;
        ++_curr_parameters.max_iterations;
        _curr_parameters.tolerance = std::sqrt(r_squaredNorm) / b_norm;
    }

    return x;
}

}

#endif