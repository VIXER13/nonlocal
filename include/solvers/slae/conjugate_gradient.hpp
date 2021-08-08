#ifndef NONLOCAL_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_CONJUGATE_GRADIENT_HPP

#include <cmath>
#undef I // for new version GCC, when use I macros
#include <eigen3/Eigen/Sparse>
#include "MPI_utils.hpp"

namespace nonlocal::slae {

template<class T>
struct conjugate_gradient_parameters final {
    T eps = 2.220446e-16;
    uintmax_t max_iterations = 10000;
};

class _conjugate_gradient final {
    explicit _conjugate_gradient() noexcept = default;

    template<class T>
    static Eigen::Matrix<T, Eigen::Dynamic, 1> calc_b_full(const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, const MPI_utils::MPI_ranges& ranges) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> b_full(ranges.ranges().back().back());
        for(size_t i = ranges.range().front(), j = 0; j < b.size(); ++i, ++j)
            b_full[i] = b[j];
        return MPI_utils::all_to_all<T>(b_full, ranges);
    }

    template<class T>
    static void reduction(Eigen::Matrix<T, Eigen::Dynamic, 1>& Ap,
                          Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& threadedAp) {
        for(size_t i = 1; i < threadedAp.cols(); ++i)
            threadedAp.col(0) += threadedAp.col(i);
#if MPI_USE
        Ap.setZero();
        MPI_Allreduce(threadedAp.col(0).data(), Ap.data(), Ap.size(), std::is_same_v<T, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        Ap = threadedAp.col(0);
#endif
    }

    template<class T, class I>
    static void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& Ap,
                                      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& threadedAp,
                                      const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
                                      const size_t shift) {
        threadedAp.setZero();
#pragma omp parallel default(none) shared(Ap, threadedAp, A, p, shift)
        {
            const int thread = omp_get_thread_num();
#pragma omp for schedule(dynamic)
            for(I row = 0; row < A.rows(); ++row) {
                for(I i = A.outerIndexPtr()[row]; i < A.outerIndexPtr()[row+1]; ++i)
                    threadedAp(row + shift, thread) += A.valuePtr()[i] * p[A.innerIndexPtr()[i]];
                for(I i = A.outerIndexPtr()[row]+1; i < A.outerIndexPtr()[row+1]; ++i)
                    threadedAp(A.innerIndexPtr()[i], thread) += A.valuePtr()[i] * p[row + shift];
            }
        }
        reduction(Ap, threadedAp);
    }

public:
    template<class T, class I>
    friend Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                                  const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                                                                  const conjugate_gradient_parameters<T>& parameters,
                                                                  const MPI_utils::MPI_ranges& ranges);
};

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                                                       const conjugate_gradient_parameters<T>& parameters,
                                                       const MPI_utils::MPI_ranges& ranges) {
#if MPI_USE
    const Eigen::Matrix<T, Eigen::Dynamic, 1> b_full = _conjugate_gradient::calc_b_full(b, ranges);
#else
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b_full = b;
#endif

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> threadedAp(b_full.size(), omp_get_max_threads());
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b_full.size()),
                                        r = b_full, p = b_full, z = r;
    uintmax_t iterations = 0;
    const T b_norm = b_full.norm();
    T r_squaredNorm = r.squaredNorm();
    const size_t shift = ranges.range().front();
    for(; iterations < parameters.max_iterations && std::sqrt(r_squaredNorm) / b_norm > parameters.eps; ++iterations) {
        _conjugate_gradient::matrix_vector_product(z, threadedAp, A, p, shift);
        const T nu = r_squaredNorm / p.dot(z);
        x += nu * p;
        r -= nu * z;
        const T r_squaredNorm_prev = std::exchange(r_squaredNorm, r.squaredNorm()),
        mu = r_squaredNorm / r_squaredNorm_prev;
        p = r + mu * p;
    }
    std::cout << "iterations = " << iterations << std::endl;
    return x;
}

//template<class T, class I>
//Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient_jacobi(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
//                                                              const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
//    static constexpr T eps = 2.220446e-16;
//    static constexpr uintmax_t maxiter = 10000;
//
//    Eigen::MatrixXd threadedAp(A.rows(), omp_get_max_threads());
//    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b.size()),
//                                        M = A.diagonal(), r = b, p = b, y = r, z = r;
//    for(size_t i = 0; i < b.size(); ++i) {
//        M[i] = T{1} / M[i];
//        p[i] *= M[i];
//        y[i] *= M[i];
//    }
//
//    uintmax_t iterations = 0;
//    const T b_norm = b.norm();
//    for(; iterations < maxiter && r.norm() / b_norm > eps; ++iterations) {
//        _conjugate_gradient::matrix_vector_product(z, threadedAp, A, p);
//        const T y_dot_r = y.dot(r),
//                nu = y_dot_r / p.dot(z);
//        x += nu * p;
//        r -= nu * z;
//        for(size_t i = 0; i < b.size(); ++i)
//            y[i] = M[i] * r[i];
//        const T mu = y.dot(r) / y_dot_r;
//        p *= mu;
//        p += y;
//    }
//    std::cout << "iterations = " << iterations << std::endl;
//
//    return x;
//}

}

#endif