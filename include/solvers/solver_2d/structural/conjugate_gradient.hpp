#ifndef NONLOCAL_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_CONJUGATE_GRADIENT_HPP

#include <eigen3/Eigen/Sparse>
#include <cmath>
#undef I // for new version GCC, when use I macros

namespace nonlocal {

template<class T, class I>
void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& Ap,
                           Eigen::MatrixXd& threadedAp,
                           const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                           const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) {
    threadedAp.setZero();
#pragma omp parallel default(none) shared(Ap, threadedAp, A, p)
{
    const int thread = omp_get_thread_num();
#pragma omp for schedule(dynamic)
    for(I row = 0; row < A.rows(); ++row) {
        T sum = T{0};
        for(I i = A.outerIndexPtr()[row]; i < A.outerIndexPtr()[row+1]; ++i)
            sum += A.valuePtr()[i] * p[A.innerIndexPtr()[i]];
        Ap[row] = sum;
        for(I i = A.outerIndexPtr()[row]+1; i < A.outerIndexPtr()[row+1]; ++i)
            threadedAp(A.innerIndexPtr()[i], thread) += A.valuePtr()[i] * p[row];
    }
}
    for(I col = 0; col < threadedAp.cols(); ++col)
        Ap += threadedAp.col(col);
}

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
    static constexpr T eps = 2.220446e-16;
    static constexpr uintmax_t maxiter = 10000;

    Eigen::MatrixXd threadedAp(A.rows(), omp_get_max_threads());
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b.size()),
                                        r = b, p = b, z = r;
    uintmax_t iterations = 0;
    const T b_norm = b.norm();
    T r_squaredNorm = r.squaredNorm();
    for(; iterations < maxiter && std::sqrt(r_squaredNorm) / b_norm > eps; ++iterations) {
        matrix_vector_product(z, threadedAp, A, p);
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

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient_jacobi(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                              const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
    static constexpr T eps = 2.220446e-16;
    static constexpr uintmax_t maxiter = 10000;

    Eigen::MatrixXd threadedAp(A.rows(), omp_get_max_threads());
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b.size()),
                                        M = A.diagonal(), r = b, p = b, y = r, z = r;
    for(size_t i = 0; i < b.size(); ++i) {
        M[i] = T{1} / M[i];
        p[i] *= M[i];
        y[i] *= M[i];
    }

    uintmax_t iterations = 0;
    const T b_norm = b.norm();
    for(; iterations < maxiter && r.norm() / b_norm > eps; ++iterations) {
        matrix_vector_product(z, threadedAp, A, p);
        const T y_dot_r = y.dot(r),
                nu = y_dot_r / p.dot(z);
        x += nu * p;
        r -= nu * z;
        for(size_t i = 0; i < b.size(); ++i)
            y[i] = M[i] * r[i];
        const T mu = y.dot(r) / y_dot_r;
        p *= mu;
        p += y;
    }
    std::cout << "iterations = " << iterations << std::endl;

    return x;
}

}

#endif