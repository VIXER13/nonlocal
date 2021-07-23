#ifndef NONLOCAL_CONJUGATE_GRADIENT_HPP
#define NONLOCAL_CONJUGATE_GRADIENT_HPP

#include <eigen3/Eigen/Sparse>
#include <cmath>
#undef I // for new version GCC, when use I macros

namespace nonlocal {

template<class T, class I>
void matrix_vector_product(Eigen::Matrix<T, Eigen::Dynamic, 1>& Az,
                           Eigen::MatrixXd& threadedAz,
                           const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                           const Eigen::Matrix<T, Eigen::Dynamic, 1>& z) {
    threadedAz.setZero();
#pragma omp parallel default(none) shared(Az, threadedAz, A, z)
{
#pragma omp for schedule(dynamic)
    for(I row = 0; row < A.rows(); ++row) {
        const int thread = omp_get_thread_num();
        T sum = T{0};
        for(I i = A.outerIndexPtr()[row]; i < A.outerIndexPtr()[row+1]; ++i)
            sum += A.valuePtr()[i] * z[A.innerIndexPtr()[i]];
        Az[row] = sum;
        for(I i = A.outerIndexPtr()[row]+1; i < A.outerIndexPtr()[row+1]; ++i)
            threadedAz(A.innerIndexPtr()[i], thread) += A.valuePtr()[i] * z[row];
    }

#pragma omp for
    for(I row = 0; row < threadedAz.rows(); ++row)
        for(I col = 0; col < threadedAz.cols(); ++col)
            Az[row] += threadedAz(row, col);
}
}

template<class T, class I>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugate_gradient(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& A,
                                                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
    static constexpr T eps = 2.220446e-16;
    static constexpr uintmax_t maxiter = 10000;

    Eigen::Matrix<T, Eigen::Dynamic, 1> r = b, z = r, Az = z, x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(b.size());
    Eigen::MatrixXd threadedAz(A.rows(), omp_get_max_threads());

    uintmax_t iterations = 0;
    const T b_norm = b.norm();
    T r_squaredNorm = r.squaredNorm();
    for(; iterations < maxiter && std::sqrt(r_squaredNorm) / b_norm > eps; ++iterations) {
        matrix_vector_product(Az, threadedAz, A, z);
        const T alpha = r_squaredNorm / Az.dot(z);
        x += alpha * z;
        r -= alpha * Az;
        const T r_squaredNorm_prev = std::exchange(r_squaredNorm, r.squaredNorm()),
                betta = r_squaredNorm / r_squaredNorm_prev;
        z = r + betta * z;
    }
    std::cout << "iterations = " << iterations << std::endl;
    return x;
}

}

#endif