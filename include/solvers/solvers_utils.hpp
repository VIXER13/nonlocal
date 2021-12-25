#ifndef SOLVERS_UTILS_HPP
#define SOLVERS_UTILS_HPP

#include "solvers_constants.hpp"
#include <eigen3/Eigen/Sparse>
#include <algorithm>
#include <ranges>

namespace nonlocal::utils {

template<class T, class I>
void prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
    for(const size_t i : std::views::iota(size_t{0}, size_t(K.rows())))
        K.outerIndexPtr()[i+1] += K.outerIndexPtr()[i];
    K.data().resize(K.outerIndexPtr()[K.rows()]);
    for(const size_t i : std::views::iota(size_t{0}, size_t(K.outerIndexPtr()[K.rows()]))) {
        K.innerIndexPtr()[i] = K.cols()-1;
        K.valuePtr()[i] = T{0};
    }
}

template<class T, class I>
void sort_indices(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
#pragma omp parallel for default(none) shared(K) schedule(dynamic)
    for(size_t i = 0; i < K.rows(); ++i)
        std::sort(&K.innerIndexPtr()[K.outerIndexPtr()[i]], &K.innerIndexPtr()[K.outerIndexPtr()[i+1]]);
}

}

#endif