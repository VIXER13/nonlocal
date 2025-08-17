#pragma once

#include <Eigen/Sparse>

#include <algorithm>
#include <cstring>
#include <ranges>
#include <variant>

namespace nonlocal::utils {

using nodes_sequence = std::variant<
    std::ranges::iota_view<size_t, size_t>,
    std::vector<size_t>
>;

template<class Callback>
void iterate(const nodes_sequence& sequence, Callback&& callback, const int threads = 1) {
    std::visit([&callback, threads](const auto& sequence) {
#pragma omp parallel for firstprivate(callback) num_threads(threads)
        for(size_t i = 0; i < sequence.size(); ++i)
            callback(sequence[i]); 
    }, sequence);
}

template<class T, class I>
void accumulate_shifts(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
    for(const size_t i : std::ranges::iota_view{0u, size_t(matrix.rows())})
        matrix.outerIndexPtr()[i + 1] += matrix.outerIndexPtr()[i];
}

// Warning: the outer indices shall be accumulated before calling!
template<class T, class I>
void allocate_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
    matrix.data().resize(matrix.outerIndexPtr()[matrix.rows()]);
    std::memset(matrix.innerIndexPtr(), '\0', sizeof(I) * size_t(matrix.nonZeros()));
    std::memset(matrix.valuePtr(), '\0', sizeof(T) * size_t(matrix.nonZeros()));
}

template<class T, class I>
void sort_indices(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
#pragma omp parallel for default(none) shared(matrix) schedule(dynamic)
    for(size_t i = 0; i < size_t(matrix.rows()); ++i)
        std::sort(&matrix.innerIndexPtr()[matrix.outerIndexPtr()[i]], &matrix.innerIndexPtr()[matrix.outerIndexPtr()[i + 1]]);
}

}