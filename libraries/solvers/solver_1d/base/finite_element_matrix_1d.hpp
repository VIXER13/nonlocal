#pragma once

#include <Eigen/Sparse>

namespace nonlocal {

template<class T, class I>
class finite_element_matrix_1d final {
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> _matrix_inner;
    std::array<std::unordered_map<size_t, T>, 2> _matrix_bound;

public:
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& inner() noexcept { return _matrix_inner; }
    std::array<std::unordered_map<size_t, T>, 2>& bound() noexcept { return _matrix_bound; }
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& inner() const noexcept { return _matrix_inner; }
    const std::array<std::unordered_map<size_t, T>, 2>& bound() const noexcept { return _matrix_bound; }

    void clear() {
        inner() = {};
        bound() = {};
    }
};

}