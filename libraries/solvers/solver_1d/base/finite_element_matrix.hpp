#pragma once

#include <Eigen/Sparse>

#include <array>
#include <unordered_map>

namespace nonlocal {

template<class T, class I>
struct finite_element_matrix_1d final {
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> inner;
    std::array<std::unordered_map<size_t, T>, 2> bound;

    void clear() {
        inner = {};
        bound = {};
    }
};

}