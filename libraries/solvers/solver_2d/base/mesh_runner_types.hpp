#ifndef NONLOCAL_MESH_RUNNER_TYPES_HPP
#define NONLOCAL_MESH_RUNNER_TYPES_HPP

#include <eigen3/Eigen/Sparse>

#include <array>

namespace nonlocal {

enum class matrix_part : size_t {
    INNER,
    BOUND,
    NO
};

template<class T, class I>
using matrix_parts_t = std::array<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, 2>;

}

#endif