#ifndef NONLOCAL_MESH_RUNNER_TYPES_HPP
#define NONLOCAL_MESH_RUNNER_TYPES_HPP

#include <Eigen/Sparse>

#include <array>

namespace nonlocal {

enum class matrix_part : size_t {
    INNER,
    BOUND,
    NO
};

template<class T, class I>
class matrix_parts final {
    std::array<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, 2> _part;

public:
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& operator[](const matrix_part part) {
        return _part[size_t(part)];
    }

    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& operator[](const matrix_part part) const {
        return _part[size_t(part)];
    }

    void clear() {
        _part.front() = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{};
        _part.back() = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{};
    }
};

}

#endif