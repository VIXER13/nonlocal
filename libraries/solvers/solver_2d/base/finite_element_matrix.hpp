#pragma once

#include <Eigen/Sparse>

#include <array>

namespace nonlocal {

enum class matrix_part : size_t {
    INNER,
    BOUND,
    NO
};

template<class T, class I>
class finite_element_matrix final {
    std::array<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, 2> _part;

public:
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& inner() noexcept { return _part[size_t(matrix_part::INNER)]; }
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& bound() noexcept { return _part[size_t(matrix_part::BOUND)]; }
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& inner() const noexcept { return _part[size_t(matrix_part::INNER)]; }
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& bound() const noexcept { return _part[size_t(matrix_part::BOUND)]; }

    void clear() {
        _part.front() = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{};
        _part.back() = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{};
    }

    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& operator[](const matrix_part part) { return _part[size_t(part)]; }
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& operator[](const matrix_part part) const { return _part[size_t(part)]; }
};

}