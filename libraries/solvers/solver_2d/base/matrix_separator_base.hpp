#ifndef NONLOCAL_MATRIX_SEPARATOR_HPP
#define NONLOCAL_MATRIX_SEPARATOR_HPP

#include "mesh_runner_types.hpp"

#include <vector>

namespace nonlocal {

template<class T, class I>
class matrix_separator_base {
    matrix_parts_t<T, I>& _matrix;
    const std::vector<bool>& _is_inner;
    const size_t _node_shift;

protected:
    explicit matrix_separator_base(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift);

    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix(const matrix_part part) noexcept;
    size_t node_shift() const noexcept;

    matrix_part part(const size_t row, const size_t col);

public:
    virtual ~matrix_separator_base() noexcept = default;
};

template<class T, class I>
matrix_separator_base<T, I>::matrix_separator_base(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift)
    : _matrix{matrix}, _is_inner{is_inner}, _node_shift{node_shift} {}

template<class T, class I>
Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix_separator_base<T, I>::matrix(const matrix_part part) noexcept {
    return _matrix[size_t(part)];
}

template<class T, class I>
size_t matrix_separator_base<T, I>::node_shift() const noexcept {
    return _node_shift;
}

template<class T, class I>
matrix_part matrix_separator_base<T, I>::part(const size_t row, const size_t col) {
    if (_is_inner[col]) {
        if (row <= col && _is_inner[row])
            return matrix_part::INNER;
    } else if (row != col)
        return matrix_part::BOUND;
    return matrix_part::NO;
}

}

#endif