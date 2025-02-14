#pragma once

#include "finite_element_matrix.hpp"

#include <vector>

namespace nonlocal {

template<class T, class I>
class matrix_separator_base {
    finite_element_matrix<T, I>& _matrix;
    const std::vector<bool>& _is_inner;
    const size_t _node_shift;
    const bool _is_symmetric;

protected:
    explicit matrix_separator_base(finite_element_matrix<T, I>& matrix, const std::vector<bool>& is_inner,
                                   const size_t node_shift, const bool is_symmetric);

    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix(const matrix_part part) noexcept;
    size_t node_shift() const noexcept;
    bool is_symmetric() const noexcept;

    matrix_part part(const size_t row, const size_t col);

public:
    virtual ~matrix_separator_base() noexcept = default;
};

template<class T, class I>
matrix_separator_base<T, I>::matrix_separator_base(finite_element_matrix<T, I>& matrix, const std::vector<bool>& is_inner,
                                                   const size_t node_shift, const bool is_symmetric)
    : _matrix{matrix}, _is_inner{is_inner}, _node_shift{node_shift}, _is_symmetric{is_symmetric} {}

template<class T, class I>
Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix_separator_base<T, I>::matrix(const matrix_part part) noexcept {
    return _matrix[part];
}

template<class T, class I>
size_t matrix_separator_base<T, I>::node_shift() const noexcept {
    return _node_shift;
}

template<class T, class I>
bool matrix_separator_base<T, I>::is_symmetric() const noexcept {
    return _is_symmetric;
}

template<class T, class I>
matrix_part matrix_separator_base<T, I>::part(const size_t row, const size_t col) {
    if (_is_inner[col]) {
        if ((!is_symmetric() || row <= col) && _is_inner[row])
            return matrix_part::INNER;
    } else if (row != col)
        return matrix_part::BOUND;
    return matrix_part::NO;
}

}