#pragma once

#include <metamath/linear/linear.hpp>

#include <array>

namespace nonlocal::solver_2d {

enum class matrix_part : size_t {
    INNER,
    BOUND,
    NO
};

template<class T>
class finite_element_matrix final {
    std::array<metamath::linear::sparse_matrix<T>, 2> _part;

public:
    metamath::linear::sparse_matrix<T>& inner() noexcept { return _part[size_t(matrix_part::INNER)]; }
    metamath::linear::sparse_matrix<T>& bound() noexcept { return _part[size_t(matrix_part::BOUND)]; }
    const metamath::linear::sparse_matrix<T>& inner() const noexcept { return _part[size_t(matrix_part::INNER)]; }
    const metamath::linear::sparse_matrix<T>& bound() const noexcept { return _part[size_t(matrix_part::BOUND)]; }

    void clear() { _part = {}; }

    metamath::linear::sparse_matrix<T>& operator[](const matrix_part part) { return _part[size_t(part)]; }
    const metamath::linear::sparse_matrix<T>& operator[](const matrix_part part) const { return _part[size_t(part)]; }
};

}