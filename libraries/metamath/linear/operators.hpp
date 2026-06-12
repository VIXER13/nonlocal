#pragma once

#include "sparse_matrix.hpp"

#include <metamath/utils/operators.hpp>

namespace metamath::linear {

template<class T, class U, std::integral I, std::integral J>
sparse_matrix<T, I, J>& operator*=(sparse_matrix<T, I, J>& matrix, const U scalar) {
    for(T& value : matrix.values)
         value *= scalar;
     return matrix;
}

template<class T, class U, std::integral I, std::integral J>
sparse_matrix<T, I, J>& operator/=(sparse_matrix<T, I, J>& matrix, const U scalar) {
    for(T& value : matrix.values)
        value /= scalar;
    return matrix;
}

template<class T, class U, std::integral I, std::integral J>
std::vector<U> operator*(const sparse_matrix<T, I, J>& matrix, const std::vector<U>& vector) {
    if (matrix.cols() != vector.size())
        throw std::invalid_argument{"Matrix columns count must match vector size for multiplication."};

    using namespace operators;
    std::vector<U> result(matrix.rows(), U{});
#pragma omp parallel for
    for(size_t row = 0; row < matrix.rows(); ++row)
        for(const size_t shift : matrix.portrait.shifts_range(row))
            result[row] += matrix.values[shift] * vector[matrix.portrait.indices[shift]];
    return result;
}

template<matrix_part Part, class T, class U, std::integral I, std::integral J>
std::vector<U> operator*(const self_adjoint_view<Part, T, I, J>& view, const std::vector<U>& vector) {
    if (view.matrix.cols() != vector.size())
        throw std::invalid_argument{"Matrix columns count must match vector size for multiplication."};

    using namespace operators;
    static constexpr std::conditional_t<Part == matrix_part::Upper, std::greater<>, std::less<>> comparator{};
    std::vector<U> result(view.matrix.rows(), U{});
    for(const size_t row : std::ranges::iota_view{0zu, view.matrix.rows()})
        for(const size_t shift : view.matrix.portrait.shifts_range(row))
            if (const size_t col = view.matrix.portrait.indices[shift]; row == col)
                result[row] += view.matrix.values[shift] * vector[col];
            else if (comparator(col, row)) {
                result[row] += view.matrix.values[shift] * vector[col];
                result[col] += view.matrix.values[shift] * vector[row];
            }
    return result;
}

}