#pragma once

#include <cstdint>
#include <cstddef>
#include <concepts>
#include <stdexcept>
#include <string>

namespace metamath::linear {

enum class matrix_part : bool {
    Upper,
    Lower
};

template<class T, std::integral I = uint32_t, std::integral J = size_t>
class sparse_matrix;

template<matrix_part Part, class T, std::integral I, std::integral J>
struct self_adjoint_view final {
    const sparse_matrix<T, I, J>& matrix;

    explicit self_adjoint_view(const sparse_matrix<T, I, J>& matrix) : matrix{matrix} {
        if (matrix.rows() != matrix.cols())
            throw std::invalid_argument{"Matrix must be square for self-adjoint view. " + 
                                        std::to_string(matrix.rows()) + "x" + std::to_string(matrix.cols()) + " matrix is not square."};
    }
};

}