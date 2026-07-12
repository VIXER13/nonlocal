#pragma once

#include <cstdint>
#include <cstddef>
#include <concepts>

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
};

}