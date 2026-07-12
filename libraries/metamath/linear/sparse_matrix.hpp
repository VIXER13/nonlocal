#pragma once

#include "self_adjoint_view.hpp"
#include "sparse_matrix_portrait.hpp"

namespace metamath::linear {

template<class T, std::integral I, std::integral J>
struct sparse_matrix final {
    sparse_matrix_portrait<I, J> portrait;
    std::vector<T> values;

    sparse_matrix() = default;
    explicit sparse_matrix(const size_t rows, const size_t cols)
        : portrait{rows, cols} {}

    template<matrix_part Part>
    self_adjoint_view<Part, T, I, J> self_adjoint() const noexcept {
        return {*this};
    }

    size_t rows() const {
        return portrait.rows();
    }

    size_t cols() const {
        return portrait.cols();
    }

    size_t non_zeros() const {
        return portrait.non_zeros();
    }

    void allocate_values() {
        values.resize(portrait.non_zeros(), T{});
    }

    T& operator()(const size_t row, const size_t col) {
        return values[portrait.shift(row, col)];
    }

    const T& operator()(const size_t row, const size_t col) const {
        return values[portrait.shift(row, col)];
    }
};

}