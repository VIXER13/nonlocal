#pragma once

#include "self_adjoint_view.hpp"
#include "sparse_matrix_portrait.hpp"

namespace metamath::linear {

template<class T, std::integral I, std::integral J>
class sparse_matrix final {
    sparse_matrix_portrait<I, J> _portrait;
    std::vector<T> _values;

public:
    sparse_matrix() = default;
    explicit sparse_matrix(const size_t rows, const size_t cols)
        : _portrait{rows, cols} {}

    sparse_matrix_portrait<I, J>& portrait() noexcept { return _portrait; }
    std::vector<T>& values() noexcept { return _values; }
    const sparse_matrix_portrait<I, J>& portrait() const noexcept { return _portrait; }
    const std::vector<T>& values() const noexcept { return _values; }

    template<matrix_part Part>
    self_adjoint_view<Part, T, I, J> self_adjoint() const noexcept {
        return {*this};
    }

    size_t rows() const {
        return portrait().rows();
    }

    size_t cols() const {
        return portrait().cols();
    }

    size_t non_zeros() const {
        return portrait().non_zeros();
    }

    void allocate_values() {
        values().resize(portrait().non_zeros(), T{});
    }

    T& operator()(const size_t row, const size_t col) {
        return values()[portrait().shift(row, col)];
    }

    const T& operator()(const size_t row, const size_t col) const {
        return values()[portrait().shift(row, col)];
    }
};

}