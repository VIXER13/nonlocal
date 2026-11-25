#pragma once

#include "self_adjoint_view.hpp"
#include "sparse_matrix_portrait.hpp"

namespace metamath::linear {

template<class T, std::integral I, std::integral J>
struct sparse_matrix final {
    sparse_matrix_portrait<I, J> portrait;
    std::vector<T> values;

    sparse_matrix() = default;
    explicit sparse_matrix(const size_t rows, const size_t cols) : portrait{ rows, cols } {}

    explicit sparse_matrix(sparse_matrix_portrait<I, J> portrait, const bool allocate = true) : portrait{ std::move(portrait) } {
        if (allocate)
            allocate_values();
    }

    template<matrix_part Part>
    self_adjoint_view<Part, T, I, J> self_adjoint() const {
        return self_adjoint_view<Part, T, I, J>{ *this };
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

    void clear() {
        portrait.shifts.clear();
        portrait.indices.clear();
        values.clear();
    }

    void set_zero() {
        std::fill(values.begin(), values.end(), T{});
    }

    template<std::ranges::input_range Rows>
    void set_zero(const Rows& rows) {
        if (values.size() != non_zeros())
            throw std::logic_error{
                "Cannot set zero values because the values vector size does not match the non-zero elements number."
            };

        for (const size_t row : rows) {
            portrait.check_row(row);
            std::fill(std::next(values.begin(), portrait.shifts[row]),
                      std::next(values.begin(), portrait.shifts[row + 1]), T{}); ///
        }
    }

    T& operator()(const size_t row, const size_t col) {
        return values[portrait.shift(row, col)];
    }

    const T& operator()(const size_t row, const size_t col) const {
        return values[portrait.shift(row, col)];
    }
};

} // namespace metamath::linear