#pragma once

#include <algorithm>
#include <cstdint>
#include <concepts>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace metamath::linear {

template<std::integral I = uint32_t, std::integral J = size_t>
struct sparse_matrix_portrait final {
    // Sparse matrix portrait in compressed sparse row (CSR) format.
    // Shifts vector has size (rows + 1) and contains the starting index of each row in the indices vector.
    // The last element of shifts is equal to the size of indices vector.
    // Indices vector contains the column indices of non-zero elements in the matrix, sorted within each row.
    // The number of non-zero elements in the matrix is equal to the last element of shifts vector.
    std::vector<J> shifts;
    std::vector<I> indices;
    size_t columns = 0zu;

    sparse_matrix_portrait() = default;
    sparse_matrix_portrait(const size_t rows, const size_t cols)
        : shifts(rows + 1zu, 0zu), columns{cols} {}

    size_t rows() const {
        return shifts.empty() ? 0zu : shifts.size() - 1;
    }

    size_t cols() const {
        return shifts.empty() ? 0zu : columns;
    }

    size_t non_zeros() const {
        return shifts.empty() ? 0zu : shifts.back();
    }

    bool contains(const size_t row, const size_t col) const {
        if (row >= rows())
            return false;
        return std::binary_search(&indices[shifts[row]], &indices[shifts[row + 1]], col);
    }

    void check_row(const size_t row) const {
        if (row >= rows())
            throw std::out_of_range{"Row index " + std::to_string(row) + " is out of range."};
    }

    std::ranges::iota_view<J, J> shifts_range(const size_t row) const {
        check_row(row);
        return std::ranges::iota_view<J, J>{shifts[row], shifts[row + 1]};
    }

    size_t shift(const size_t row, const size_t col) const {
        check_row(row);
        const auto range = shifts_range(row);
        const auto it = std::lower_bound(&indices[*range.begin()], &indices[*range.end()], col);
        if (it == &indices[*range.end()] || *it != col)
            throw std::out_of_range{"Column index " + std::to_string(col) + " is out of range on the row " + std::to_string(row) + "."};
        return std::distance(indices.data(), it);
    }

    void set_size(const size_t rows, const size_t cols) {
        shifts.resize(rows + 1zu, 0zu);
        columns = cols;
    }

    void accumulate_shifts() {
        if (!shifts.empty())
            for(const size_t row : std::ranges::iota_view{0u, rows()})
                shifts[row + 1] += shifts[row];
    }

    void allocate_indices() {
        indices.resize(non_zeros(), 0zu);
    }

    void sort_indices() {
        if (!shifts.empty()) {
#pragma omp parallel for schedule(dynamic)
            for(const size_t row : std::ranges::iota_view{0u, rows()})
                std::sort(&indices[shifts[row]], &indices[shifts[row + 1]]);
        }
    }
};

}