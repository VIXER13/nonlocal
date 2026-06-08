#pragma once

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace metamath::linear {

template<class I = uint32_t, class J = size_t>
class sparce_matrix_portrait final {
    std::vector<J> _shifts;
    std::vector<I> _indices;

public:
    sparce_matrix_portrait() = default;
    sparce_matrix_portrait(const size_t rows)
        : _shifts(rows + 1zu, 0zu) {}

    std::vector<J>& shifts() noexcept { return _shifts; }
    std::vector<I>& indices() noexcept { return _indices; }
    const std::vector<J>& shifts() const noexcept { return _shifts; }
    const std::vector<I>& indices() const noexcept { return _indices; }

    size_t rows() const {
        return _shifts.empty() ? 0zu : _shifts.size() - 1;
    }

    size_t cols() const {
        return _indices.empty() ? 0zu : *std::ranges::max_element(_indices) + 1;
    }

    size_t non_zeros() const {
        return _indices.empty() ? 0zu : _indices.back();
    }

    bool contains(const size_t row, const size_t col) const {
        if (row >= rows())
            return false;
        return std::binary_search(&_indices[_shifts[row]], &_indices[_shifts[row + 1]], col);
    }

    size_t index(const size_t row, const size_t col) const {
        if (row >= rows())
            throw std::out_of_range{"Row index is out of range."};
        const auto it = std::lower_bound(&_indices[_shifts[row]], &_indices[_shifts[row + 1]], col);
        if (it == &_indices[_shifts[row + 1]] || *it != col)
            throw std::out_of_range{"Column index is out of range."};
        return std::distance(_indices.data(), it);
    }

    void set_rows_count(const size_t rows) {
        _shifts.resize(rows + 1zu, 0zu);
    }

    void accumulate_shifts() {
        if (!_shifts.empty())
            for(const size_t row : std::ranges::iota_view{0u, _shifts.size() - 1})
                _shifts[row + 1] += _shifts[row];
    }

    void allocate_indices() {
        indices().resize(non_zeros(), 0zu);
    }

    void sort_indices() {
        if (!_shifts.empty()) {
#pragma omp parallel for schedule(dynamic)
            for(const size_t row : std::ranges::iota_view{0u, _shifts.size() - 1})
                std::sort(&_indices[_shifts[row]], &_indices[_shifts[row + 1]]);
        }
    }
};

}