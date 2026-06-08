#pragma once

#include "sparce_matrix_portrait.hpp"

namespace metamath::linear {

template<class T, class I = uint32_t, class J = size_t>
class sparce_matrix final {
    sparce_matrix_portrait<I, J> _portrait;
    std::vector<T> _values;

public:
    sparce_matrix() = default;
    explicit sparce_matrix(const size_t rows)
        : _portrait{rows} {}

    sparce_matrix_portrait<I, J>& portrait() noexcept { return _portrait; }
    std::vector<T>& values() noexcept { return _values; }
    const sparce_matrix_portrait<I, J>& portrait() const noexcept { return _portrait; }
    const std::vector<T>& values() const noexcept { return _values; }

    void allocate_values() {
        values().resize(portrait().non_zeros(), T{});
    }

    T& operator()(const size_t row, const size_t col) {
        return values()[portrait().index(row, col)];
    }

    const T& operator()(const size_t row, const size_t col) const {
        return values()[portrait().index(row, col)];
    }
};

}