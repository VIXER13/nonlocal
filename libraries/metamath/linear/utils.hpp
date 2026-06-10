#pragma once

#include "sparse_matrix.hpp"

namespace metamath::linear {

template<std::integral I>
void validate_shifts(const std::vector<I>& shifts) {
    if (shifts.empty())
        return;
    if (shifts.size() < 2)
        throw std::logic_error{"Shifts vector shall have at least two elements."};
    if (shifts.front() != 0)
        throw std::logic_error{"Shifts vector shall start with 0."};
    for(const size_t i : std::ranges::iota_view{1zu, shifts.size()}) {
        if constexpr (std::is_signed_v<I>)
            if (shifts[i] < 0)
                throw std::logic_error{"Shifts vector cannot contain negative values."};
        if (shifts[i - 1] > shifts[i])
            throw std::logic_error{"Shifts vector shall be non-decreasing."};
    }
}

template<std::integral I, std::integral J>
void validate_sparse_matrix_portrait(const sparse_matrix_portrait<I, J>& portrait) {
    validate_shifts(portrait.shifts());
    if (portrait.indices().size() != portrait.non_zeros())
        throw std::logic_error{"The last element of shifts shall be equal to the size of indices."};
    for(const size_t row : std::ranges::iota_view{0zu, portrait.rows()}) {
        const auto shifts_range = portrait.shifts(row);
        if (shifts_range.size() > portrait.cols())
            throw std::logic_error{"Number of non-zero elements in a row cannot be greater than the number of columns."};
        for(const size_t shift : shifts_range) {
            if constexpr (std::is_signed_v<I>)
                if (portrait.indices()[shift] < 0)
                    throw std::logic_error{"Column index in indices vector cannot be negative."};
            if (portrait.indices()[shift] >= portrait.cols())
                throw std::logic_error{"Column index in indices vector is out of range."};
            if (shift < shifts_range.back()) {
                if (portrait.indices()[shift] == portrait.indices()[shift + 1])
                    throw std::logic_error{"Column indices in each row shall be unique."};
                if (portrait.indices()[shift + 1] < portrait.indices()[shift])
                    throw std::logic_error{"Column indices in each row shall be sorted."};
            }
        }
    }
}

template<class T, std::integral I, std::integral J>
void validate_sparse_matrix(const sparse_matrix<T, I, J>& matrix) {
    validate_sparse_matrix_portrait(matrix.portrait());
    if (matrix.values().size() != matrix.portrait().non_zeros())
        throw std::logic_error{"Values vector size shall be equal to the number of non-zero elements in the portrait."};
}

}