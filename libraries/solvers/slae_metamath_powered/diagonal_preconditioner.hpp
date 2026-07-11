#pragma once

#include "preconditioner_base.hpp"

#include <metamath/linear/sparse_matrix.hpp>
#include <metamath/linear/fixed_matrix.hpp>

namespace nonlocal::slae {

template<class T>
class diagonal_preconditioner final : public preconditioner_base<T> {
    std::vector<T> _inverse_diagonal;

public:
    using typename preconditioner_base<T>::entity_t;

    template<std::integral I, std::integral J>
    explicit diagonal_preconditioner(const metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _inverse_diagonal(matrix.cols()) {
        if (matrix.rows() != matrix.cols())
            throw std::invalid_argument{"Diagonal preconditioner requires a square matrix."};
        for(const size_t i : std::ranges::iota_view{0zu, _inverse_diagonal.size()})
            _inverse_diagonal[i] = metamath::linear::inverse(matrix(i, i));
    }

    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        using namespace metamath::linear;
        if (rhs.size() != _inverse_diagonal.size())
            throw std::invalid_argument{"Diagonal preconditioner requires rhs vector of the same size as the matrix."};
        std::vector<entity_t> result(rhs.size());
        for (const size_t i : std::ranges::iota_view{0zu, rhs.size()})
            result[i] = _inverse_diagonal[i] * rhs[i];
        return result;
    }
};

}
