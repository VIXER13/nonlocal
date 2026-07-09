#pragma once

#include "preconditioner_base.hpp"

#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/types/traits.hpp>

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
class diagonal_preconditioner final : public preconditioner_base<T, I, J> {
    std::vector<T> _inverse_diagonal;

public:
    using typename preconditioner_base<T, I, J>::entity_t;

    void compute(metamath::linear::sparse_matrix<T, I, J>&& matrix) override {
        _inverse_diagonal.resize(matrix.cols());
        for(const size_t i : std::ranges::iota_view{0zu, _inverse_diagonal.size()})
            _inverse_diagonal[i] = metamath::linear::inverse(matrix(i, i));
    }

    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        using namespace metamath::linear;
        std::vector<entity_t> result(rhs.size());
        for (const size_t i : std::ranges::iota_view{0zu, rhs.size()})
            result[i] = _inverse_diagonal[i] * rhs[i];
        return result;
    }
};

}
