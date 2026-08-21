#pragma once

#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/linear/sparse_matrix.hpp>
#include <metamath/types/traits.hpp>
#include <metamath/utils/operators.hpp>

#include <stdexcept>
#include <vector>

namespace nonlocal::slae {

// Full LDLT direct solver for symmetric sparse matrices stored in upper triangular format.
// T may be scalar or square_matrix<floating_point, N> for block structure.
// The sparse input is expanded to a dense working buffer so all fill-in is handled automatically.
// Upper triangle of the buffer stores U_ij = inv(D_i) * (A_ij - fill); diagonal stores inv(D_i).
template<class T, std::integral I = uint32_t, std::integral J = size_t>
class ldlt_solver final {
    size_t _n;
    std::vector<std::vector<T>> _dense;

    void compute(const metamath::linear::sparse_matrix<T, I, J>& src) {
        using namespace metamath::linear;
        _dense.assign(_n, std::vector<T>(_n, T{}));
        for (const size_t i : std::ranges::iota_view{0zu, _n})
            for (const size_t s : src.portrait.shifts_range(i))
                _dense[i][src.portrait.indices[s]] = src.values[s];

        std::vector<T> diagonal(_n);
        for (const size_t i : std::ranges::iota_view{0zu, _n}) {
            diagonal[i] = _dense[i][i];
            for (const size_t k : std::ranges::iota_view{0zu, i}) {
                const T& uki = _dense[k][i];
                diagonal[i] -= transpose(uki) * diagonal[k] * uki;
            }
            _dense[i][i] = inverse(diagonal[i]);

            for (const size_t j : std::ranges::iota_view{i + 1, _n}) {
                for (const size_t k : std::ranges::iota_view{0zu, i})
                    _dense[i][j] -= transpose(_dense[k][i]) * diagonal[k] * _dense[k][j];
                _dense[i][j] = _dense[i][i] * _dense[i][j];
            }
        }
    }

public:
    using entity_t = metamath::types::container_type_t<T>;

    explicit ldlt_solver(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _n{matrix.cols()} {
        if (matrix.rows() != matrix.cols())
            throw std::invalid_argument{"LDLT solver requires a square matrix."};
        compute(matrix);
    }

    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const {
        using namespace metamath::linear;
        using metamath::operators::operator-=;
        if (rhs.size() != _n)
            throw std::invalid_argument{"LDLT solver requires rhs vector of the same size as the matrix."};
        std::vector<entity_t> result = rhs;

        for (const size_t i : std::ranges::iota_view{0zu, _n})
            for (const size_t k : std::ranges::iota_view{0zu, i})
                result[i] -= transpose(_dense[k][i]) * result[k];

        for (const size_t i : std::ranges::iota_view{0zu, _n})
            result[i] = _dense[i][i] * result[i];

        for (const size_t i : std::ranges::iota_view{0zu, _n} | std::views::reverse)
            for (const size_t j : std::ranges::iota_view{i + 1, _n})
                result[i] -= _dense[i][j] * result[j];

        return result;
    }
};

}
