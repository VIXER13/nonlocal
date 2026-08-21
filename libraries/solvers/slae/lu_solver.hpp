#pragma once

#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/linear/sparse_matrix.hpp>
#include <metamath/types/traits.hpp>
#include <metamath/utils/operators.hpp>

#include <stdexcept>
#include <vector>

namespace nonlocal::slae {

// Full LU direct solver for general sparse matrices.
// T may be scalar or square_matrix<floating_point, N> for block structure.
// The sparse input is expanded to a dense working buffer so all fill-in is handled automatically.
// Lower triangle stores L (unit diagonal implicit), upper triangle + diagonal stores U (inv(U_ii) at diagonal).
template<class T, std::integral I = uint32_t, std::integral J = size_t>
class lu_solver final {
    size_t _n;
    std::vector<std::vector<T>> _dense;

    void compute(const metamath::linear::sparse_matrix<T, I, J>& src) {
        using namespace metamath::linear;
        _dense.assign(_n, std::vector<T>(_n, T{}));
        for (const size_t i : std::ranges::iota_view{0zu, _n})
            for (const size_t s : src.portrait.shifts_range(i))
                _dense[i][src.portrait.indices[s]] = src.values[s];

        for (const size_t i : std::ranges::iota_view{0zu, _n}) {
            for (const size_t k : std::ranges::iota_view{0zu, i}) {
                for (const size_t p : std::ranges::iota_view{0zu, k})
                    _dense[i][k] -= _dense[i][p] * _dense[p][k];
                _dense[i][k] *= _dense[k][k]; // right-multiply by stored inv(U_kk)
            }
            for (const size_t j : std::ranges::iota_view{i, _n})
                for (const size_t p : std::ranges::iota_view{0zu, i})
                    _dense[i][j] -= _dense[i][p] * _dense[p][j];
            _dense[i][i] = inverse(_dense[i][i]);
        }
    }

public:
    using entity_t = metamath::types::container_type_t<T>;

    explicit lu_solver(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _n{matrix.cols()} {
        if (matrix.rows() != matrix.cols())
            throw std::invalid_argument{"LU solver requires a square matrix."};
        compute(matrix);
    }

    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const {
        using namespace metamath::linear;
        using metamath::operators::operator-=;
        if (rhs.size() != _n)
            throw std::invalid_argument{"LU solver requires rhs vector of the same size as the matrix."};
        std::vector<entity_t> result = rhs;

        for (const size_t i : std::ranges::iota_view{0zu, _n})
            for (const size_t k : std::ranges::iota_view{0zu, i})
                result[i] -= _dense[i][k] * result[k];

        for (const size_t i : std::ranges::iota_view{0zu, _n} | std::views::reverse) {
            for (const size_t j : std::ranges::iota_view{i + 1, _n})
                result[i] -= _dense[i][j] * result[j];
            result[i] = _dense[i][i] * result[i];
        }

        return result;
    }
};

}
