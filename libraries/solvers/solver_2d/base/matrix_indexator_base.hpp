#pragma once

#include "finite_element_matrix.hpp"

#include "indexator_base.hpp"

#include <array>
#include <vector>

namespace nonlocal {

template<size_t DoF>
class matrix_indexator_base : public mesh::indexator_base {
    std::array<std::array<std::vector<bool>, DoF>, 2> _flags;

    std::vector<bool>& flags(const matrix_part part, const size_t dof);

protected:
    explicit matrix_indexator_base(const size_t size, const bool is_symmetric);
    ~matrix_indexator_base() noexcept override = default;

    template<class Callback>
    void check_flag(const matrix_part part, size_t dof, const size_t col, const Callback& callback);

public:
    void reset(const size_t node) override;
};

template<size_t DoF>
std::vector<bool>& matrix_indexator_base<DoF>::flags(const matrix_part part, const size_t dof) {
    return _flags[size_t(part)][dof];
}

template<size_t DoF>
template<class Callback>
void matrix_indexator_base<DoF>::check_flag(const matrix_part part, size_t dof, const size_t col, const Callback& callback) {
    if (part == matrix_part::NO)
        return;
    if (auto& flag = flags(part, dof); !flag[col]) {
        callback();
        flag[col] = true;
    }
}

template<size_t DoF>
matrix_indexator_base<DoF>::matrix_indexator_base(const size_t size, const bool is_symmetric)
    : indexator_base{is_symmetric} {
    for(auto& part : _flags)
        for(auto& dof : part)
            dof.resize(size);
}

template<size_t DoF>
void matrix_indexator_base<DoF>::reset(const size_t node) {
    for(const size_t dof : std::ranges::iota_view{0u, DoF}) {
        auto& bound = flags(matrix_part::BOUND, dof);
        auto& inner = flags(matrix_part::INNER, dof);
        std::fill(bound.begin(), bound.end(), false);
        std::fill(std::next(inner.begin(), is_symmetric() ? DoF * node : 0), inner.end(), false);
    }
}

}