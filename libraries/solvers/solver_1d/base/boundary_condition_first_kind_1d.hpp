#ifndef NONLOCAL_BOUNDARY_CONDITION_FIRST_KIND_1D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_FIRST_KIND_1D_HPP

#include "boundary_conditions_1d.hpp"

#include <ranges>
#include <unordered_map>

namespace nonlocal {

template<class T, class Vector>
void boundary_condition_first_kind_1d(Vector& f,
                                      const std::unordered_map<size_t, T>& matrix_bound,
                                      const boundary_condition_1d<T>& boundary_condition,
                                      const size_t index) {
    if (const auto* const condition_ptr = dynamic_cast<const first_kind_1d<T>*>(&boundary_condition)) {
        const auto& condition = *condition_ptr;
        for(const auto& [i, val] : matrix_bound)
            f[i] -= val * condition();
        f[index] = condition();
    }
}

template<class T, class Vector, class Condition>
void boundary_condition_first_kind_1d(Vector& f,
                                      const std::array<std::unordered_map<size_t, T>, 2>& matrix_bound,
                                      const std::array<Condition, 2>& boundaries_conditions) {
    const std::array<size_t, 2> indices = {0, size_t(f.size() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        boundary_condition_first_kind_1d(f, matrix_bound[b], *boundaries_conditions[b], indices[b]);
}

}

#endif