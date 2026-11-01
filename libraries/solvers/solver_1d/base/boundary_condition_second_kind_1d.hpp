#pragma once

#include "boundary_conditions_1d.hpp"

#include <ranges>

namespace nonlocal::solver_1d {

template<class T, class Vector>
void boundary_condition_second_kind_1d(Vector& f, const boundary_data_1d<T>& condition, const size_t index) {
    if (!condition.prescribed_value)
        f[index] += condition.load;
}

template<class T, class Vector>
void boundary_condition_second_kind_1d(Vector& f, const boundaries_data_1d<T>& conditions,
                                       const bool is_neumann = false) {
    const std::array<size_t, 2> indices = { 0, size_t(f.size() - 1 - is_neumann) };
    for (const size_t b : std::ranges::iota_view{ 0zu, 2zu })
        boundary_condition_second_kind_1d(f, conditions[b], indices[b]);
}

template<class T, physics_t Physics, class Vector>
void boundary_condition_second_kind_1d(Vector& f, const boundary_condition_1d<T, Physics>& condition,
                                       const size_t index) {
    boundary_condition_second_kind_1d(f, condition.data(), index);
}

template<class T, physics_t Physics, class Vector>
void boundary_condition_second_kind_1d(Vector& f, const boundaries_conditions_1d<T, Physics>& conditions,
                                       const bool is_neumann = false) {
    boundary_condition_second_kind_1d(f, normalize_boundaries_1d(conditions), is_neumann);
}

} // namespace nonlocal::solver_1d
