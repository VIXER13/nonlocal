#pragma once

#include "boundary_conditions_1d.hpp"

#include <ranges>
#include <unordered_map>

namespace nonlocal::solver_1d {

template<class T, class Vector>
void boundary_condition_first_kind_1d(Vector& f, const std::unordered_map<size_t, T>& matrix_bound,
                                      const boundary_data_1d<T>& condition, const size_t index) {
    if (condition.prescribed_value) {
        for (const auto& [i, val] : matrix_bound)
            f[i] -= val * *condition.prescribed_value;
        f[index] = *condition.prescribed_value;
    }
}

template<class T, class Vector>
void boundary_condition_first_kind_1d(Vector& f, const std::array<std::unordered_map<size_t, T>, 2>& matrix_bound,
                                      const boundaries_data_1d<T>& conditions) {
    const std::array<size_t, 2> indices = { 0, size_t(f.size() - 1) };
    for (const size_t b : std::ranges::iota_view{ 0zu, 2zu })
        boundary_condition_first_kind_1d(f, matrix_bound[b], conditions[b], indices[b]);
}

template<class T, physics_t Physics, class Vector>
void boundary_condition_first_kind_1d(Vector& f, const std::unordered_map<size_t, T>& matrix_bound,
                                      const boundary_condition_1d<T, Physics>& condition, const size_t index) {
    boundary_condition_first_kind_1d(f, matrix_bound, condition.data(), index);
}

template<class T, physics_t Physics, class Vector>
void boundary_condition_first_kind_1d(Vector& f, const std::array<std::unordered_map<size_t, T>, 2>& matrix_bound,
                                      const boundaries_conditions_1d<T, Physics>& conditions) {
    boundary_condition_first_kind_1d(f, matrix_bound, normalize_boundaries_1d(conditions));
}

} // namespace nonlocal::solver_1d
