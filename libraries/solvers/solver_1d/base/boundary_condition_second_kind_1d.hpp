#ifndef NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_1D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_1D_HPP

#include "boundary_conditions_1d.hpp"

#include <ranges>

namespace nonlocal {

template<class T, class Vector>
void boundary_condition_second_kind_1d(Vector& f, 
                                       const boundary_condition_1d<T>& boundary_condition,
                                       const size_t index) {
    if (const auto* const condition_ptr = dynamic_cast<const second_kind_1d<T>*>(&boundary_condition)) {
        const auto& condition = *condition_ptr;
        f[index] += condition();
    }
}

template<class T, class Vector, class Condition>
void boundary_condition_second_kind_1d(Vector& f,
                                       const std::array<Condition, 2>& boundaries_conditions,
                                       const bool is_neumann = false) {
    const std::array<size_t, 2> indices = {0, size_t(f.size() - 1 - is_neumann)};
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        boundary_condition_second_kind_1d(f, *boundaries_conditions[b], indices[b]);
}

}

#endif