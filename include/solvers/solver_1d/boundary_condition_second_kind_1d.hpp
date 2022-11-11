#ifndef NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_1D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_1D_HPP

#include "boundary_conditions_1d.hpp"

namespace nonlocal {

template<class T, class Vector>
void boundary_condition_second_kind_1d(Vector& f, 
                                       const stationary_boundary_condition_1d& boundary_condition,
                                       const size_t index) {
    if (const auto* const condition_ptr = dynamic_cast<const stationary_second_kind_1d<T>*>(&boundary_condition)) {
        const auto& condition = *condition_ptr;
        f[index] += condition();
    }
}

}

#endif