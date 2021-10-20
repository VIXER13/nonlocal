#ifndef HEAT_CAPACITY_MATRIX_1D_HPP
#define HEAT_CAPACITY_MATRIX_1D_HPP

#include "finite_element_matrix_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class heat_capacity_matrix_1d : public finite_element_matrix_1d<T, I> {
    using _base = finite_element_matrix_base_1d<T, I>;

protected:

public:
    ~heat_capacity_matrix_1d() override = default;
};

}

#endif