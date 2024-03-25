#ifndef TESTS_PROBLEMS_UTILS_1D_HPP
#define TESTS_PROBLEMS_UTILS_1D_HPP

#include "make_element_1d.hpp"
#include "nonlocal_config.hpp"
#include "mesh_1d.hpp"


namespace nonstat_1d_tests {

using namespace nonlocal;
using namespace nonlocal::thermal;

template<std::floating_point T>
std::shared_ptr<mesh::mesh_1d<T>> make_mesh_1d(
    const std::vector<mesh::segment_data<T>>& segments,
    const config::order_t element_order, const config::order_t& quadrature_order) {
    return std::make_shared<mesh::mesh_1d<T>>(
        make_element<T>(element_order, quadrature_order),
        segments
    );
}

}

#endif