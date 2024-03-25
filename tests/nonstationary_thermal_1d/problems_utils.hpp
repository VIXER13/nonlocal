#ifndef TESTS_PROBLEMS_UTILS_1D_HPP
#define TESTS_PROBLEMS_UTILS_1D_HPP

#include "make_element_1d.hpp"
#include "mesh_1d.hpp"


namespace nonstat_1d_tests {

using namespace nonlocal;

template<class T, std::signed_integral I>
std::shared_ptr<mesh::mesh_1d<T>> make_mesh_1d(const std::vector<mesh::segment_data<T>>& segments, I element_order, I quadrature_order) {
    return std::make_shared<mesh::mesh_1d<T>>(
        make_element<T, I>(element_order, quadrature_order),
        segments
    );
}

}

#endif