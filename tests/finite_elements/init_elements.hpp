#pragma once

#include "metamath.hpp"

namespace unit_tests {

// Element initialization is separated into a separate translation unit to avoid the SIOF problem.

template<class T>
std::vector<std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>> init_lagrangian_elements_1d();

template<class T>
std::vector<std::unique_ptr<metamath::finite_element::element_2d_integrate_base<T>>> init_triangle_elements_2d();

template<class T>
std::vector<std::unique_ptr<metamath::finite_element::element_2d_integrate_base<T>>> init_serendipity_elements_2d();

template<class T>
std::vector<std::unique_ptr<metamath::finite_element::element_2d_integrate_base<T>>> init_lagrangian_elements_2d();

}