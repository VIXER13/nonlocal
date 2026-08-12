#pragma once

#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include <ranges>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <limits>
#include <array>

// Suppress re-instantiation of heavy templates for T=double in all test TUs.
// Explicit instantiation definitions live in metamath_instantiations/metamath_instantiations.cpp.
namespace metamath::finite_element {
    extern template class element_1d_integrate<double>;
    extern template std::unique_ptr<element_1d_base<double>>     make_element_1d<double>(size_t);
    extern template element_1d_integrate<double>                 make_element_1d_integrated<double>(size_t, size_t);
    extern template std::unique_ptr<quadrature_1d_base<double>>  make_quadrature_1d<double>(size_t);

    extern template class element_2d_integrate<double>;
    extern template std::unique_ptr<element_2d_base<double>>     make_triangle_element_2d<double>(size_t);
    extern template std::unique_ptr<element_2d_base<double>>     make_serendipity_element_2d<double>(size_t);
    extern template std::unique_ptr<element_2d_base<double>>     make_lagrangian_element_2d<double>(size_t);
    extern template std::unique_ptr<quadrature_2d_base<double>>  make_barycentric_quadrature_2d<double>(size_t);
    extern template std::unique_ptr<quadrature_2d_base<double>>  make_quadrature_2d<double>(geometry_t, size_t);
    extern template std::unique_ptr<quadrature_2d_base<double>>  make_quadrature_2d<double>(size_t, size_t);
    extern template element_2d_integrate<double>                 make_element_2d_integrated<double>(element_t, size_t, size_t);
}
