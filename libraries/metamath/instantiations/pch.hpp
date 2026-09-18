#pragma once

#include <metamath/metamath.hpp>

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

// Suppress re-instantiation of heavy templates for T=double in all TUs.
// Explicit instantiation definitions live in libraries/metamath/instantiations/metamath_instantiations.cpp.
namespace metamath::linear {
    extern template struct sparse_matrix_portrait<uint32_t, size_t>;

    extern template struct sparse_matrix<double, uint32_t, size_t>;
    extern template sparse_matrix<double, uint32_t, size_t>& operator*=(sparse_matrix<double, uint32_t, size_t>&, double);
    extern template sparse_matrix<double, uint32_t, size_t>& operator/=(sparse_matrix<double, uint32_t, size_t>&, double);
    extern template std::vector<double> operator*(const sparse_matrix<double, uint32_t, size_t>&, const std::vector<double>&);
    extern template sparse_matrix<double, uint32_t, size_t>& operator+=(sparse_matrix<double, uint32_t, size_t>&, const sparse_matrix<double, uint32_t, size_t>&);

    extern template struct sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>;
    extern template sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>& operator*=(sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, double);
    extern template sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>& operator/=(sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, double);
    extern template std::vector<std::array<double, 2>> operator*(const sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, const std::vector<std::array<double, 2>>&);
    extern template sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>& operator+=(sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, const sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&);
}

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
