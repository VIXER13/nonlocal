#include <metamath/metamath.hpp>

namespace metamath::linear {
    template struct sparse_matrix_portrait<uint32_t, size_t>;

    template struct sparse_matrix<double, uint32_t, size_t>;
    template sparse_matrix<double, uint32_t, size_t>& operator*=(sparse_matrix<double, uint32_t, size_t>&, double);
    template sparse_matrix<double, uint32_t, size_t>& operator/=(sparse_matrix<double, uint32_t, size_t>&, double);
    template std::vector<double> operator*(const sparse_matrix<double, uint32_t, size_t>&, const std::vector<double>&);
    template sparse_matrix<double, uint32_t, size_t>& operator+=(sparse_matrix<double, uint32_t, size_t>&, const sparse_matrix<double, uint32_t, size_t>&);

    template struct sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>;
    template sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>& operator*=(sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, double);
    template sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>& operator/=(sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, double);
    template std::vector<std::array<double, 2>> operator*(const sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, const std::vector<std::array<double, 2>>&);
    template sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>& operator+=(sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&, const sparse_matrix<square_matrix<double, 2>, uint32_t, size_t>&);
}

namespace metamath::finite_element {
    template class element_1d_integrate<double>;
    template std::unique_ptr<element_1d_base<double>>     make_element_1d<double>(size_t);
    template element_1d_integrate<double>                 make_element_1d_integrated<double>(size_t, size_t);
    template std::unique_ptr<quadrature_1d_base<double>>  make_quadrature_1d<double>(size_t);

    template class element_2d_integrate<double>;
    template std::unique_ptr<element_2d_base<double>>     make_triangle_element_2d<double>(size_t);
    template std::unique_ptr<element_2d_base<double>>     make_serendipity_element_2d<double>(size_t);
    template std::unique_ptr<element_2d_base<double>>     make_lagrangian_element_2d<double>(size_t);
    template std::unique_ptr<quadrature_2d_base<double>>  make_barycentric_quadrature_2d<double>(size_t);
    template std::unique_ptr<quadrature_2d_base<double>>  make_quadrature_2d<double>(geometry_t, size_t);
    template std::unique_ptr<quadrature_2d_base<double>>  make_quadrature_2d<double>(size_t, size_t);
    template element_2d_integrate<double> make_element_2d_integrated<double>(element_t, size_t, size_t);
}
