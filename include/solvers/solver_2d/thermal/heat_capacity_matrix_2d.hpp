#ifndef NONLOCAL_HEAT_CAPACITY_MATRIX_2D_HPP
#define NONLOCAL_HEAT_CAPACITY_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"
#include "heat_equation_parameters.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class heat_capacity_matrix_2d : public finite_element_matrix_2d<1, T, I, Matrix_Index> {
    using _base = finite_element_matrix_2d<1, T, I, Matrix_Index>;

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

public:
    explicit heat_capacity_matrix_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);
    ~heat_capacity_matrix_2d() noexcept override = default;

    template<material_t Material>
    void calc_matrix(const equation_parameters<T, Material>& eq_parameters,
                     const std::vector<bool>& is_inner);
};

template<class T, class I, class Matrix_Index>
heat_capacity_matrix_2d<T, I, Matrix_Index>::heat_capacity_matrix_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy)
    : _base{mesh_proxy} {}

template<class T, class I, class Matrix_Index>
T heat_capacity_matrix_2d<T, I, Matrix_Index>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q) * _base::jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<material_t Material>
void heat_capacity_matrix_2d<T, I, Matrix_Index>::calc_matrix(const equation_parameters<T, Material>& eq_parameters,
                                                              const std::vector<bool>& is_inner) {
    const size_t rows = _base::mesh_proxy()->last_node() - _base::mesh_proxy()->first_node(),
                 cols = _base::mesh().nodes_count();
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    _base::create_matrix_portrait<theory_t::LOCAL>(is_inner);
    static constexpr auto NO_INFLUENCE = []() constexpr noexcept {};
    _base::calc_matrix(is_inner, theory_t::LOCAL, NO_INFLUENCE,
        [this, factor = eq_parameters.c * eq_parameters.rho](const size_t e, const size_t i, const size_t j) {
            return factor * integrate_basic_pair(e, i, j);
        },
        [](const size_t, const size_t, const size_t, const size_t, const decltype(NO_INFLUENCE)) constexpr noexcept { return 0; }
    );
}

}

#endif