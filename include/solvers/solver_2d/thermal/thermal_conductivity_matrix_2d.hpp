#ifndef THERMAL_CONDUCTIVITY_MATRIX_2D_HPP
#define THERMAL_CONDUCTIVITY_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"
#include "heat_equation_parameters_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class thermal_conductivity_matrix_2d : public finite_element_matrix_2d<1, T, I, Matrix_Index> {
    using _base = finite_element_matrix_2d<1, T, I, Matrix_Index>;

protected:
    T integrate_basic(const size_t e, const size_t i) const;
    template<material_t Material>
    T integrate_loc([[maybe_unused]] const decltype(equation_parameters_2d<T, Material>::lambda)& lambda,
                    const size_t e, const size_t i, const size_t j) const;
    template<material_t Material, class Influence_Function>
    T integrate_nonloc([[maybe_unused]] const decltype(equation_parameters_2d<T, Material>::lambda)& lambda,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(const std::vector<bool>& is_inner, const theory_t theory, const bool is_neumann);

    void neumann_problem_col_fill();

public:
    explicit thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);
    ~thermal_conductivity_matrix_2d() noexcept override = default;

    template<material_t Material, class Influence_Function>
    void calc_matrix(const decltype(equation_parameters_2d<T, Material>::lambda)& lambda,
                     const std::vector<bool>& is_inner,
                     const T p1, const Influence_Function& influence_fun,
                     const bool is_neumann = false);
};

template<class T, class I, class Matrix_Index>
thermal_conductivity_matrix_2d<T, I, Matrix_Index>::thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy)
    : _base{mesh_proxy} {}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * _base::mesh_proxy()->jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<material_t Material>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc([[maybe_unused]] const decltype(equation_parameters_2d<T, Material>::lambda)& lambda,
                                                                    const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    if constexpr (Material == material_t::ISOTROPIC) {
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj)
            integral += el->weight(q) * ((*dNdi)[X] * (*dNdj)[X] + (*dNdi)[Y] * (*dNdj)[Y]) / _base::mesh_proxy()->jacobian(*J);
    } else if constexpr (Material == material_t::ORTHOTROPIC) {
        std::array<T, 2> integral_part = {};
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj) {
            const T factor = el->weight(q) / _base::mesh_proxy()->jacobian(*J);
            integral_part[X] += factor * (*dNdi)[X] * (*dNdj)[X];
            integral_part[Y] += factor * (*dNdi)[Y] * (*dNdj)[Y];
        }
        integral = lambda[X] * integral_part[X] + lambda[Y] * integral_part[Y];
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<material_t Material, class Influence_Function>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc([[maybe_unused]] const decltype(equation_parameters_2d<T, Material>::lambda)& lambda,
                                                                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                                                       const Influence_Function& influence_function) const {
    T integral = 0;
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
    const auto  qcoordNL_begin = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_begin    = _base::mesh_proxy()->dNdX(eNL, jNL);
    std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1> integral_part = {};
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto dNdNL    = dNdNL_begin;
        auto qcoordNL = qcoordNL_begin;
        std::array<T, 2> inner_integral_part = {};
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
            const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
            inner_integral_part[X] += influence_weight * (*dNdNL)[X];
            inner_integral_part[Y] += influence_weight * (*dNdNL)[Y];
        }
        if constexpr (Material == material_t::ISOTROPIC)
            integral += elL->weight(qL) * (inner_integral_part[X] * (*dNdL)[X] + inner_integral_part[Y] * (*dNdL)[Y]);
        else if constexpr (Material == material_t::ORTHOTROPIC) {
            integral_part[X] += elL->weight(qL) * inner_integral_part[X] * (*dNdL)[X];
            integral_part[Y] += elL->weight(qL) * inner_integral_part[Y] * (*dNdL)[Y];
        }
    }
    if constexpr (Material == material_t::ORTHOTROPIC)
        integral = lambda[X] * integral_part[X] + lambda[Y] * integral_part[Y];
    return integral;
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::create_matrix_portrait(const std::vector<bool>& is_inner,
                                                                                const theory_t theory, const bool is_neumann) {
    if (is_neumann)
        for(const size_t row : std::views::iota(size_t{0}, size_t(_base::matrix_inner().rows())))
            _base::matrix_inner().outerIndexPtr()[row+1] = 1;
    if (theory == theory_t::LOCAL)
        _base::template create_matrix_portrait<theory_t::LOCAL>(is_inner);
    else if (theory == theory_t::NONLOCAL)
        _base::template create_matrix_portrait<theory_t::NONLOCAL>(is_inner);
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::neumann_problem_col_fill() {
#pragma omp parallel for default(none)
    for(size_t node = _base::mesh_proxy()->first_node(); node < _base::mesh_proxy()->last_node(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node - _base::mesh_proxy()->first_node(), _base::mesh().nodes_count());
        for(const I e : _base::mesh_proxy()->nodes_elements_map(node))
            val += integrate_basic(e, _base::mesh_proxy()->global_to_local_numbering(e).find(node)->second);
    }
}

template<class T, class I, class Matrix_Index>
template<material_t Material, class Influence_Function>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::calc_matrix(const decltype(equation_parameters_2d<T, Material>::lambda)& lambda,
                                                                     const std::vector<bool>& is_inner,
                                                                     const T p1, const Influence_Function& influence_fun,
                                                                     const bool is_neumann) {
    const theory_t theory = p1 < MAX_NONLOCAL_WEIGHT<T> ? theory_t::NONLOCAL : theory_t::LOCAL;
    const size_t rows = _base::mesh_proxy()->last_node() - _base::mesh_proxy()->first_node() +
                        (is_neumann && MPI_utils::MPI_rank() == MPI_utils::MPI_size() - 1),
                 cols = _base::mesh().nodes_count() + is_neumann;
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    create_matrix_portrait(is_inner, theory, is_neumann);

    T factor_loc = p1, factor_nonloc = T{1} - p1;
    if constexpr (Material == material_t::ISOTROPIC) {
        factor_loc *= lambda;
        factor_nonloc *= lambda;
    }
    _base::template calc_matrix(is_inner, theory, influence_fun,
        [this, factor_loc, &lambda = lambda](const size_t e, const size_t i, const size_t j) {
            return factor_loc * integrate_loc<Material>(lambda, e, i, j);
        },
        [this, factor_nonloc, &lambda = lambda]
        (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            return factor_nonloc * integrate_nonloc<Material>(lambda, eL, eNL, iL, jNL, influence_function);
        });
    if (is_neumann)
        neumann_problem_col_fill();
}

}

#endif