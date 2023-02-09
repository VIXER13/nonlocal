#ifndef NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_2D_HPP
#define NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_2D_HPP

#include "finite_element_matrix_2d.hpp"
#include "thermal_parameters_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class thermal_conductivity_matrix_2d : public finite_element_matrix_2d<1, T, I, Matrix_Index> {
    using _base = finite_element_matrix_2d<1, T, I, Matrix_Index>;

protected:
    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_loc(const parameter_2d<T>& parameter, const size_t e, const size_t i, const size_t j) const;
    template<class Influence_Function>
    T integrate_nonloc(const parameter_2d<T>& parameter, 
                       const size_t eL, const size_t eNL, 
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence) const;

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const std::vector<bool>& is_inner, const bool is_neumann);

    void neumann_problem_col_fill();

public:
    explicit thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~thermal_conductivity_matrix_2d() noexcept override = default;

    void compute(const parameters_2d<T>& parameters, const std::vector<bool>& is_inner, const bool is_neumann = false);
};

template<class T, class I, class Matrix_Index>
thermal_conductivity_matrix_2d<T, I, Matrix_Index>::thermal_conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : std::ranges::iota_view{0u, el.nodes_count()})
        integral += el.weight(q) * el.qN(i, q) * mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
    return integral;
}

template<class T, class I, class Matrix_Index>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_loc(
    const parameter_2d<T>& parameter, const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& conductivity = parameter.conductivity;
    const auto& el = _base::mesh().container().element_2d(e);
    if (parameter.material == material_t::ISOTROPIC) {
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
            const std::array<T, 2>& dNi = _base::mesh().derivatives(e, i, q);
            const std::array<T, 2>& dNj = _base::mesh().derivatives(e, j, q);
            integral += el.weight(q) * (dNi[X] * dNj[X] + dNi[Y] * dNj[Y]) /
                        mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
        }
        integral *= conductivity[X][X];
    } else if (parameter.material == material_t::ORTHOTROPIC) {
        std::array<T, 2> integral_part = {};
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
            const std::array<T, 2>& dNi = _base::mesh().derivatives(e, i, q);
            const std::array<T, 2>& dNj = _base::mesh().derivatives(e, j, q);
            const T factor = el.weight(q) / mesh::jacobian(_base::mesh().jacobi_matrix(e, q));
            integral_part[X] += factor * dNi[X] * dNj[X];
            integral_part[Y] += factor * dNi[Y] * dNj[Y];
        }
        integral = conductivity[X][X] * integral_part[X] + conductivity[Y][Y] * integral_part[Y];
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
T thermal_conductivity_matrix_2d<T, I, Matrix_Index>::integrate_nonloc(const parameter_2d<T>& parameter, 
                                                                       const size_t eL, const size_t eNL, 
                                                                       const size_t iL, const size_t jNL,
                                                                       const Influence_Function& influence) const {
    std::array<T, 2> integral = {};
    const auto& conductivity = parameter.conductivity;
    const auto& elL  = _base::mesh().container().element_2d(eL );
    const auto& elNL = _base::mesh().container().element_2d(eNL);
    for(const size_t qL : std::ranges::iota_view{0u, elL.qnodes_count()}) {
        std::array<T, 2> inner_integral = {};
        for(const size_t qNL : std::ranges::iota_view{0u, elNL.qnodes_count()}) {
            const std::array<T, 2>& dNjNL = _base::mesh().derivatives(eNL, jNL, qNL);
            const T influence_weight = elNL.weight(qNL) * influence(_base::mesh().quad_coord(eL,  qL ), 
                                                                    _base::mesh().quad_coord(eNL, qNL));
            inner_integral[X] += influence_weight * dNjNL[X];
            inner_integral[Y] += influence_weight * dNjNL[Y];
        }
        const std::array<T, 2>& dNiL = _base::mesh().derivatives(eL, iL, qL);
        if (parameter.material == material_t::ISOTROPIC)
            integral[X] += elL.weight(qL) * (inner_integral[X] * dNiL[X] + inner_integral[Y] * dNiL[Y]);
        else if (parameter.material == material_t::ORTHOTROPIC) {
            integral[X] += elL.weight(qL) * inner_integral[X] * dNiL[X];
            integral[Y] += elL.weight(qL) * inner_integral[Y] * dNiL[Y];
        }
    }
    if (parameter.material == material_t::ORTHOTROPIC)
        return conductivity[X][X] * integral[X] + conductivity[Y][Y] * integral[Y];
    return conductivity[X][X] * integral[X];
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                                                                const std::vector<bool>& is_inner, const bool is_neumann) {
    if (is_neumann)
        for(const size_t row : std::views::iota(0u, size_t(_base::matrix_inner().rows())))
            _base::matrix_inner().outerIndexPtr()[row + 1] = 1;
    _base::init_shifts(theories, is_inner);
    _base::init_indices(theories, is_inner, false);
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, size_t(_base::matrix_inner().rows())})
            _base::matrix_inner().innerIndexPtr()[_base::matrix_inner().outerIndexPtr()[row + 1] - 1] = _base::mesh().container().nodes_count();
    utils::sort_indices(_base::matrix_inner());
    utils::sort_indices(_base::matrix_bound());
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::neumann_problem_col_fill() {
    const auto process_nodes = _base::mesh().process_nodes();
#pragma omp parallel for default(none) shared(process_nodes)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node - process_nodes.front(), _base::mesh().container().nodes_count());
        for(const I e : _base::mesh().elements(node))
            val += integrate_basic(e, _base::mesh().global_to_local(e, node));
    }
}

template<class T, class I, class Matrix_Index>
void thermal_conductivity_matrix_2d<T, I, Matrix_Index>::compute(const parameters_2d<T>& parameters,
                                                                 const std::vector<bool>& is_inner,
                                                                 const bool is_neumann) {
    const size_t rows = _base::mesh().process_nodes().size() + (is_neumann && parallel_utils::MPI_rank() == parallel_utils::MPI_size() - 1);
    const size_t cols = _base::mesh().container().nodes_count() + is_neumann;
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    create_matrix_portrait(theories, is_inner, is_neumann);
    _base::calc_coeffs(theories, is_inner,
        [this, &parameters](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const auto& parameter = parameters.at(group);
            //std::cout << "e = " << e << " cond = " << parameter.physical.conductivity[X][X] << std::endl;
            return parameter.model.local_weight * integrate_loc(parameter.physical, e, i, j);
        },
        [this, &parameters](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            const auto& parameter = parameters.at(group);
            return nonlocal_weight(parameter.model.local_weight) * 
                   integrate_nonloc(parameter.physical, eL, eNL, iL, jNL, parameter.model.influence);
        });
    if (is_neumann)
        neumann_problem_col_fill();
}

}

#endif