#pragma once

#include "assebmler_base.hpp"

#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_1d {

template<std::floating_point T, std::integral I>
class spatial_assembler_1d : public assembler_base_1d<T, I> {
    using _base = assembler_base_1d<T, I>;

    std::vector<T> _solution;

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_local(const coefficient_t<T, 1>& coefficient, const size_t e, const size_t i, const size_t j) const;
    T integrate_nonlocal(const coefficient_t<T, 1>& coefficient, const std::function<T(const T&, const T&)>& influence,
                         const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void integral_condition(); // for Neumann problem

  public:
    explicit spatial_assembler_1d(finite_element_matrix_1d<T, I>& matrix, const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                  const std::optional<utils::nodes_sequence>& nodes_to_assemble = std::nullopt);
    ~spatial_assembler_1d() override = default;

    template<class Parameters, class Get_Coefficient>
    void calc_matrix(const Parameters& parameters, const problem_settings& settings, const Get_Coefficient get_coefficient,
                     const std::optional<std::vector<T>>& solution = std::nullopt);
};

template<std::floating_point T, std::integral I>
spatial_assembler_1d<T, I>::spatial_assembler_1d(
    finite_element_matrix_1d<T, I>& matrix, const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
    const std::optional<utils::nodes_sequence>& nodes_to_assemble)
    : _base{ matrix, mesh, nodes_to_assemble } {}

template<std::floating_point T, std::integral I>
T spatial_assembler_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{ 0 };
    const auto& el = _base::mesh().element();
    for (const size_t q : std::ranges::iota_view{ 0u, el.qnodes_count() })
        integral += el.weight(q) * el.qN(i, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<std::floating_point T, std::integral I>
T spatial_assembler_1d<T, I>::integrate_local(const coefficient_t<T, 1>& coefficient, const size_t e, const size_t i,
                                              const size_t j) const {
    T integral = T{ 0 };
    const auto& el = _base::mesh().element();
    for (const size_t q : _base::mesh().element().qnodes())
        integral += _base::evaluate(coefficient, e, q, _solution) * el.weight(q) * el.qNxi(i, q) * el.qNxi(j, q);
    return integral / _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<std::floating_point T, std::integral I>
T spatial_assembler_1d<T, I>::integrate_nonlocal(const coefficient_t<T, 1>& coefficient,
                                                 const std::function<T(const T&, const T&)>& influence,
                                                 const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    T integral = T{ 0 };
    const auto& el = _base::mesh().element();
    const auto qnodes = el.qnodes();
    for (const size_t qL : qnodes) {
        T inner_integral = T{ 0 };
        const T qcoordL = _base::mesh().qnode_coord(eL, qL);
        for (const size_t qNL : qnodes) {
            const T qcoordNL = _base::mesh().qnode_coord(eNL, qNL);
            inner_integral += influence(qcoordL, qcoordNL) * _base::evaluate(coefficient, eNL, qNL, _solution) *
                              el.weight(qNL) * el.qNxi(jNL, qNL);
        }
        integral += el.weight(qL) * el.qNxi(iL, qL) * inner_integral;
    }
    return integral;
}

template<std::floating_point T, std::integral I>
void spatial_assembler_1d<T, I>::integral_condition() {
#pragma omp parallel for default(none)
    for (size_t node = 0; node < _base::mesh().nodes_count(); ++node) {
        T& val = _base::matrix().inner.coeffRef(node, _base::mesh().nodes_count());
        for (const auto& [e, i] : _base::mesh().node_elements(node).to_array())
            if (e)
                val += integrate_basic(e, i);
    }
}

template<std::floating_point T, std::integral I>
template<class Parameters, class Get_Coefficient>
void spatial_assembler_1d<T, I>::calc_matrix(const Parameters& parameters, const problem_settings& settings,
                                             const Get_Coefficient get_coefficient,
                                             const std::optional<std::vector<T>>& solution) {
    if (parameters.size() != _base::mesh().segments_count())
        throw std::domain_error{ "The number of segments and the number of material parameters do not match." };
    if (settings.is_neumann)
        integral_condition();
    if (!settings.is_symmetric())
        _base::matrix().inner =
            Eigen::SparseMatrix<T, Eigen::RowMajor, I>(_base::matrix().inner.template selfadjointView<Eigen::Upper>());
    if (solution)
        _solution = *solution;

    _base::template calc_matrix(
        settings,
        [this, &parameters, get_coefficient](const size_t segment, const size_t e, const size_t i, const size_t j) {
            const auto& [model, physic] = parameters[segment];
            return model.local_weight * integrate_local(std::invoke(get_coefficient, physic), e, i, j);
        },
        [this, &parameters, get_coefficient](const size_t segment, const size_t eL, const size_t eNL,
                                             const size_t iL, const size_t jNL) {
            const auto& [model, physic] = parameters[segment];
            return nonlocal_weight(model.local_weight) *
                   integrate_nonlocal(std::invoke(get_coefficient, physic), model.influence, eL, eNL, iL, jNL);
        });
}

} // namespace nonlocal::solver_1d
