#pragma once

#include "thermal_parameters_1d.hpp"

#include <solvers/base/equation_parameters.hpp>
#include <solvers/solver_1d/base/assebmler_base.hpp>

#include <optional>

namespace nonlocal::thermal {

template<class T, class I>
class thermal_conductivity_assembler_1d final : public assembler_base_1d<T, I> {
    using _base = assembler_base_1d<T, I>;

    std::vector<T> _solution; // stub for nonlinear problems

    T evaluate(const coefficient_t<T, 1>& conductivity, const size_t e, const size_t q) const;

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_local(const coefficient_t<T, 1>& conductivity, const size_t e, const size_t i, const size_t j) const;
    T integrate_nonlocal(const coefficient_t<T, 1>& conductivity, const std::function<T(const T&, const T&)>& influence,
                         const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

    void integral_condition(); // for Neumann problem

public:
    explicit thermal_conductivity_assembler_1d(finite_element_matrix_1d<T, I>& matrix, const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                               const std::optional<utils::nodes_sequence>& nodes_to_assemble = std::nullopt);
    ~thermal_conductivity_assembler_1d() override = default;

    void calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2> is_first_kind,
                     const bool is_neumann = false, const bool is_symmetric = true,
                     const std::optional<std::vector<T>>& solution = std::nullopt);
};

template<class T, class I>
thermal_conductivity_assembler_1d<T, I>::thermal_conductivity_assembler_1d(finite_element_matrix_1d<T, I>& matrix, 
                                                                           const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                           const std::optional<utils::nodes_sequence>& nodes_to_assemble)
    : _base{matrix, mesh, nodes_to_assemble} {}

template<class T, class I>
T thermal_conductivity_assembler_1d<T, I>::evaluate(const coefficient_t<T, 1>& conductivity, const size_t e, const size_t q) const {
    return std::visit(metamath::visitor{
        [](const T value) noexcept { return value; },
        [this, e, q](const spatial_dependency<T, 1u>& value) { return value(_base::mesh().qnode_coord(e, q)); },
        [this, e, q](const solution_dependency<T, 1u>& value) { 
            const size_t qshift = _base::mesh().qnode_number(e, q);
            return value(_base::mesh().qnode_coord(e, q), _solution[qshift]); 
        }
    }, conductivity);
}

template<class T, class I>
T thermal_conductivity_assembler_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
T thermal_conductivity_assembler_1d<T, I>::integrate_local(const coefficient_t<T, 1>& conductivity, const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : _base::mesh().element().qnodes())
        integral += evaluate(conductivity, e, q) * el.weight(q) * el.qNxi(i, q) * el.qNxi(j, q);
    return integral / _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
T thermal_conductivity_assembler_1d<T, I>::integrate_nonlocal(const coefficient_t<T, 1>& conductivity, 
                                                              const std::function<T(const T&, const T&)>& influence,
                                                              const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    const auto qnodes = el.qnodes();
    for(const size_t qL : qnodes) {
        T inner_integral = T{0};
        const size_t qcoordL = _base::mesh().qnode_coord(eL, qL);
        for (const size_t qNL : qnodes) {
            const size_t qcoordNL = _base::mesh().qnode_coord(eNL, qNL);
            inner_integral += influence(qcoordL, qcoordNL) * evaluate(conductivity, eNL, qNL) * el.weight(qNL) * el.qNxi(jNL, qNL);
        }
        integral += el.weight(qL) * el.qNxi(iL, qL) * inner_integral;
    }
    return integral;
}

template<class T, class I>
void thermal_conductivity_assembler_1d<T, I>::integral_condition() {
#pragma omp parallel for default(none)
    for(size_t node = 0; node < _base::mesh().nodes_count(); ++node) {
        T& val = _base::matrix().inner.coeffRef(node, _base::mesh().nodes_count());
        for(const auto& [e, i] : _base::mesh().node_elements(node).to_array())
            if (e) 
                val += integrate_basic(e, i);
    }
}

template<class T, class I>
void thermal_conductivity_assembler_1d<T, I>::calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2> is_first_kind,
                                                       const bool is_neumann, const bool is_symmetric,
                                                       const std::optional<std::vector<T>>& solution) {
    if (parameters.size() != _base::mesh().segments_count())
       throw std::domain_error{"The number of segments and the number of material parameters do not match."};
    if (is_neumann)
        integral_condition();
    if (!is_symmetric)
        _base::matrix().inner = Eigen::SparseMatrix<T, Eigen::RowMajor, I>(_base::matrix().inner.template selfadjointView<Eigen::Upper>());

    _base::template calc_matrix(theories_types(parameters), is_first_kind, is_symmetric,
        [this, &parameters, &solution](const size_t segment, const size_t e, const size_t i, const size_t j) {
            const auto& [model, physic] = parameters[segment];
            return model.local_weight * integrate_local(physic.conductivity, e, i, j);
        },
        [this, &parameters, &solution](const size_t segment, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            const auto& [model, physic] = parameters[segment];
            return nonlocal_weight(model.local_weight) * integrate_nonlocal(physic.conductivity, model.influence, eL, eNL, iL, jNL);
        }
    );
}

}