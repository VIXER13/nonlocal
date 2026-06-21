#pragma once

#include <solvers/solver_1d/mechanical/mechanical_parameters_1d.hpp>
#include <solvers/solver_1d/base/assebmler_base.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T, std::integral I>
class mass_assembler_1d final : public assembler_base_1d<T, I> {
    using _base = assembler_base_1d<T, I>;

    std::vector<T> _solution;

    T evaluate(const coefficient_t<T, 1>& density, const size_t e, const size_t q) const;
    T integrate_mass(const coefficient_t<T, 1>& density, const size_t e, const size_t i, const size_t j) const;

public:
    explicit mass_assembler_1d(finite_element_matrix_1d<T, I>& matrix, const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                               const std::optional<utils::nodes_sequence>& nodes_to_assemble = std::nullopt);
    ~mass_assembler_1d() override = default;

    void calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2>& is_first_kind,
                     const std::optional<std::vector<T>>& solution = std::nullopt);
};

template<std::floating_point T, std::integral I>
mass_assembler_1d<T, I>::mass_assembler_1d(finite_element_matrix_1d<T, I>& matrix,
                                           const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                           const std::optional<utils::nodes_sequence>& nodes_to_assemble)
    : _base{matrix, mesh, nodes_to_assemble} {}

template<std::floating_point T, std::integral I>
T mass_assembler_1d<T, I>::evaluate(const coefficient_t<T, 1>& density, const size_t e, const size_t q) const {
    return std::visit(metamath::visitor{
        [](const T value) noexcept { return value; },
        [this, e, q](const spatial_dependency<T, 1u>& value) { return value(_base::mesh().qnode_coord(e, q)); },
        [this, e, q](const solution_dependency<T, 1u>& value) {
            const size_t qshift = _base::mesh().qnode_number(e, q);
            return value(_base::mesh().qnode_coord(e, q), _solution[qshift]);
        }
    }, density);
}

template<std::floating_point T, std::integral I>
T mass_assembler_1d<T, I>::integrate_mass(const coefficient_t<T, 1>& density, const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : el.qnodes())
        integral += evaluate(density, e, q) * el.weight(q) * el.qN(i, q) * el.qN(j, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<std::floating_point T, std::integral I>
void mass_assembler_1d<T, I>::calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2>& is_first_kind,
                                          const std::optional<std::vector<T>>& solution) {
    if (parameters.size() != _base::mesh().segments_count())
        throw std::domain_error{"The number of segments and the number of material parameters do not match."};
    if (solution)
        _solution = *solution;
    const problem_settings settings = {
        .theories = std::vector<theory_t>(_base::mesh().segments_count(), theory_t::LOCAL),
        .is_first_kind = is_first_kind
    };
    _base::template calc_matrix(settings,
        [this, &parameters](const size_t segment, const size_t e, const size_t i, const size_t j) {
            return integrate_mass(parameters[segment].physical.density, e, i, j);
        }
    );
}

}
