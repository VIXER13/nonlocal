#pragma once

#include "assebmler_base.hpp"

namespace nonlocal::solver_1d {

template<std::floating_point T, std::integral I>
class weighted_mass_assembler_1d : public assembler_base_1d<T, I> {
    using _base = assembler_base_1d<T, I>;

    std::vector<T> _solution;

    T integrate_mass(const coefficient_t<T, 1>& coefficient, const size_t e, const size_t i, const size_t j) const;
    T integrate_mass(const T factor, const size_t e, const size_t i, const size_t j) const;

  public:
    explicit weighted_mass_assembler_1d(finite_element_matrix_1d<T, I>& matrix,
                                        const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                        const std::optional<utils::nodes_sequence>& nodes_to_assemble = std::nullopt);
    ~weighted_mass_assembler_1d() override = default;

    template<class Parameters, class Get_Coefficient>
    void calc_matrix(const Parameters& parameters, const std::array<bool, 2>& is_first_kind,
                     const Get_Coefficient get_coefficient,
                     const std::optional<std::vector<T>>& solution = std::nullopt);
};

template<std::floating_point T, std::integral I>
weighted_mass_assembler_1d<T, I>::weighted_mass_assembler_1d(
    finite_element_matrix_1d<T, I>& matrix, const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
    const std::optional<utils::nodes_sequence>& nodes_to_assemble)
    : _base{ matrix, mesh, nodes_to_assemble } {}

template<std::floating_point T, std::integral I>
T weighted_mass_assembler_1d<T, I>::integrate_mass(const coefficient_t<T, 1>& coefficient, const size_t e, const size_t i,
                                                   const size_t j) const {
    T integral = T{ 0 };
    const auto& el = _base::mesh().element();
    for (const size_t q : el.qnodes())
        integral += _base::evaluate(coefficient, e, q, _solution) * el.weight(q) * el.qN(i, q) * el.qN(j, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<std::floating_point T, std::integral I>
T weighted_mass_assembler_1d<T, I>::integrate_mass(const T factor, const size_t e, const size_t i, const size_t j) const {
    return factor * integrate_mass(coefficient_t<T, 1>{ T{ 1 } }, e, i, j);
}

template<std::floating_point T, std::integral I>
template<class Parameters, class Get_Coefficient>
void weighted_mass_assembler_1d<T, I>::calc_matrix(const Parameters& parameters,
                                                    const std::array<bool, 2>& is_first_kind,
                                                    const Get_Coefficient get_coefficient,
                                                    const std::optional<std::vector<T>>& solution) {
    if (parameters.size() != _base::mesh().segments_count())
        throw std::domain_error{ "The number of segments and the number of material parameters do not match." };
    if (solution)
        _solution = *solution;
    const problem_settings settings = { .theories = std::vector<theory_t>(_base::mesh().segments_count(), theory_t::LOCAL),
                                        .is_first_kind = is_first_kind };
    _base::template calc_matrix(
        settings, [this, &parameters, get_coefficient](const size_t segment, const size_t e, const size_t i, const size_t j) {
            return integrate_mass(std::invoke(get_coefficient, parameters[segment].physical), e, i, j);
        });
}

} // namespace nonlocal::solver_1d
