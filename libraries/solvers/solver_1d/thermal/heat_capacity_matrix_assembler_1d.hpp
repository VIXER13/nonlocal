#pragma once

#include "../../equation_parameters.hpp"

#include "finite_element_matrix_1d.hpp"
#include "thermal_parameters_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class heat_capacity_matrix_assembler_1d : public matrix_assembler_base<T, I> {
    using _base = finite_element_matrix_1d<T, I>;

    static std::vector<T> calc_factors(const nonlocal::thermal::parameters_1d<T>& parameters);

protected:
    T integrate(const size_t e, const size_t i, const size_t j) const;

public:
    explicit heat_capacity_matrix_assembler_1d(finite_element_matrix_1d<T, I>& matrix,
                                               const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~heat_capacity_matrix_assembler_1d() override = default;

    void compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                     const std::array<bool, 2> is_first_kind);
};

template<class T, class I>
std::vector<T> heat_capacity_matrix_assembler_1d<T, I>::calc_factors(const nonlocal::thermal::parameters_1d<T>& parameters) {
    std::vector<T> factors(parameters.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.size()})
        factors[i] = parameters[i].physical->capacity * parameters[i].physical->density;
    return factors;
}

template<class T, class I>
heat_capacity_matrix_assembler_1d<T, I>::heat_capacity_matrix_assembler_1d(finite_element_matrix_1d<T, I>& matrix,
                                                                           const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _base{matrix, mesh} {}

template<class T, class I>
T heat_capacity_matrix_assembler_1d<T, I>::integrate(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q) * el.qN(j, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
void heat_capacity_matrix_assembler_1d<T, I>::compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                                                      const std::array<bool, 2> is_first_kind) {
    if (parameters.size() != _base::mesh().segments_count())
        throw std::domain_error{"The number of segments and the number of material parameters do not match."};
    const std::vector<theory_t> local_theories(_base::mesh().segments_count(), theory_t::LOCAL);
    static constexpr bool SYMMETRIC = true;
    static constexpr bool NEUMANN = false;
    _base::template calc_matrix(local_theories, is_first_kind, SYMMETRIC, NEUMANN,
        [this, factors = calc_factors(parameters)](const size_t segment, const size_t e, const size_t i, const size_t j) {
            return factors[segment] * integrate(e, i, j);
        },
        [](const size_t, const size_t, const size_t, const size_t, const size_t) constexpr noexcept { return 0; }
    );
}

}