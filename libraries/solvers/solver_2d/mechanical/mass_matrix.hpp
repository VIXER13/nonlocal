#pragma once

#include <solvers/solver_2d/base/matrix_assembler.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
class mass_matrix : public matrix_assembler_base<metamath::linear::square_matrix<T, 2>> {
    using _base = matrix_assembler_base<metamath::linear::square_matrix<T, 2>>;

    static void throw_if_monostate_density(const evaluated_mechanical_parameters<T>& parameters) {
        for (const auto& [group, parameter] : parameters)
            if (std::holds_alternative<std::monostate>(parameter.physical.density))
                throw std::domain_error{"Density for group \"" + group + "\" is not defined!"};
    }

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
        T integral = 0;
        const auto& el = _base::mesh().container().element_2d(e);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * el.qN(j, q) * _base::mesh().jacobian(e, q);
        return integral;
    }

    T integrate_basic_pair(const metamath::types::vector_with_shifted_index<T>& density,
                           const size_t e, const size_t i, const size_t j) const {
        T integral = 0;
        const auto& el = _base::mesh().container().element_2d(e);
        const size_t qshift = _base::mesh().quad_shift(e);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += density[qshift + q] * el.weight(q) * el.qN(i, q) * el.qN(j, q) * _base::mesh().jacobian(e, q);
        return integral;
    }

public:
    explicit mass_matrix(const std::shared_ptr<mesh::mesh_2d<T>>& mesh) : _base{mesh} {}
    ~mass_matrix() noexcept override = default;

    void compute(const evaluated_mechanical_parameters<T>& parameters, const problem_settings& settings) {
        logger::info() << "Mass matrix assembly started" << std::endl;
        throw_if_monostate_density(parameters);
        _base::matrix().clear();
        _base::matrix().portrait.set_size(_base::rows(), _base::mesh().container().nodes_count());
        _base::init_shifts(settings);
        _base::init_indices(settings);
        _base::calc_coeffs(settings,
            [this, &parameters](const std::string& group, const size_t e, const size_t i, const size_t j) {
                const auto& density = std::get<evaluated_parameters<T>>(parameters.at(group).physical.density);
                const T integral = std::visit(metamath::types::visitor{
                    [this, e, i, j](const T density) { return density * integrate_basic_pair(e, i, j); },
                    [this, e, i, j](const auto& density) { return integrate_basic_pair(density, e, i, j); }
                }, density);
                return metamath::linear::square_matrix<T, 2>{integral, 0, 0, integral};
            },
            [](const std::string&, const size_t, const size_t, const size_t, const size_t) constexpr noexcept { return metamath::linear::square_matrix<T, 2>{}; }
        );
        logger::info() << "Mass matrix assembly finished" << std::endl;
    }
};

}