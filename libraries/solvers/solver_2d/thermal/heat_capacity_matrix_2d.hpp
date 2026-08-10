#pragma once

#include <solvers/solver_2d/base/matrix_assembler.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
class heat_capacity_matrix_2d : public matrix_assembler_base<T> {
    using _base = matrix_assembler_base<T>;

    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
        T integral = 0;
        const auto& el = _base::mesh().container().element_2d(e);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * el.qN(j, q) * _base::mesh().jacobian(e, q);
        return integral;
    }

public:
    explicit heat_capacity_matrix_2d(const mesh::mesh_2d<T>& mesh) : _base{mesh} {}
    ~heat_capacity_matrix_2d() noexcept override = default;

    void compute(const parameters_2d<T>& parameters, const problem_settings& settings) {
        logger::info() << "Capacity matrix assembly started" << std::endl;
        _base::matrix().clear();
        _base::matrix().portrait.set_size(_base::rows(), _base::mesh().container().nodes_count());
        _base::init_shifts(settings);
        _base::init_indices(settings);
        _base::calc_coeffs(settings,
            [this, &parameters](const std::string& group, const size_t e, const size_t i, const size_t j) {
                const auto& parameter = parameters.at(group).physical;
                return parameter.density * parameter.capacity * integrate_basic_pair(e, i, j); 
            },
            [](const std::string&, size_t, size_t, size_t, size_t) { return T{0}; }
        );
        logger::info() << "Capacity matrix assembly finished" << std::endl;
    }
};

}