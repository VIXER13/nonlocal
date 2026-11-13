#pragma once

#include "thermal_parameters_1d.hpp"

#include <solvers/solver_1d/base/weighted_mass_assembler_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<std::floating_point T, std::integral I>
class heat_capacity_assembler_1d : public weighted_mass_assembler_1d<T, I> {
    using _base = weighted_mass_assembler_1d<T, I>;

  public:
    using _base::_base;

    void calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2> is_first_kind) {
        if (parameters.size() != _base::mesh().segments_count())
            throw std::runtime_error{ "The number of segments and the number of material parameters do not match." };
        static constexpr auto Factor = [](const parameter_1d<T>& parameter) constexpr noexcept {
            return parameter.capacity * parameter.density;
        };
        _base::calc_matrix(parameters, is_first_kind, Factor);
    }
};

} // namespace nonlocal::solver_1d::thermal
