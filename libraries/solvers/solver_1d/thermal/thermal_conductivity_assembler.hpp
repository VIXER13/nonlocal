#pragma once

#include "thermal_parameters_1d.hpp"

#include <solvers/solver_1d/base/spatial_assembler_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<std::floating_point T, std::integral I>
class thermal_conductivity_assembler_1d final : public spatial_assembler_1d<T, I> {
    using _base = spatial_assembler_1d<T, I>;

  public:
    using _base::_base;

    void calc_matrix(const parameters_1d<T>& parameters, const problem_settings& settings,
                     const std::optional<std::vector<T>>& solution = std::nullopt) {
        _base::calc_matrix(parameters, settings, &parameter_1d<T>::conductivity, solution);
    }
};

} // namespace nonlocal::solver_1d::thermal
