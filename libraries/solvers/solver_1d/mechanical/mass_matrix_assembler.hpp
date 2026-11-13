#pragma once

#include "mechanical_parameters_1d.hpp"

#include <solvers/solver_1d/base/weighted_mass_assembler_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T, std::integral I>
class mass_assembler_1d final : public weighted_mass_assembler_1d<T, I> {
    using _base = weighted_mass_assembler_1d<T, I>;

  public:
    using _base::_base;

    void calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2>& is_first_kind,
                     const std::optional<std::vector<T>>& solution = std::nullopt) {
        _base::calc_matrix(parameters, is_first_kind, &parameter_1d<T>::density, solution);
    }
};

} // namespace nonlocal::solver_1d::mechanical
