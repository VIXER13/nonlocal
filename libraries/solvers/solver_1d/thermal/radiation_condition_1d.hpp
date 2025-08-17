#pragma once

#include "thermal_boundary_conditions_1d.hpp"
#include <Eigen/Sparse>
#include <ranges>

namespace nonlocal::solver_1d::thermal {

template<bool Is_Stationary, class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& right_part,
                            const thermal_boundaries_conditions_1d<T>& boundaries_conditions, 
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev, 
                            const T time_step = T{1}) {
    const std::array<size_t, 2> indices = {0, size_t(matrix.outerSize() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u}) {
        const size_t index = indices[b];
        if (const auto* const condition = dynamic_cast<const radiation_1d<T>*>(boundaries_conditions[b].get())) {
            const T value = time_step * condition->emissivity() * metamath::functions::power<3>(temperature_prev[index]);
            if constexpr (Is_Stationary) {
                matrix.coeffRef(index, index) -= T{4} * Stefan_Boltzmann_Constant<T> * value;
                // In a stationary calculation, the right-hand side refers to the residual
                right_part[index] -= value * Stefan_Boltzmann_Constant<T> * temperature_prev[index];
            } else {
                matrix.coeffRef(index, index) += T{4} * Stefan_Boltzmann_Constant<T> * value;
                right_part[index] += T{3} * Stefan_Boltzmann_Constant<T> * value * temperature_prev[index];
            }
        }
    }
}

}