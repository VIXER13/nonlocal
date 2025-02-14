#pragma once

#include "thermal_boundary_conditions_1d.hpp"
#include <Eigen/Sparse>
#include <ranges>

namespace nonlocal::thermal {

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                            const thermal_boundary_condition_1d<T>& boundary_condition,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev,
                            const T time_step, const size_t index) {
    static constexpr T STEFAN_BOLTZMANN_CONSTANT_X4 = T{4} * STEFAN_BOLTZMANN_CONSTANT<T>;
    static constexpr T STEFAN_BOLTZMANN_CONSTANT_X3 = T{3} * STEFAN_BOLTZMANN_CONSTANT<T>;
    if (const auto* const condition = dynamic_cast<const radiation_1d<T>*>(&boundary_condition)) {
        const T temperature_pow_3 = metamath::functions::power<3>(temperature_prev[index]);
        matrix.coeffRef(index, index) += time_step * condition->emissivity() * STEFAN_BOLTZMANN_CONSTANT_X4 * temperature_pow_3;
        f[index] += time_step * condition->emissivity() * STEFAN_BOLTZMANN_CONSTANT_X3 * temperature_pow_3 * temperature_prev[index];
    }
}

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                            const thermal_boundaries_conditions_1d<T>& boundaries_conditions, 
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev, 
                            const T time_step) {
    const std::array<size_t, 2> indices = {0, size_t(matrix.outerSize() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        radiation_condition_1d(matrix, f, *boundaries_conditions[b], temperature_prev, time_step, indices[b]);
}

}