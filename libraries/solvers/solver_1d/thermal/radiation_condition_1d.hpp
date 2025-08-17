#pragma once

#include "thermal_boundary_conditions_1d.hpp"
#include <Eigen/Sparse>
#include <ranges>

namespace nonlocal::solver_1d::thermal {

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                            const thermal_boundary_condition_1d<T>& boundary_condition,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev,
                            const T time_step, const size_t index) {
    static constexpr T Stefan_Boltzmann_Constant_X4 = T{4} * Stefan_Boltzmann_Constant<T>;
    static constexpr T Stefan_Boltzmann_Constant_X3 = T{3} * Stefan_Boltzmann_Constant<T>;
    if (const auto* const condition = dynamic_cast<const radiation_1d<T>*>(&boundary_condition)) {
        const T temperature_pow_3 = metamath::functions::power<3>(temperature_prev[index]);
        matrix.coeffRef(index, index) += time_step * condition->emissivity() * Stefan_Boltzmann_Constant_X4 * temperature_pow_3;
        f[index] += time_step * condition->emissivity() * Stefan_Boltzmann_Constant_X3 * temperature_pow_3 * temperature_prev[index];
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

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& residual,
                            const thermal_boundary_condition_1d<T>& boundary_condition,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev,
                            const size_t index) {
    if (const auto* const condition = dynamic_cast<const radiation_1d<T>*>(&boundary_condition)) {
        static constexpr T Stefan_Boltzmann_Constant_X4 = T{4} * Stefan_Boltzmann_Constant<T>;
        const T temperature_pow_3 = metamath::functions::power<3>(temperature_prev[index]);
        matrix.coeffRef(index, index) -= Stefan_Boltzmann_Constant_X4 * condition->emissivity() * temperature_pow_3;
        residual[index] -= Stefan_Boltzmann_Constant<T> * condition->emissivity() * temperature_pow_3 * temperature_prev[index];
    }
}

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& residual,
                            const thermal_boundaries_conditions_1d<T>& boundaries_conditions, 
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev) {
    const std::array<size_t, 2> indices = {0, size_t(matrix.outerSize() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u}) {
        radiation_condition_1d(matrix, residual, *boundaries_conditions[b], temperature_prev, indices[b]);
    }
}

}