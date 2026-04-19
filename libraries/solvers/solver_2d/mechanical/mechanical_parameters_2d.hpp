#pragma once

#include "elastic_parameters.hpp"
#include "thermal_expansion_parameters.hpp"

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
struct raw_mechanical_parameters_t final {
    elastic_parameters_t<T> elastic;
    std::optional<raw_thermal_expansion_t<T>> thermal_expansion;
};

template<std::floating_point T>
struct evaluated_mechanical_parameters_t final {
    evaluated_hook_matrix_t<T> elastic;
    evaluated_thermal_expansion_t<T> thermal_expansion;
    std::optional<evaluated_thermal_strain<T>> thermal_strain;
};

template<std::floating_point T>
using raw_mechanical_parameters = std::unordered_map<std::string, equation_parameters<2, T, raw_mechanical_parameters_t>>;

template<std::floating_point T>
using evaluated_mechanical_parameters = std::unordered_map<std::string, equation_parameters<2, T, evaluated_mechanical_parameters_t>>;

template<std::floating_point T>
std::array<T, 3> calc_stress(const hooke_matrix_t<T>& hooke, const std::array<T, 3>& strain) noexcept {
    return std::visit(metamath::types::visitor{
        [&strain](const isotropic_hook_matrix_t<T>& hooke) -> std::array<T, 3> {
            using namespace isotropic_indices;
            return {hooke[_11] * strain[XX] + hooke[_12] * strain[YY],
                    hooke[_12] * strain[XX] + hooke[_11] * strain[YY],
                2 * hooke[_66] * strain[XY]};
        },
        [&strain](const orthotropic_hook_matrix_t<T>& hooke) -> std::array<T, 3> {
            using namespace orthotropic_indices;
            return {hooke[_11] * strain[XX] + hooke[_12] * strain[YY],
                    hooke[_12] * strain[XX] + hooke[_22] * strain[YY],
                2 * hooke[_66] * strain[XY]};
        },
        [&strain](const anisotropic_hook_matrix_t<T>& hooke) -> std::array<T, 3> {
            using namespace anisotropic_indices;
            return {hooke[_11] * strain[XX] + hooke[_12] * strain[YY] + 2 * hooke[_16] * strain[XY],
                    hooke[_12] * strain[XX] + hooke[_22] * strain[YY] + 2 * hooke[_26] * strain[XY],
                    hooke[_16] * strain[XX] + hooke[_26] * strain[YY] + 2 * hooke[_66] * strain[XY]};
        }
    }, hooke);
}

}