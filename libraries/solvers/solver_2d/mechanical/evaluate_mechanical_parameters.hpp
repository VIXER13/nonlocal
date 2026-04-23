#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T, std::integral I>
evaluated_mechanical_parameters<T> evaluate_mechanical_parameters(const mesh::mesh_2d<T, I>& mesh, 
                                                                  const raw_mechanical_parameters<T>& parameters,
                                                                  const std::vector<T>& delta_temperature = {}) {
    evaluated_mechanical_parameters<T> result;
    const auto temperature = delta_temperature.empty() ? delta_temperature : nonlocal::mesh::utils::nodes_to_qnodes(mesh, delta_temperature);
    for (const auto& [name, parameter] : parameters) {
        auto& [_, physical] = result[name] = { .model = parameter.model };

        physical.elastic = std::visit([&mesh, &name](const auto& elastic) -> evaluated_hook_matrix_t<T> {
            if (elastic.is_constant())
                return elastic.hooke({});
            const auto qshifts = mesh.quad_shifts(name);
            using Hooke = std::remove_cvref_t<decltype(elastic.hooke({}))>;
            metamath::types::vector_with_shifted_index<Hooke> result = {
                .container = std::vector<Hooke>(qshifts.size()),
                .shift = qshifts.front()
            };
            for(const size_t qshift : qshifts)
                result[qshift] = elastic.hooke(mesh.quad_coord(qshift));
            return result;
        }, parameter.physical.elastic);

        if (temperature.empty())
            continue;

        physical.thermal_strain = std::visit(metamath::types::visitor{
            [&mesh, &name, &temperature](const raw_isotropic_thermal_expansion_t<T>& thermal_expansion) -> evaluated_thermal_strain_t<T> {
                const auto qshifts = mesh.quad_shifts(name);
                isotropic_thermal_strain<T> result{.strain = {std::vector<T>(qshifts.size(), T{0}), qshifts.front()}};
                for(const size_t qshift : qshifts)
                    result.strain[qshift] = temperature[qshift] * evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(qshift), {});
                return result;
            },
            [&mesh, &name, &temperature](const auto& thermal_expansion) -> evaluated_thermal_strain_t<T> {
                using Expansion = std::remove_cvref_t<decltype(thermal_expansion)>;
                if constexpr (std::is_same_v<Expansion, std::monostate>)
                    return std::monostate{};
                else {
                    const auto qshifts = mesh.quad_shifts(name);
                    if (is_constant(thermal_expansion)) {
                        std::conditional_t<std::is_same_v<Expansion, raw_orthotropic_thermal_expansion_t<T>>,
                            orthotropic_constant_thermal_strain<T>,
                            anisotropic_constant_thermal_strain<T>
                        > result = {
                            .delta_temperature = metamath::types::vector_with_shifted_index<T>{
                                .container = std::vector<T>(std::next(temperature.begin(), qshifts.front()), 
                                                            std::next(temperature.begin(), qshifts.back() + 1)),
                                .shift = qshifts.front()
                            },
                            .thermal_expansion = evaluate<T, 2u>(thermal_expansion, {}, {})
                        };
                        return result;
                    }
                    using namespace metamath::functions;
                    std::conditional_t<std::is_same_v<Expansion, raw_orthotropic_thermal_expansion_t<T>>,
                        orthotropic_thermal_strain<T>,
                        anisotropic_thermal_strain<T>
                    > result = {
                        .strain = {
                            .container = std::vector<std::array<T, std::tuple_size_v<Expansion>>>(qshifts.size()),
                            .shift = qshifts.front()
                        } 
                    };
                    for(const size_t qshift : qshifts)
                        result.strain[qshift] = temperature[qshift] * evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(qshift), {});
                    return result;
                }
            }
        }, parameter.physical.thermal_expansion);
    }
    return result;
}

}