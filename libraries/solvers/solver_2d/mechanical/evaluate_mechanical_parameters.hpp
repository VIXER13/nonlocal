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
            for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                result[q] = elastic.hooke(mesh.quad_coord(q));
            return result;
        }, parameter.physical.elastic);

        if (temperature.empty() || std::holds_alternative<std::monostate>(parameter.physical.thermal_expansion))
            continue;

        physical.thermal_strain = { std::visit(metamath::types::visitor{
            [&mesh, &name, &temperature](const raw_isotropic_thermal_expansion_t<T>& thermal_expansion) -> evaluated_thermal_strain_t<T> {
                const auto qshifts = mesh.quad_shifts(name);
                isotropic_thermal_strain_t<T> result{std::vector<T>(qshifts.size(), T{0}), qshifts.front()};
                for(const size_t qshift : qshifts)
                    result[qshift] = temperature[qshift];
                if (is_constant(thermal_expansion)) {
                    using namespace metamath::functions;
                    result.container *= std::get<T>(thermal_expansion);
                } else {
                    for(const size_t qshift : qshifts)
                        result[qshift] *= evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(qshift), {});
                }
                return result;
            },
            [&mesh, &name, &temperature](const raw_orthotropic_thermal_expansion_t<T>& thermal_expansion) -> evaluated_thermal_strain_t<T> {
                const auto qshifts = mesh.quad_shifts(name);
                if (is_constant(thermal_expansion)) {
                    orthotropic_constant_thermal_strain_t<T> result{
                        metamath::types::vector_with_shifted_index<T>{std::vector<T>(qshifts.size(), T{0}), qshifts.front()},
                        evaluate<T, 2u>(thermal_expansion, {}, {})
                    };
                    for(const size_t qshift : qshifts)
                        result.first[qshift] = temperature[qshift];
                    return result;
                }
                using namespace metamath::functions;
                orthotropic_thermal_strain_t<T> result{ std::vector<std::array<T, 2>>(qshifts.size()), qshifts.front() };
                for(const size_t qshift : qshifts)
                    result[qshift] = temperature[qshift] * evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(qshift), {});
                return result;
            },
            [&mesh, &name, &temperature](const raw_anisotropic_thermal_expansion_t<T>& thermal_expansion) -> evaluated_thermal_strain_t<T> {
                const auto qshifts = mesh.quad_shifts(name);
                if (is_constant(thermal_expansion)) {
                    anisotropic_constant_thermal_strain_t<T> result{
                        metamath::types::vector_with_shifted_index<T>{std::vector<T>(qshifts.size(), T{0}), qshifts.front()},
                        evaluate<T, 2u>(thermal_expansion, {}, {})
                    };
                    for(const size_t qshift : qshifts)
                        result.first[qshift] = temperature[qshift];
                    return result;
                }
                using namespace metamath::functions;
                anisotropic_thermal_strain_t<T> result{ std::vector<std::array<T, 3>>(qshifts.size()), qshifts.front() };
                for(const size_t qshift : qshifts)
                    result[qshift] = temperature[qshift] * evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(qshift), {});
                return result;
            },
            [](auto&&) -> evaluated_thermal_strain_t<T> { throw std::domain_error{"Unexpected thermal_expansion variant state."}; }
        }, parameter.physical.thermal_expansion) };
    }
    return result;
}

}