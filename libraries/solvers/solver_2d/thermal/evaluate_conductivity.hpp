#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T, std::integral I>
evaluated_conductivity_2d<T> evaluate_conductivity(const mesh::mesh_2d<T, I>& mesh, 
                                                   const parameters_2d<T>& parameters,
                                                   const std::vector<T>& solution) {
    evaluated_conductivity_2d<T> conductivity;
    for (const auto& pair : parameters) {
        const auto& name = pair.first;
        const auto& parameter = pair.second;

        auto conduct = std::visit(metamath::types::visitor{
            [&mesh, &name, &solution](const raw_isotropic_conductivity_t<T>& conductivity) -> evaluated_conductivity_t<T> {
                if (is_constant(conductivity))
                    return std::get<T>(conductivity);
                std::vector<T> result;
                const auto qshifts = mesh.quad_shifts(name);
                result.reserve(qshifts.size());
                for(const size_t qshift: qshifts)
                    result.push_back(evaluate<T, 2u>(conductivity, mesh.quad_coord(qshift), solution[qshift]));
                return metamath::types::vector_with_shifted_index<T>{std::move(result), qshifts.front()};
            },
            [&mesh, &name, &solution](const auto& conductivity) -> evaluated_conductivity_t<T> {
                if (is_constant(conductivity))
                    return evaluate<T, 2u>(conductivity, {}, {});
                const auto qshifts = mesh.quad_shifts(name);
                static constexpr size_t N = decltype(conductivity){}.size();
                metamath::types::vector_with_shifted_index<std::array<T, N>> result = {
                    .container = std::vector<std::array<T, N>>(qshifts.size()),
                    .shift = qshifts.front()
                };
                for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                    result[q] = evaluate<T, 2u>(conductivity, mesh.quad_coord(q), solution[q]);
                return result;
            }
        }, parameter.physical.conductivity);

        conductivity[name] = {
            .model = parameter.model,
            .physical = std::move(conduct)
        };
    }
    return conductivity;
}

}