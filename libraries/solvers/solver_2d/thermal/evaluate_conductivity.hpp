#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T, std::integral I>
evaluated_conductivity_2d<T> evaluate_conductivity(const mesh::mesh_2d<T, I>& mesh, 
                                                   const parameters_2d<T>& parameters,
                                                   const std::vector<T>& solution) {
    evaluated_conductivity_2d<T> conductivity;
    for (const auto& [name, parameter] : parameters)
        conductivity[name] = {
            .model = parameter.model,
            .physical = std::visit([&mesh, &name, &solution](const auto& conductivity) -> evaluated_conductivity_t<T> {
                using Type = decltype(evaluate<T, 2u>(conductivity, {}, {}));
                if (is_constant(conductivity))
                    return evaluate<T, 2u>(conductivity, {}, {});
                const auto qshifts = mesh.quad_shifts(name);
                metamath::types::vector_with_shifted_index<Type> result = {
                    .container = std::vector<Type>(qshifts.size()),
                    .shift = qshifts.front()
                };
                for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                    result[q] = evaluate<T, 2u>(conductivity, mesh.quad_coord(q), solution[q]);
                return result;
            }, parameter.physical.conductivity)
        };
    return conductivity;
}

}