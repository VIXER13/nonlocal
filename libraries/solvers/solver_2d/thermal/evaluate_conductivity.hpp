#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
evaluated_thermal_parameters<T> evaluate_conductivity(const mesh::mesh_2d<T>& mesh, 
                                                      const raw_thermal_parameters<T>& parameters,
                                                      const std::vector<T>& solution) {
    evaluated_thermal_parameters<T> conductivity;
    for (const auto& [name, parameter] : parameters) {
        auto& param = conductivity[name] = {.model = parameter.model};

        param.physical.conductivity = std::visit([&mesh, &name, &solution](const auto& conductivity) -> evaluated_conductivity_t<T> {
            using Type = decltype(evaluate<T, 2u>(conductivity, {}, {}));
            if (is_constant(conductivity))
                return evaluate<T, 2u>(conductivity, {}, {});
            const auto qshifts = mesh.quad_shifts(name);
            metamath::types::vector_with_shifted_index<Type> result = {
                .container = std::vector<Type>(qshifts.size()),
                .shift = qshifts.front()
            };
            if (is_spatial(conductivity)) {
                for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                    result[q] = evaluate<T, 2u>(conductivity, mesh.quad_coord(q), {});
            } else {
                for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                    result[q] = evaluate<T, 2u>(conductivity, mesh.quad_coord(q), solution[q]);
            }
            return result;
        }, parameter.physical.conductivity);

        param.physical.capacity = parameter.physical.capacity;
        param.physical.density = parameter.physical.density;
        param.physical.relaxation_time = parameter.physical.relaxation_time;
    }
    return conductivity;
}

}