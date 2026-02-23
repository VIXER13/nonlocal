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
                if (is_constant<T, 2u>(conductivity))
                    return std::get<T>(conductivity);
    
                std::vector<T> result;
                const auto qshifts = mesh.quad_shifts(name);
                result.reserve(qshifts.size());
                for(const size_t qshift: qshifts)
                    result.push_back(evaluate<T, 2u>(conductivity, mesh.quad_coord(qshift), solution[qshift]));
                return metamath::types::vector_with_shifted_index<T>{std::move(result), qshifts.front()};
            },
            [&mesh, &name, &solution](const raw_orthotropic_conductivity_t<T>& conductivity) -> evaluated_conductivity_t<T> {
                if (is_constant<T, 2u>(conductivity[X]) && is_constant<T, 2u>(conductivity[Y]))
                    return std::array{std::get<T>(conductivity[X]), std::get<T>(conductivity[Y])};
    
                std::vector<std::array<T, 2>> result;
                const auto qshifts = mesh.quad_shifts(name);
                result.reserve(qshifts.size());
                for(const size_t qshift: qshifts)
                    result.push_back({
                        evaluate<T, 2u>(conductivity[X], mesh.quad_coord(qshift), solution[qshift]),
                        evaluate<T, 2u>(conductivity[Y], mesh.quad_coord(qshift), solution[qshift])
                    });
                return metamath::types::vector_with_shifted_index<std::array<T, 2>>{std::move(result), qshifts.front()};
            },
            [&mesh, &name, &solution](const raw_anisotropic_conductivity_t<T>& conductivity) -> evaluated_conductivity_t<T> {
                if (is_constant<T, 2u>(conductivity[XX]) && is_constant<T, 2u>(conductivity[YY]) && is_constant<T, 2u>(conductivity[XY]))
                    return std::array<T, 3>{std::get<T>(conductivity[XX]), std::get<T>(conductivity[YY]), std::get<T>(conductivity[XY])};
    
                std::vector<std::array<T, 3>> result;
                const auto qshifts = mesh.quad_shifts(name);
                result.reserve(qshifts.size());
                for(const size_t qshift: qshifts)
                    result.push_back({
                        evaluate<T, 2u>(conductivity[XX], mesh.quad_coord(qshift), solution[qshift]),
                        evaluate<T, 2u>(conductivity[YY], mesh.quad_coord(qshift), solution[qshift]),
                        evaluate<T, 2u>(conductivity[XY], mesh.quad_coord(qshift), solution[qshift])
                    });
                return metamath::types::vector_with_shifted_index<std::array<T, 3>>{std::move(result), qshifts.front()};
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