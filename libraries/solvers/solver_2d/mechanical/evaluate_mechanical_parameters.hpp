#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T, std::integral I>
evaluated_mechanical_parameters<T> evaluate_mechanical_parameters(const mesh::mesh_2d<T, I>& mesh, 
                                                                  const raw_mechanical_parameters<T>& parameters) {
    evaluated_mechanical_parameters<T> result;
    for (const auto& [name, parameter] : parameters) {
        auto& [_, physical] = result[name] = { .model = parameter.model };

        physical.elastic = std::visit([&mesh, &name](const auto& elastic) -> evaluated_hook_matrix_t<T> {
            if (elastic.is_constant())
                return elastic.hooke({});
            const auto qshifts = mesh.quad_shifts(name);
            static constexpr size_t N = std::tuple_size_v<std::remove_cvref_t<decltype(elastic.hooke({}))>>;
            metamath::types::vector_with_shifted_index<std::array<T, N>> result = {
                .container = std::vector<std::array<T, N>>(qshifts.size()),
                .shift = qshifts.front()
            };
            for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                result[q] = elastic.hooke(mesh.quad_coord(q));
            return result;
        }, parameter.physical.elastic);

        if (parameter.physical.thermal_expansion.valueless_by_exception())
            continue;
        
        physical.thermal_expansion = std::visit(metamath::types::visitor{
            [&mesh, &name](const raw_isotropic_thermal_expansion_t<T>& thermal_expansion) -> evaluated_thermal_expansion_t<T> {
                if (is_constant(thermal_expansion))
                    return std::get<T>(thermal_expansion);
                std::vector<T> result;
                const auto qshifts = mesh.quad_shifts(name);
                result.reserve(qshifts.size());
                for(const size_t qshift: qshifts)
                    result.push_back(evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(qshift), {}));
                return metamath::types::vector_with_shifted_index<T>{std::move(result), qshifts.front()};
            },
            [&mesh, &name](const auto& thermal_expansion) -> evaluated_thermal_expansion_t<T> {
                if (is_constant(thermal_expansion))
                    return evaluate<T, 2u>(thermal_expansion, {}, {});
                const auto qshifts = mesh.quad_shifts(name);
                static constexpr size_t N = std::tuple_size_v<std::remove_cvref_t<decltype(thermal_expansion)>>;
                metamath::types::vector_with_shifted_index<std::array<T, N>> result = {
                    .container = std::vector<std::array<T, N>>(qshifts.size()),
                    .shift = qshifts.front()
                };
                for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                    result[q] = evaluate<T, 2u>(thermal_expansion, mesh.quad_coord(q), {});
                return result;
            }
        }, parameter.physical.thermal_expansion);
    }
    return result;
}

}