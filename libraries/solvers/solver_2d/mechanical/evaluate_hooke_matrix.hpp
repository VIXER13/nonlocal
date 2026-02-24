#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T, std::integral I>
evaluated_hook_matrix_t<T> evaluate_hooke_matrix(const mesh::mesh_2d<T, I>& mesh, 
                                                 const parameters_2d<T>& parameters) {
    evaluated_hook_matrix_t<T> result;
    for (const auto& pair : parameters.materials) {
        const auto& name = pair.first;
        const auto& parameter = pair.second;

        auto hook_matrices = std::visit([&mesh, &name](const auto& hook_matrix) -> evaluated_hook_matrix_t<T> {
            if (is_constant<T, 2u>(hook_matrix))
                return evaluate<T, 2u>(hook_matrix, {}, {});
            const auto qshifts = mesh.quad_shifts(name);
            static constexpr size_t N = decltype(hook_matrix){}.size();
            metamath::types::vector_with_shifted_index<std::array<T, N>> result = {
                .container = std::vector<std::array<T, N>>(qshifts.size()),
                .shift = qshifts.front()
            };
            for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                result[q] = evaluate<T, 2u>(hook_matrix, mesh.quad_coord(q), T{0});
            return result;
        }, parameter.physical.hook());

        result[name] = {
            .model = parameter.model;
            .physical = std::move(hook_matrices)
        };
    }
    return result;
}

}