#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T, std::integral I>
evaluated_hook_matrices_2d<T> evaluate_hooke_matrices(const mesh::mesh_2d<T, I>& mesh, 
                                                      const elastic_parameters<T>& parameters) {
    evaluated_hook_matrices_2d<T> result;
    for (const auto& [name, parameter] : parameters) {

        auto hook_matrices = std::visit([&mesh, &name](const auto& elastic) -> evaluated_hook_matrix_t<T> {
            if (elastic.is_constant())
                return elastic.hooke({});
            const auto qshifts = mesh.quad_shifts(name);
            static constexpr size_t N = decltype(elastic.hooke({})){}.size();
            metamath::types::vector_with_shifted_index<std::array<T, N>> result = {
                .container = std::vector<std::array<T, N>>(qshifts.size()),
                .shift = qshifts.front()
            };
            for(size_t q = qshifts.front(); q <= qshifts.back(); ++q)
                result[q] = elastic.hooke(mesh.quad_coord(q));
            return result;
        }, parameter.physical);

        result[name] = {
            .model = parameter.model,
            .physical = std::move(hook_matrices)
        };
    }
    return result;
}

}