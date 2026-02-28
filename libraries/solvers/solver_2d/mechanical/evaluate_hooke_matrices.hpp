#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T, std::integral I>
evaluated_hook_matrices_2d<T> evaluate_hooke_matrices(const mesh::mesh_2d<T, I>& mesh, 
                                                      const elastic_parameters<T>& parameters) {
    evaluated_hook_matrices_2d<T> result;
    for (const auto& pair : parameters) {
        const auto& name = pair.first;
        const auto& parameter = pair.second;

        auto hook_matrices = std::visit([&mesh, &name](const auto& elastic) -> evaluated_hook_matrix_t<T> {
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(elastic)>, isotropic_elastic_parameters<T>>) {
                if (is_constant(elastic.young_modulus) && is_constant(elastic.poissons_ratio))
                    return elastic.hooke({});
            } else if constexpr (std::is_same_v<std::remove_cvref_t<decltype(elastic)>, orthotropic_elastic_parameters<T>>) {
                if (is_constant(elastic.young_modulus) && 
                    is_constant(elastic.poissons_ratio) && 
                    is_constant(elastic.shear_modulus))
                    return elastic.hooke({});
            } else if constexpr (std::is_same_v<std::remove_cvref_t<decltype(elastic)>, anisotropic_elastic_parameters<T>>) {
                if (is_constant(elastic.main_parameters.young_modulus) && 
                    is_constant(elastic.main_parameters.poissons_ratio) && 
                    is_constant(elastic.main_parameters.shear_modulus))
                    return elastic.hooke({});
            } else
                static_assert(false, "Unknown material type");

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