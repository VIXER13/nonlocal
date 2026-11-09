#pragma once

#include "boundary_condition_second_kind_1d.hpp"

#include <mesh/mesh_1d/mesh_1d.hpp>

#include <Eigen/Core>

#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <ranges>
#include <stdexcept>

namespace nonlocal::solver_1d {

template<class T, class Vector, class Right_Part>
void integrate_right_part(Vector& f, const mesh::mesh_1d<T>& mesh, const Right_Part& right_part) {
    const auto integrate_function_on_element = [&mesh, &func = right_part](const size_t e, const size_t i) {
        T integral = T{ 0 };
        const auto& el = mesh.element();
        for (const size_t q : std::ranges::iota_view{ size_t{ 0 }, el.qnodes_count() })
            integral += el.weight(q) * el.qN(i, q) * func(mesh.qnode_coord(e, q));
        return integral * mesh.jacobian(mesh.segment_number(e));
    };
#pragma omp parallel for default(none) shared(f, mesh, integrate_function_on_element)
    for (size_t node = 0; node < mesh.nodes_count(); ++node)
        for (const auto node_data : mesh.node_elements(node).to_array())
            if (node_data)
                f[node] += integrate_function_on_element(node_data.element, node_data.node);
}

template<class T, physics_t Physics>
Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_right_part_1d(
    const mesh::mesh_1d<T>& mesh, const boundaries_conditions_1d<T, Physics>& boundaries,
    const std::optional<std::function<T(const T)>>& source, const bool is_neumann = false) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> right_part = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh.nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(right_part, boundaries, is_neumann);
    if (source)
        integrate_right_part(right_part, mesh, *source);
    if (is_neumann && std::abs(std::reduce(right_part.begin(), right_part.end(), T{ 0 })) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{ "It's unsolvable Neumann problem." };
    return right_part;
}

} // namespace nonlocal::solver_1d
