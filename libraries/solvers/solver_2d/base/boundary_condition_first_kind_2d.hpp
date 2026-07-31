#pragma once

#include "boundary_conditions_2d.hpp"
#include "problem_settings.hpp"
#include "solvers_utils.hpp"

#include <metamath/linear/linear.hpp>
#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d {

class _boundary_condition_first_kind_2d final {
    explicit constexpr _boundary_condition_first_kind_2d() noexcept = default;

    template<class T, physics_t Physics, size_t DoF>
    static std::vector<T> calc_vector(const mesh::mesh_container_2d<T>& mesh,
                                      const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
        std::vector<T> x(DoF * mesh.nodes_count(), T{0});
        utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
            [&x, &mesh](const first_kind_2d<T, Physics>& condition, const size_t, const size_t node, const size_t degree) {
                if (T& val = x[DoF * node + degree]; val == T{0})
                    val = condition(mesh.node_coord(node));
            });
        return x;
    }

    template<class T, physics_t Physics, size_t DoF>
    static void set_values(std::vector<T>& f, const std::vector<T>& x, const mesh::mesh_2d<T>& mesh,
                           const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
        utils::run_by_boundaries<first_kind_2d, Physics>(mesh.container(), boundaries_conditions,
            [&f, &x, process_nodes = mesh.process_nodes()](const first_kind_2d<T, Physics>&, const size_t, const size_t node, const size_t degree) {
                if (node >= process_nodes.front() && node <= process_nodes.back()) {
                    const size_t index = DoF * node + degree;
                    f[index] = x[index];
                }
            });
    }

public:
    template<class T, physics_t Physics, size_t DoF>
    friend void boundary_condition_first_kind_2d(std::vector<T>& f, const mesh::mesh_2d<T>& mesh,
                                                 const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions,
                                                 const metamath::linear::sparse_matrix<T>& K_bound);
};

template<class T, physics_t Physics, size_t DoF>
void boundary_condition_first_kind_2d(std::vector<T>& f,
                                      const mesh::mesh_2d<T>& mesh,
                                      const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions,
                                      const metamath::linear::sparse_matrix<T>& K_bound) {
    const auto x = _boundary_condition_first_kind_2d::calc_vector(mesh.container(), boundaries_conditions);
    const auto process_nodes = mesh.process_nodes();
    const auto result = K_bound * x;
    for (const size_t row : std::ranges::iota_view(DoF * process_nodes.front(), DoF * process_nodes.back() + 1))
        f[row] -= result[row - DoF * process_nodes.front()];
    _boundary_condition_first_kind_2d::set_values(f, x, mesh, boundaries_conditions);
}

template<class T>
metamath::linear::sparse_matrix<T> get_first_kind_matrix(const metamath::linear::sparse_matrix<T>& matrix, 
                                                         const std::vector<bool>& is_inner_nodes, const bool is_symmetric) {
    metamath::linear::sparse_matrix<T> boundary_matrix(matrix.cols(), matrix.cols());

    for(const size_t row : std::ranges::iota_view{0zu, boundary_matrix.rows()})
        for(const size_t col : matrix.portrait.indices_range(row))
            if (row != col) {
                if (!is_inner_nodes[col] && is_inner_nodes[row])
                    ++boundary_matrix.portrait.shifts[row + 1];
                if (is_symmetric && !is_inner_nodes[row] && is_inner_nodes[col])
                    ++boundary_matrix.portrait.shifts[col + 1];
            }
    boundary_matrix.portrait.accumulate_shifts();

    boundary_matrix.portrait.allocate_indices();
    auto current_shifts = boundary_matrix.portrait.shifts;
    for(const size_t row : std::ranges::iota_view(0zu, boundary_matrix.rows()))
        for(const size_t col : matrix.portrait.indices_range(row))
            if (row != col) {
                if (!is_inner_nodes[col] && is_inner_nodes[row])
                    boundary_matrix.portrait.indices[current_shifts[row]++] = col;
                if (is_symmetric && !is_inner_nodes[row] && is_inner_nodes[col])
                    boundary_matrix.portrait.indices[current_shifts[col]++] = row;
            }
    current_shifts.clear();
    current_shifts.shrink_to_fit();
    if (is_symmetric) // for non-symmetric matrices it doesn't needed, because matrix is already sorted
        boundary_matrix.portrait.sort_indices();

    boundary_matrix.allocate_values();
    for(const size_t row : std::ranges::iota_view(0zu, boundary_matrix.rows()))
        for(const size_t shift : matrix.portrait.shifts_range(row))
            if (const size_t col = matrix.portrait.indices[shift]; row != col) {
                if (!is_inner_nodes[col] && is_inner_nodes[row])
                    boundary_matrix(row, col) = matrix.values[shift];
                if (is_symmetric && !is_inner_nodes[row] && is_inner_nodes[col])
                    boundary_matrix(col, row) = matrix.values[shift];
            }

    return boundary_matrix;
}

template<class T>
void remove_first_kind_elements(metamath::linear::sparse_matrix<T>& matrix,
                                const std::vector<bool>& is_inner_nodes, const bool set_diagonal = true) {
    for(const size_t row : std::ranges::iota_view(0zu, matrix.rows()))
        for (const size_t shift : matrix.portrait.shifts_range(row))
            if (const size_t col = matrix.portrait.indices[shift]; !is_inner_nodes[row] || !is_inner_nodes[col])
                matrix.values[shift] = set_diagonal && row == col ? T{1} : T{0};
}

template<class T, physics_t Physics, size_t DoF>
std::vector<T> calc_first_kind_vector(const mesh::mesh_container_2d<T>& mesh,
                                      const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    std::vector<T> x(DoF * mesh.nodes_count(), T{0});
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&x, &mesh](const first_kind_2d<T, Physics>& condition, const size_t, const size_t node, const size_t degree) {
            if (T& val = x[DoF * node + degree]; val == T{0})
                val = condition(mesh.node_coord(node));
        });
    return x;
}

template<class T, physics_t Physics, size_t DoF>
void first_kind_fill_2d(std::vector<T>& right_part, const mesh::mesh_container_2d<T>& mesh,
                        const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions,
                        const bool include_value = true) {
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&mesh, &right_part, include_value](const  first_kind_2d<T, Physics>& condition, const size_t, const size_t row, const size_t) {
            right_part[row] = include_value ? condition(mesh.node_coord(row)) : T{0};
        });
}

template<class T, physics_t Physics, size_t DoF>
void boundary_condition_first_kind_2d(metamath::linear::sparse_matrix<T>& matrix, std::vector<T>& right_part,
                                      const problem_settings& settings, const mesh::mesh_container_2d<T>& mesh,
                                      const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    const auto boundary_matrix = get_first_kind_matrix(matrix, settings.is_inner_nodes, settings.is_symmetric());
    const auto boundary_vector = calc_first_kind_vector(mesh, boundaries_conditions);
    remove_first_kind_elements(matrix, settings.is_inner_nodes);
    using namespace metamath::operators;
    right_part -= boundary_matrix * boundary_vector;
    first_kind_fill_2d(right_part, mesh, boundaries_conditions);
}
    
}