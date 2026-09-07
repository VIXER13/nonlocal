#pragma once

#include "boundary_conditions_2d.hpp"
#include "problem_settings.hpp"
#include "solvers_utils.hpp"

#include <metamath/linear/linear.hpp>
#include <mesh/mesh_2d/mesh_2d.hpp>
#include <solvers/base/degree_of_freedom.hpp>

namespace nonlocal::solver_2d {

inline bool has_first_kind(const std::vector<bool>& is_inner_nodes, const size_t node, const size_t DoF) {
    for (const size_t degree : std::ranges::iota_view{0zu, DoF})
        if (!is_inner_nodes[DoF * node + degree])
            return true;
    return false;
}

template<class T>
metamath::linear::sparse_matrix<T> get_first_kind_matrix(const metamath::linear::sparse_matrix<T>& matrix, 
                                                         const std::vector<bool>& is_inner_nodes, const bool is_symmetric) {
    if (is_inner_nodes.size() != DoF<T> * matrix.cols())
        throw std::invalid_argument("The size of the vector of inner nodes does not match the number of degrees of freedom.");
    metamath::linear::sparse_matrix<T> boundary_matrix(matrix.cols(), matrix.cols());

    const auto need_to_add = [&is_inner_nodes](const size_t row, const size_t col) {
        for (const size_t row_global : std::ranges::iota_view{DoF<T> * row, DoF<T> * (row + 1)})
            for (const size_t col_global : std::ranges::iota_view{DoF<T> * col, DoF<T> * (col + 1)})
                if (row_global != col_global && !is_inner_nodes[col_global] && is_inner_nodes[row_global])
                    return true;
        return false;
    };

    for(const size_t row : std::ranges::iota_view{0zu, boundary_matrix.rows()})
        for(const size_t col : matrix.portrait.indices_range(row)) {
            if (need_to_add(row, col))
                ++boundary_matrix.portrait.shifts[row + 1];
            if (is_symmetric && need_to_add(col, row))
                ++boundary_matrix.portrait.shifts[col + 1];
        }
    boundary_matrix.portrait.accumulate_shifts();

    boundary_matrix.portrait.allocate_indices();
    auto current_shifts = boundary_matrix.portrait.shifts;
    for(const size_t row : std::ranges::iota_view(0zu, boundary_matrix.rows()))
        for(const size_t col : matrix.portrait.indices_range(row)) {
            if (need_to_add(row, col))
                boundary_matrix.portrait.indices[current_shifts[row]++] = col;
            if (is_symmetric && need_to_add(col, row))
                boundary_matrix.portrait.indices[current_shifts[col]++] = row;
        }
    current_shifts.clear();
    current_shifts.shrink_to_fit();
    if (is_symmetric) // for non-symmetric matrices it doesn't needed, because matrix is already sorted
        boundary_matrix.portrait.sort_indices();

    boundary_matrix.allocate_values();
    for(const size_t row : std::ranges::iota_view(0zu, boundary_matrix.rows()))
        for(const size_t shift : matrix.portrait.shifts_range(row)) {
            const size_t col = matrix.portrait.indices[shift];
            for(const size_t row_loc : std::ranges::iota_view{0zu, DoF<T>}) {
                const size_t row_global = DoF<T> * row + row_loc;
                for(const size_t col_loc : std::ranges::iota_view{0zu, DoF<T>}) {
                    const size_t col_global = DoF<T> * col + col_loc;
                    if (row_global != col_global) {
                        if (!is_inner_nodes[col_global] && is_inner_nodes[row_global]) {
                            if constexpr (DoF<T> == 1)
                                boundary_matrix(row, col) = matrix.values[shift];
                            else
                                boundary_matrix(row, col)[row_loc][col_loc] = matrix.values[shift][row_loc][col_loc];
                        }
                        if (is_symmetric && !is_inner_nodes[row_global] && is_inner_nodes[col_global]) {
                            if constexpr (DoF<T> == 1)
                                boundary_matrix(col, row) = matrix.values[shift];
                            else
                                boundary_matrix(col, row)[col_loc][row_loc] = matrix.values[shift][row_loc][col_loc];
                        }
                    }
                }
            }
        }

    return boundary_matrix;
}

template<class T>
void remove_first_kind_elements(metamath::linear::sparse_matrix<T>& matrix,
                                const std::vector<bool>& is_inner_nodes, const bool set_diagonal = true) {
    if (is_inner_nodes.size() != DoF<T> * matrix.cols())
        throw std::invalid_argument("The size of the vector of inner nodes does not match the number of degrees of freedom.");
    for(const size_t row : std::ranges::iota_view(0zu, matrix.rows()))
        for (const size_t shift : matrix.portrait.shifts_range(row)) {
            using U = metamath::types::container_type_t<metamath::types::container_type_t<T>>;
            const size_t col = matrix.portrait.indices[shift];
            for(const size_t row_loc : std::ranges::iota_view{0zu, DoF<T>})
                for(const size_t col_loc : std::ranges::iota_view{0zu, DoF<T>})
                    if (!is_inner_nodes[DoF<T> * row + row_loc] || !is_inner_nodes[DoF<T> * col + col_loc]) {
                        if constexpr (DoF<T> == 1)
                            matrix.values[shift] = set_diagonal && row == col ? T{1} : T{0};
                        else
                            matrix.values[shift][row_loc][col_loc] = set_diagonal && row == col && row_loc == col_loc ? U{1} : U{0};
                    }
        }
}

template<std::floating_point T, physics_t Physics, size_t DoF>
auto calc_first_kind_vector(const mesh::mesh_container_2d<T>& mesh,
                            const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    using entity_t = std::conditional_t<DoF == 1, T, std::array<T, DoF>>;
    std::vector<entity_t> first_kind_vector(mesh.nodes_count(), entity_t{});
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&first_kind_vector, &mesh](const first_kind_2d<T, Physics>& condition, const size_t, const size_t node, const size_t degree) {
            if constexpr (DoF == 1) {
                if (auto& val = first_kind_vector[node]; val == T{0})
                    val = condition(mesh.node_coord(node));
            } else {
                if (auto& val = first_kind_vector[node][degree]; val == T{0})
                    val = condition(mesh.node_coord(node));
            }
        });
    return first_kind_vector;
}

template<class T, std::floating_point U, physics_t Physics, size_t DoF>
void first_kind_fill_2d(std::vector<T>& right_part, const mesh::mesh_container_2d<U>& mesh,
                        const boundaries_conditions_2d<U, Physics, DoF>& boundaries_conditions,
                        const bool include_value = true) {
    static_assert(std::is_same_v<U, metamath::types::container_type_t<T>>, "The floating point type of the mesh and the vector shall be the same");
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&mesh, &right_part, include_value](const  first_kind_2d<U, Physics>& condition, const size_t, const size_t row, const size_t degree) {
            if constexpr (DoF == 1)
                right_part[row] = include_value ? condition(mesh.node_coord(row)) : U{0};
            else
                right_part[row][degree] = include_value ? condition(mesh.node_coord(row)) : U{0};
        });
}

template<class T, class U, std::floating_point V, physics_t Physics, size_t DoF>
void boundary_condition_first_kind_2d(metamath::linear::sparse_matrix<T>& matrix, std::vector<U>& right_part,
                                      const problem_settings& settings, const mesh::mesh_container_2d<V>& mesh,
                                      const boundaries_conditions_2d<V, Physics, DoF>& boundaries_conditions) {
    static_assert(std::is_same_v<U, metamath::types::container_type_t<T>> &&
                  std::is_same_v<V, metamath::types::container_type_t<U>>,
                  "The floating point type of the mesh, matrix and the vector shall be the same");
    const auto boundary_matrix = get_first_kind_matrix(matrix, settings.is_inner_nodes, settings.is_symmetric());
    const auto boundary_vector = calc_first_kind_vector(mesh, boundaries_conditions);
    remove_first_kind_elements(matrix, settings.is_inner_nodes);
    using namespace metamath::operators;
    right_part -= boundary_matrix * boundary_vector;
    first_kind_fill_2d(right_part, mesh, boundaries_conditions);
}
    
}