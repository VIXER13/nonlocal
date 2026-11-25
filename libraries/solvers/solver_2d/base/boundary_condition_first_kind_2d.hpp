#pragma once

#include "boundary_conditions_2d.hpp"
#include "problem_settings.hpp"
#include "solvers_utils.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>
#include <metamath/linear/linear.hpp>
#include <solvers/base/degree_of_freedom.hpp>

namespace nonlocal::solver_2d {

template<class T>
void check_matrices_sizes(const metamath::linear::sparse_matrix<T>& matrix,
                          const metamath::linear::sparse_matrix<T>& boundary_matrix) {
    if (matrix.rows() != boundary_matrix.rows())
        throw std::invalid_argument(
            "The number of rows in the main matrix does not match the number of rows in the boundary matrix.");
    if (matrix.cols() != boundary_matrix.cols())
        throw std::invalid_argument(
            "The number of columns in the main matrix does not match the number of columns in the boundary matrix.");
}

template<class T>
void check_matrix_and_inner_nodes_sizes(const metamath::linear::sparse_matrix<T>& matrix,
                                        const std::vector<bool>& is_inner_nodes) {
    if (is_inner_nodes.size() != DoF<T> * matrix.cols())
        throw std::invalid_argument("The size of the vector of inner nodes does not match the number of degrees of freedom.");
}

template<class T, std::integral I, std::integral J>
metamath::linear::sparse_matrix_portrait<I, J> get_first_kind_matrix_portrait(
    const metamath::linear::sparse_matrix<T, I, J>& matrix, const std::vector<bool>& is_inner_nodes, const bool is_symmetric) {
    check_matrix_and_inner_nodes_sizes(matrix, is_inner_nodes);
    const auto need_to_add = [&is_inner_nodes](const size_t row, const size_t col) {
        for (const size_t row_global : std::ranges::iota_view{ DoF<T> * row, DoF<T> * (row + 1) })
            for (const size_t col_global : std::ranges::iota_view{ DoF<T> * col, DoF<T> * (col + 1) })
                if (!is_inner_nodes[col_global] && is_inner_nodes[row_global])
                    return true;
        return false;
    };

    metamath::linear::sparse_matrix_portrait<I, J> portrait(matrix.rows(), matrix.cols());

    for (const size_t row : std::ranges::iota_view{ 0zu, matrix.rows() })
        for (const size_t col : matrix.portrait.indices_range(row)) {
            if (need_to_add(row, col))
                ++portrait.shifts[row + 1];
            if (is_symmetric && need_to_add(col, row))
                ++portrait.shifts[col + 1];
        }
    portrait.accumulate_shifts();

    portrait.allocate_indices();
    auto current_shifts = portrait.shifts;
    for (const size_t row : std::ranges::iota_view(0zu, matrix.rows()))
        for (const size_t col : matrix.portrait.indices_range(row)) {
            if (need_to_add(row, col))
                portrait.indices[current_shifts[row]++] = col;
            if (is_symmetric && need_to_add(col, row))
                portrait.indices[current_shifts[col]++] = row;
        }
    current_shifts.clear();
    current_shifts.shrink_to_fit();
    if (is_symmetric) // for non-symmetric matrices it doesn't needed, because matrix is already sorted
        portrait.sort_indices();

    return portrait;
}

template<class T>
void get_first_kind_matrix_values(metamath::linear::sparse_matrix<T>& boundary_matrix,
                                  const metamath::linear::sparse_matrix<T>& matrix, const std::vector<bool>& is_inner_nodes,
                                  const bool is_symmetric) {
    check_matrices_sizes(matrix, boundary_matrix);
    check_matrix_and_inner_nodes_sizes(matrix, is_inner_nodes);

    for (const size_t row : std::ranges::iota_view(0zu, boundary_matrix.rows())) {
        for (const size_t shift : matrix.portrait.shifts_range(row)) {
            const size_t col = matrix.portrait.indices[shift];
            for (const size_t row_loc : std::ranges::iota_view{ 0zu, DoF<T> }) {
                const size_t row_global = DoF<T> * row + row_loc;
                for (const size_t col_loc : std::ranges::iota_view{ 0zu, DoF<T> }) {
                    const size_t col_global = DoF<T> * col + col_loc;
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
}

template<class T>
metamath::linear::sparse_matrix<T> get_first_kind_matrix(const metamath::linear::sparse_matrix<T>& matrix,
                                                         const std::vector<bool>& is_inner_nodes, const bool is_symmetric) {
    metamath::linear::sparse_matrix<T> boundary_matrix{ get_first_kind_matrix_portrait(matrix, is_inner_nodes, is_symmetric) };
    get_first_kind_matrix_values(boundary_matrix, matrix, is_inner_nodes, is_symmetric);
    return boundary_matrix;
}

template<class T>
void remove_first_kind_elements(metamath::linear::sparse_matrix<T>& matrix, const std::vector<bool>& is_inner_nodes,
                                const bool set_diagonal = true) {
    check_matrix_and_inner_nodes_sizes(matrix, is_inner_nodes);
    for (const size_t row : std::ranges::iota_view(0zu, matrix.rows()))
        for (const size_t shift : matrix.portrait.shifts_range(row)) {
            using U = metamath::types::container_type_t<metamath::types::container_type_t<T>>;
            const size_t col = matrix.portrait.indices[shift];
            for (const size_t row_loc : std::ranges::iota_view{ 0zu, DoF<T> })
                for (const size_t col_loc : std::ranges::iota_view{ 0zu, DoF<T> })
                    if (!is_inner_nodes[DoF<T> * row + row_loc] || !is_inner_nodes[DoF<T> * col + col_loc]) {
                        if constexpr (DoF<T> == 1)
                            matrix.values[shift] = set_diagonal && row == col ? T{ 1 } : T{ 0 };
                        else
                            matrix.values[shift][row_loc][col_loc] = set_diagonal && row == col && row_loc == col_loc ? U{ 1 } : U{ 0 }; ///
                    }
        }
}

template<class T>
void revert_first_kind_elements(metamath::linear::sparse_matrix<T>& matrix,
                                const metamath::linear::sparse_matrix<T>& boundary_matrix, const bool is_symmetric) {
    static_assert(DoF<T> == 1, "Reverting first kind elements is only supported for scalar problems.");
    check_matrices_sizes(matrix, boundary_matrix);
    for (const size_t row : std::ranges::iota_view{ 0zu, boundary_matrix.rows() }) {
        for (const size_t shift : boundary_matrix.portrait.shifts_range(row)) {
            if (const size_t col = boundary_matrix.portrait.indices[shift]; !is_symmetric || row <= col)
                matrix(row, col) = boundary_matrix.values[shift];
            else
                matrix(col, row) = boundary_matrix.values[shift];
        }
    }
}

template<std::floating_point T, physics_t Physics, size_t DoF>
auto calc_first_kind_vector(const mesh::mesh_container_2d<T>& mesh,
                            const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    using entity_t = std::conditional_t<DoF == 1, T, std::array<T, DoF>>;
    std::vector<entity_t> first_kind_vector(mesh.nodes_count(), entity_t{});
    // clang-format off
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&first_kind_vector, &mesh](const first_kind_2d<T, Physics>& condition, const size_t, const size_t node, const size_t degree) {
            if constexpr (DoF == 1) {
                if (auto& val = first_kind_vector[node]; val == T{ 0 })
                    val = condition(mesh.node_coord(node));
            } else {
                if (auto& val = first_kind_vector[node][degree]; val == T{ 0 })
                    val = condition(mesh.node_coord(node));
            }
        });
    // clang-format on
    return first_kind_vector;
}

template<class T, std::floating_point U, physics_t Physics, size_t DoF>
void first_kind_fill_2d(std::vector<T>& right_part, const mesh::mesh_container_2d<U>& mesh,
                        const boundaries_conditions_2d<U, Physics, DoF>& boundaries_conditions, const bool include_value = true) {
    static_assert(std::is_same_v<U, metamath::types::container_type_t<T>>,
                  "The floating point type of the mesh and the vector shall be the same");
    // clang-format off
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&mesh, &right_part, include_value](const first_kind_2d<U, Physics>& condition, const size_t, const size_t row, const size_t degree) {
            if constexpr (DoF == 1)
                right_part[row] = include_value ? condition(mesh.node_coord(row)) : U{ 0 };
            else
                right_part[row][degree] = include_value ? condition(mesh.node_coord(row)) : U{ 0 };
        });
    // clang-format on
}

template<class T, class U, std::floating_point V, physics_t Physics, size_t DoF>
void boundary_condition_first_kind_2d(metamath::linear::sparse_matrix<T>& matrix, std::vector<U>& right_part,
                                      const problem_settings& settings, const mesh::mesh_container_2d<V>& mesh,
                                      const boundaries_conditions_2d<V, Physics, DoF>& boundaries_conditions) {
    static_assert(std::is_same_v<U, metamath::types::container_type_t<T>> &&
                  std::is_same_v<V, metamath::types::container_type_t<U>>, ///
                  "The floating point type of the mesh, matrix and the vector shall be the same");
    const auto boundary_matrix = get_first_kind_matrix(matrix, settings.is_inner_nodes, settings.is_symmetric());
    const auto boundary_vector = calc_first_kind_vector(mesh, boundaries_conditions);
    remove_first_kind_elements(matrix, settings.is_inner_nodes);
    using namespace metamath::operators;
    right_part -= boundary_matrix * boundary_vector;
    first_kind_fill_2d(right_part, mesh, boundaries_conditions);
}

} // namespace nonlocal::solver_2d