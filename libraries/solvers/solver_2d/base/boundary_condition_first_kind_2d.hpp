#pragma once

#include "boundary_conditions_2d.hpp"
#include "solvers_utils.hpp"

#include <metamath/linear/linear.hpp>
#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d {

class _boundary_condition_first_kind_2d final {
    explicit constexpr _boundary_condition_first_kind_2d() noexcept = default;

    template<class T, std::integral I, physics_t Physics, size_t DoF>
    static std::vector<T> calc_vector(const mesh::mesh_container_2d<T, I>& mesh,
                                      const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
        std::vector<T> x(DoF * mesh.nodes_count(), T{0});
        utils::run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
            [&x, &mesh](const first_kind_2d<T, Physics>& condition, const size_t, const size_t node, const size_t degree) {
                if (T& val = x[DoF * node + degree]; val == T{0})
                    val = condition(mesh.node_coord(node));
            });
        return x;
    }

    template<class T, std::integral I, physics_t Physics, size_t DoF>
    static void set_values(std::vector<T>& f, const std::vector<T>& x, const mesh::mesh_2d<T, I>& mesh,
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
    template<class T, std::integral I, physics_t Physics, size_t DoF>
    friend void boundary_condition_first_kind_2d(std::vector<T>& f, const mesh::mesh_2d<T, I>& mesh,
                                                 const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions,
                                                 const metamath::linear::sparse_matrix<T>& K_bound);
};

template<class T, std::integral I, physics_t Physics, size_t DoF>
void boundary_condition_first_kind_2d(std::vector<T>& f,
                                      const mesh::mesh_2d<T, I>& mesh,
                                      const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions,
                                      const metamath::linear::sparse_matrix<T>& K_bound) {
    const auto x = _boundary_condition_first_kind_2d::calc_vector(mesh.container(), boundaries_conditions);
    const auto process_nodes = mesh.process_nodes();
    const auto result = K_bound * x;
    for (const size_t row : std::ranges::iota_view(DoF * process_nodes.front(), DoF * process_nodes.back() + 1))
        f[row] -= result[row - DoF * process_nodes.front()];
    _boundary_condition_first_kind_2d::set_values(f, x, mesh, boundaries_conditions);
}

template<class T, std::integral I, physics_t Physics, size_t DoF>
void first_kind_matrix_fill_2d(metamath::linear::sparse_matrix<T>& K, std::vector<T>& residual, const mesh::mesh_2d<T, I>& mesh,
                               const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    utils::run_by_boundaries<first_kind_2d, Physics>(mesh.container(), boundaries_conditions,
        [&K, &mesh, &residual, process_nodes = mesh.process_nodes()]
        (const  first_kind_2d<T, Physics>&, const size_t be, const size_t row, const size_t) {
            if (row >= process_nodes.front() && row <= process_nodes.back()) {
                for (const size_t s : K.portrait.shifts_range(row))
                    if (K.portrait.indices[s] >= row)
                        K.values[s] = (K.portrait.indices[s] == row ? T(1) : T(0));
                residual[row] = T(0);
            }
        });
}
    
    
}