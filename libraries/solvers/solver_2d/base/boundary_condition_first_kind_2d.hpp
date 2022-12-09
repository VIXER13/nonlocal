#ifndef NONLOCAL_BOUNDARY_CONDITION_FIRST_KIND_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_FIRST_KIND_2D_HPP

#include "boundary_conditions_2d.hpp"
#include "solvers_utils.hpp"

#include "mesh_2d.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace nonlocal
{

class _boundary_condition_first_kind_2d final {
    explicit constexpr _boundary_condition_first_kind_2d() noexcept = default;

    template<size_t DoF, class T, class I, class Conditions_Map>
    static Eigen::Matrix<T, Eigen::Dynamic, 1> calc_vector(const mesh::mesh_container_2d<T, I>& mesh,
                                                           const Conditions_Map& boundaries_conditions) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh.nodes_count());
        utils::run_by_boundaries<DoF, first_kind_2d>(mesh, boundaries_conditions,
            [&x, &mesh](const auto& condition, const size_t, const size_t node, const size_t degree) {
                if (T& val = x[DoF * node + degree]; val == T{0})
                    val = condition(mesh.node_coord(node));
            });
        return x;
    }

    template<size_t DoF, class T, class I, class Conditions_Map>
    static void set_values(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                           const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
                           const mesh::mesh_2d<T, I>& mesh,
                           const Conditions_Map& boundaries_conditions) {
        utils::run_by_boundaries<DoF, first_kind_2d>(mesh.container(), boundaries_conditions,
            [&f, &x, process_nodes = mesh.process_nodes()](const auto&, const size_t, const size_t node, const size_t degree) {
                if (node >= process_nodes.front() && node <= process_nodes.back()) {
                    const size_t index = DoF * node + degree;
                    f[index] = x[index];
                }
            });
    }

public:
    template<size_t DoF, class T, class I, class Matrix_Index, class Conditions_Map>
    friend void boundary_condition_first_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                 const mesh::mesh_2d<T, I>& mesh,
                                                 const Conditions_Map& boundaries_conditions,
                                                 const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound);
};

template<size_t DoF, class T, class I, class Matrix_Index, class Conditions_Map>
void boundary_condition_first_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                      const mesh::mesh_2d<T, I>& mesh,
                                      const Conditions_Map& boundaries_conditions,
                                      const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound) {
    const auto x = _boundary_condition_first_kind_2d::calc_vector<DoF>(mesh.container(), boundaries_conditions);
    f -= K_bound * x;
    _boundary_condition_first_kind_2d::set_values<DoF>(f, x, mesh, boundaries_conditions);
}
    
}

#endif