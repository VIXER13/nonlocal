#ifndef NONLOCAL_SOLVERS_2D_UTILS_HPP
#define NONLOCAL_SOLVERS_2D_UTILS_HPP

#include "mesh_2d.hpp"
#include "boundary_conditions_2d.hpp"

namespace nonlocal::utils {

template<size_t DoF, template<class> class Condition, class T, class I, class Callback>
void run_by_boundary(const mesh::mesh_container_2d<T, I>& mesh, const std::string& bound_name,
                     const boundary_condition_2d<T>& condition, const size_t degree, 
                     const Callback& callback) {
    if (dynamic_cast<const Condition<T>*>(&condition))
        for(const size_t be : mesh.elements(bound_name))
            for(const size_t node : mesh.nodes(be))
                callback(condition, be, node, degree);
}

template<size_t DoF, template<class> class Condition, class T, class I, class Conditions_Map, class Callback>
void run_by_boundaries(const mesh::mesh_container_2d<T, I>& mesh,
                       const Conditions_Map& boundaries_conditions,
                       const Callback& callback) {
    for(const auto& [bound_name, conditions] : boundaries_conditions)
        if constexpr (DoF == 1)
            run_by_boundary<DoF, Condition>(mesh, bound_name, *conditions, 0, callback);
        else
            for(const size_t degree : std::ranges::iota_view{0u, DoF})
                run_by_boundary<DoF, Condition>(mesh, bound_name, *conditions[degree], degree, callback);
}

template<size_t DoF, class T, class I, class Conditions_Map>
std::vector<bool> inner_nodes(const mesh::mesh_container_2d<T, I>& mesh,
                              const Conditions_Map& boundaries_conditions) {
    std::vector<bool> is_inner(DoF * mesh.nodes_count(), true);
    run_by_boundaries<DoF, first_kind_2d>(mesh, boundaries_conditions,
        [&is_inner](const auto&, const size_t, const size_t node, const size_t degree) {
            is_inner[DoF * node + degree] = false;
        });
    return is_inner;
}

}

#endif