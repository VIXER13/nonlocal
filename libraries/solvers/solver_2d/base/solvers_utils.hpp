#ifndef NONLOCAL_SOLVERS_2D_UTILS_HPP
#define NONLOCAL_SOLVERS_2D_UTILS_HPP

#include "boundary_conditions_2d.hpp"

#include "mesh_2d.hpp"

namespace nonlocal::utils {

template<template<class, auto...> class Condition, auto... Args, class T, class I, physics_t Physics, class Callback>
void run_by_boundary(const mesh::mesh_container_2d<T, I>& mesh, const std::string& bound_name,
                     const boundary_condition_2d<T, Physics> *const condition, const size_t degree, 
                     const Callback& callback) {
    if (const auto *const cond = dynamic_cast<const Condition<T, Args...>*>(condition))
        for(const size_t be : mesh.elements(bound_name))
            for(const size_t node : mesh.nodes(be))
                callback(*cond, be, node, degree);
}

template<template<class, auto...> class Condition, auto... Args, class T, class I, physics_t Physics, size_t DoF, class Callback>
void run_by_boundaries(const mesh::mesh_container_2d<T, I>& mesh,
                       const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions,
                       const Callback& callback) {
    for(const auto& [bound_name, conditions] : boundaries_conditions)
        if constexpr (DoF == 1)
            run_by_boundary<Condition, Args...>(mesh, bound_name, conditions.get(), 0, callback);
        else
            for(const size_t degree : std::ranges::iota_view{0u, DoF})
                run_by_boundary<Condition, Args...>(mesh, bound_name, conditions[degree].get(), degree, callback);
}

template<class T, class I, physics_t Physics, size_t DoF>
std::vector<bool> inner_nodes(const mesh::mesh_container_2d<T, I>& mesh,
                              const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    std::vector<bool> is_inner(DoF * mesh.nodes_count(), true);
    run_by_boundaries<first_kind_2d, Physics>(mesh, boundaries_conditions,
        [&is_inner](const first_kind_2d<T, Physics>&, const size_t, const size_t node, const size_t degree) {
            is_inner[DoF * node + degree] = false;
        });
    return is_inner;
}

}

#endif