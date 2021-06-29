#ifndef FINITE_ELEMENT_SOLVER_BASE_1D_HPP
#define FINITE_ELEMENT_SOLVER_BASE_1D_HPP

#include "mesh.hpp"

namespace nonlocal {

enum class boundary_condition_t : uint8_t {
    FIRST_KIND,
    SECOND_KIND
};

template<class T>
struct equation_parameters final {
    T lambda = T{1},
      rho    = T{1},
      c      = T{1};
};

template<class T, class I>
class finite_element_solver_base_1d {
    std::shared_ptr<mesh::mesh_1d<T, I>> _mesh;

public:
    explicit finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh);

    const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh() const;

    template<class Right_Part, class Influence_Function>
    std::vector<T> stationary(const equation_parameters<T>& parameters,
                              const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                              const Right_Part& right_part,
                              const T p1, const T r, const Influence_Function& influence_function) const;
};

template<class T, class I>
finite_element_solver_base_1d<T, I>::finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T, I>>& finite_element_solver_base_1d<T, I>::mesh() const { return _mesh; }

template<class T, class I>
template<class Right_Part, class Influence_Function>
std::vector<T> finite_element_solver_base_1d<T, I>::stationary(const equation_parameters<T>& parameters,
                                                               const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                               const Right_Part& right_part,
                                                               const T p1, const T r, const Influence_Function& influence_function) const {
    std::vector<T> solution(mesh()->nodes_count(), 5);
    return std::move(solution);
}

}

#endif