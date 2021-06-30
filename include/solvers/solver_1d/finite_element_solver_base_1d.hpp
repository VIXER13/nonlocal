#ifndef FINITE_ELEMENT_SOLVER_BASE_1D_HPP
#define FINITE_ELEMENT_SOLVER_BASE_1D_HPP

#include "mesh.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace nonlocal {

enum class boundary_condition_t : uint8_t {
    FIRST_KIND,
    SECOND_KIND
};

template<class T>
struct equation_parameters final {
    T lambda = T{1},
      rho    = T{1},
      c      = T{1},
      p1     = T{1},
      r      = T{0};
};

template<class T, class I>
class finite_element_solver_base_1d {
    std::shared_ptr<mesh::mesh_1d<T, I>> _mesh;

    void boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                       const std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound) const;

    void boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                        const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond) const;

    template<class Function>
    T integrate_function(const size_t e, const size_t i, const Function& func) const;

    template<class Right_Part>
    void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const Right_Part& right_part) const;

public:
    explicit finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh);

    const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh() const;

    template<class Right_Part, class Influence_Function>
    std::vector<T> stationary(const equation_parameters<T>& parameters,
                              const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                              const Right_Part& right_part,
                              const Influence_Function& influence_function) const;
};

template<class T, class I>
finite_element_solver_base_1d<T, I>::finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T, I>>& finite_element_solver_base_1d<T, I>::mesh() const { return _mesh; }

template<class T, class I>
void finite_element_solver_base_1d<T, I>::boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                        const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                                        const std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound) const {
    std::array<T*, 2> fval = {&f[0], &f[mesh()->nodes_count()-1]};
    for(size_t b = 0; b < 2; ++b)
        if(bound_cond[b].first == boundary_condition_t::FIRST_KIND) {
            for(const auto& [i, val] : K_bound[b]) {
                std::cout << i << ' ' << val << std::endl;
                f[i] += val * bound_cond[b].second;
            }
            *fval[b] = bound_cond[b].second;
        }
}

template<class T, class I>
void finite_element_solver_base_1d<T, I>::boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                         const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond) const {
    std::array<T*, 2> fval = {&f[0], &f[mesh()->nodes_count()-1]};
    for(size_t b = 0; b < 2; ++b)
        if(bound_cond[b].first == boundary_condition_t::SECOND_KIND)
            *fval[b] += bound_cond[b].second;
}

template<class T, class I>
template<class Function>
T finite_element_solver_base_1d<T, I>::integrate_function(const size_t e, const size_t i, const Function& func) const {
    T integral = 0;
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q) * func(mesh()->quad_coord(e, q));
    return integral * mesh()->jacobian();
}

template<class T, class I>
template<class Right_Part>
void finite_element_solver_base_1d<T, I>::integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const Right_Part& right_part) const {
#pragma omp parallel for default(none) shared(f, right_part)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        for(const auto& [e, i] : mesh()->node_elements(node).arr)
            if (e != std::numeric_limits<size_t>::max())
                f[node] += integrate_function(e, i, right_part);
}

template<class T, class I>
template<class Right_Part, class Influence_Function>
std::vector<T> finite_element_solver_base_1d<T, I>::stationary(const equation_parameters<T>& parameters,
                                                               const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                               const Right_Part& right_part,
                                                               const Influence_Function& influence_function) const {
    const bool neumann_task = bound_cond[0].first == boundary_condition_t::SECOND_KIND &&
                              bound_cond[1].first == boundary_condition_t::SECOND_KIND;
    if (neumann_task && bound_cond[0].second + bound_cond[1].second > 1e-5)
        throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};
    const bool nonlocal_task = parameters.p1 < 0.999;

    const size_t size = mesh()->nodes_count() + neumann_task;
    Eigen::SparseMatrix<T, Eigen::RowMajor> K_inner(size, size);
    std::array<std::vector<std::pair<size_t, T>>, 2> K_bound;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(size);



    boundary_condition_second_kind(f, bound_cond);
    integrate_right_part(f, right_part);
    boundary_condition_first_kind(f, bound_cond, K_bound);

    std::cout << f.transpose() << std::endl;
    T sum = 0;
    for(size_t i = 0; i < f.size(); ++i)
        sum += f[i];
    std::cout << sum << std::endl;

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K_inner};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh()->nodes_count())};
}

}

#endif