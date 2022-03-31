#ifndef NONSTATIONARY_HEAT_EQUATION_SOLVER_2D_HPP
#define NONSTATIONARY_HEAT_EQUATION_SOLVER_2D_HPP

#include "conjugate_gradient.hpp"
#include "convection_condition_2d.hpp"
#include "heat_capacity_matrix_2d.hpp"
#include "heat_equation_solution_2d.hpp"
#include "right_part_2d.hpp"
#include "thermal_conductivity_matrix_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class nonstationary_heat_equation_solver_2d final {
    std::unique_ptr<slae::conjugate_gradient<T, Matrix_Index>> slae_solver;
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> _conductivity;
    heat_capacity_matrix_2d<T, I, Matrix_Index> _capacity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    T _time = 0;
    const T _tau = 1;

public:
    explicit nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy, const T tau, const T time = 0);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;

    template<material_t Material, class Init_Dist, class Influence_Function>
    void compute(const equation_parameters_2d<T, Material>& equation_param,
                 const std::unordered_map<std::string, std::array<boundary_condition_t, 1>>& boundary_types,
                 const Init_Dist& init_dist,
                 const T p1,
                 const Influence_Function& influence_function);

    template<class Right_Part>
    void calc_step(const std::unordered_map<std::string, T>& alpha,
                   const std::unordered_map<std::string, nonstationary_boundary_2d_t<boundary_condition_t, T, 1>>& boundary_condition,
                   const Right_Part& right_part);
};

template<class T, class I, class Matrix_Index>
nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy, const T tau, const T time)
    : _conductivity{mesh_proxy}
    , _capacity{mesh_proxy}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh_proxy->mesh().nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh_proxy->mesh().nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh_proxy->mesh().nodes_count())}
    , _time{time}
    , _tau{tau} {}

template<class T, class I, class Matrix_Index>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I, class Matrix_Index>
template<material_t Material, class Init_Dist, class Influence_Function>
void nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::compute(const equation_parameters_2d<T, Material>& equation_param,
                                                                        const std::unordered_map<std::string, std::array<boundary_condition_t, 1>>& boundary_types,
                                                                        const Init_Dist& init_dist,
                                                                        const T p1,
                                                                        const Influence_Function& influence_function) {
    const auto& mesh = _conductivity.mesh_proxy()->mesh();
    const std::vector<bool> is_inner = inner_nodes(mesh, boundary_types);
    _conductivity.template calc_matrix<Material>(equation_param.lambda, is_inner, p1, influence_function);
    convection_condition_matrix_part_2d(_conductivity.matrix_inner(), *_conductivity.mesh_proxy(), boundary_types, equation_param.alpha);
    _capacity.calc_matrix(equation_param.c, equation_param.rho, is_inner);

    _conductivity.matrix_inner() *= _tau;
    _conductivity.matrix_bound() *= _tau;
    _conductivity.matrix_inner() += _capacity.matrix_inner();
    for(const size_t i : std::views::iota(size_t{0}, is_inner.size()))
        if(!is_inner[i])
            _conductivity.matrix_inner().coeffRef(i, i) = T{1};

    for(const size_t i : std::views::iota(size_t{0}, mesh.nodes_count()))
        _temperature_prev[i] = init_dist(mesh.node(i));

    slae_solver = std::make_unique<slae::conjugate_gradient<T, Matrix_Index>>(_conductivity.matrix_inner());
}

template<class T, class I, class Matrix_Index>
template<class Right_Part>
void nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::calc_step(
    const std::unordered_map<std::string, T>& alpha,
    const std::unordered_map<std::string, nonstationary_boundary_2d_t<boundary_condition_t, T, 1>>& boundary_condition,
    const Right_Part& right_part) {
    _temperature_prev.swap(_temperature_curr);
    _time += _tau;
    _right_part.setZero();
    const auto& mesh_proxy = _conductivity.mesh_proxy();
    const auto stationary_bound = to_stationary(boundary_condition, _time);
    boundary_condition_second_kind_2d(_right_part, *mesh_proxy, stationary_bound);
    convection_condition_right_part_2d(_right_part, *mesh_proxy, stationary_bound, alpha);
    integrate_right_part<1>(_right_part, *mesh_proxy, [&right_part, t = _time](const std::array<T, 2>& x) { return right_part(t, x); });
    _right_part *= _tau;
    _right_part += _capacity.matrix_inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    boundary_condition_first_kind_2d(_right_part, *mesh_proxy, stationary_bound, _conductivity.matrix_bound());
    _temperature_curr = slae_solver->solve(_right_part, _temperature_prev);
}

}

#endif