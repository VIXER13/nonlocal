#pragma once

#include "conductivity_matrix_2d.hpp"
#include "heat_capacity_matrix_2d.hpp"
#include "convection_condition_2d.hpp"
#include "radiation_condition_2d.hpp"
#include "init_problem_settings.hpp"
#include "thermal_parameters_2d.hpp"

#include <solvers/slae/init_solver.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
class nonstationary_heat_equation_solver_2d final {
    std::shared_ptr<mesh::mesh_2d<T>> _mesh;
    std::unique_ptr<slae::iterative_solver_base<T>> _slae_solver;
    conductivity_matrix_2d<T> _conductivity;
    metamath::linear::sparse_matrix<T> _boundary_matrix;
    heat_capacity_matrix_2d<T> _capacity;
    std::vector<T> _right_part;
    std::vector<T> _temperature_curr;
    std::vector<T> _temperature_next;

    std::function<T(const std::array<T, 2>&)> _inner_flux;
    thermal_boundaries_conditions_2d<T> _boundaries_conditions;
    evaluated_conductivity_2d<T> _parameters;
    T _time_step = T{1};
    T _time = T{0};
    bool is_symmetric = false;

public:
    explicit nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);

    const std::vector<T>& temperature() const noexcept;
    heat_equation_solution_2d<T> solution(const bool flux = true) const;
    constexpr T time_step() const noexcept;
    constexpr T time() const noexcept;

    void compute(const parameters_2d<T>& parameters,
                 thermal_boundaries_conditions_2d<T>&& boundaries_conditions,
                 const T time_step, 
                 const std::function<T(const std::array<T, 2>&)>& right_part = nullptr,
                 const std::function<T(const std::array<T, 2>&)>& init_dist = nullptr, 
                 const T initial_time = T{0});

    void calc_step();
};

template<std::floating_point T>
nonstationary_heat_equation_solver_2d<T>::nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _mesh{mesh}
    , _conductivity{*_mesh}
    , _capacity{*_mesh}
    , _right_part(mesh->container().nodes_count(), T{0})
    , _temperature_curr(mesh->container().nodes_count(), T{0})
    , _temperature_next(mesh->container().nodes_count(), T{0}) {}

template<std::floating_point T>
const std::vector<T>& nonstationary_heat_equation_solver_2d<T>::temperature() const noexcept {
    return _temperature_next;
}

template<std::floating_point T>
heat_equation_solution_2d<T> nonstationary_heat_equation_solver_2d<T>::solution(const bool flux) const {
    heat_equation_solution_2d<T> sol{_mesh, _parameters, _temperature_next};
    if (flux)
        sol.calc_flux();
    return sol;
}

template<std::floating_point T>
constexpr T nonstationary_heat_equation_solver_2d<T>::time_step() const noexcept {
    return _time_step;
}

template<std::floating_point T>
constexpr T nonstationary_heat_equation_solver_2d<T>::time() const noexcept {
    return _time;
}

template<std::floating_point T>
void nonstationary_heat_equation_solver_2d<T>::compute(const parameters_2d<T>& parameters,
                                                       thermal_boundaries_conditions_2d<T>&& boundaries_conditions,
                                                       const T time_step, 
                                                       const std::function<T(const std::array<T, 2>&)>& right_part,
                                                       const std::function<T(const std::array<T, 2>&)>& init_dist, 
                                                       const T initial_time) {
    _time_step = time_step;
    _time = initial_time;
    _boundaries_conditions = std::move(boundaries_conditions);

    if (init_dist)
        for(const size_t node : _mesh->container().nodes())
            _temperature_next[node] = init_dist(_mesh->container().node_coord(node));
    _parameters = evaluate_conductivity(*_mesh, parameters, mesh::utils::nodes_to_qnodes<T>(*_mesh, _temperature_next));

    static constexpr bool Is_Stationary = false;
    auto settings = init_problem_settings(_mesh->container(), parameters, boundaries_conditions, Is_Stationary);
    log_problem_settings(settings);
    _conductivity.compute(_parameters, settings);
    convection_condition_2d(_conductivity.matrix(), settings, *_mesh, _boundaries_conditions);
    _conductivity.matrix() *= time_step;
    _boundary_matrix = get_first_kind_matrix(_conductivity.matrix(), settings.is_inner_nodes, settings.is_symmetric());
    remove_first_kind_elements(_conductivity.matrix(), settings.is_inner_nodes);

    _slae_solver = slae::init_iterative_solver(_conductivity.matrix(), settings.is_symmetric());
    
    const auto theroires_setter = std::views::all(_mesh->container().groups_2d()) |
                                  std::views::transform([](const std::string& group) { return std::pair{group, theory_t::LOCAL}; });
    settings.theories = std::unordered_map<std::string, theory_t>(theroires_setter.begin(), theroires_setter.end());
    _capacity.compute(parameters, settings);
    static constexpr bool Set_Diagonal = false;
    remove_first_kind_elements(_capacity.matrix(), settings.is_inner_nodes, Set_Diagonal);

    _conductivity.matrix() += _capacity.matrix();
}

template<std::floating_point T>
void nonstationary_heat_equation_solver_2d<T>::calc_step() {
    std::fill(_right_part.begin(), _right_part.end(), T{0});
    _temperature_curr.swap(_temperature_next);
    // TODO: radiation condition
    boundary_condition_second_kind_2d(_right_part, _conductivity.mesh(), _boundaries_conditions);
    if (_inner_flux)
        integrate_right_part(_right_part, _conductivity.mesh(), _inner_flux);
    using namespace metamath::operators;
    _right_part *= time_step();
    _right_part += _capacity.matrix().template self_adjoint<metamath::linear::matrix_part::Upper>() * _temperature_curr;
    _right_part -= _boundary_matrix * calc_first_kind_vector(_mesh->container(), _boundaries_conditions);
    _temperature_next = _slae_solver->solve(_right_part, _temperature_curr);
}

}