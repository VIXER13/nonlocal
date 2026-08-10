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
    std::vector<T> _temperature;

    std::function<T(const std::array<T, 2>&)> _inner_flux;
    thermal_boundaries_conditions_2d<T> _boundaries_conditions;
    evaluated_thermal_parameters<T> _parameters;
    T _time_step = T{1};
    T _time = T{0};
    bool is_symmetric = false;

public:
    explicit nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);

    const std::vector<T>& temperature() const noexcept;
    heat_equation_solution_2d<T> solution(const bool flux = true) const;
    constexpr T time_step() const noexcept;
    constexpr T time() const noexcept;

    void compute(const raw_thermal_parameters<T>& parameters,
                 thermal_boundaries_conditions_2d<T>&& boundaries_conditions,
                 const T time_step, 
                 const std::function<T(const std::array<T, 2>&)>& inner_flux = nullptr,
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
    , _temperature(mesh->container().nodes_count(), T{0}) {}

template<std::floating_point T>
const std::vector<T>& nonstationary_heat_equation_solver_2d<T>::temperature() const noexcept {
    return _temperature;
}

template<std::floating_point T>
heat_equation_solution_2d<T> nonstationary_heat_equation_solver_2d<T>::solution(const bool flux) const {
    heat_equation_solution_2d<T> sol{_mesh, _parameters, temperature()};
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
void nonstationary_heat_equation_solver_2d<T>::compute(const raw_thermal_parameters<T>& parameters,
                                                       thermal_boundaries_conditions_2d<T>&& boundaries_conditions,
                                                       const T time_step, 
                                                       const std::function<T(const std::array<T, 2>&)>& inner_flux,
                                                       const std::function<T(const std::array<T, 2>&)>& init_dist, 
                                                       const T initial_time) {
    _time_step = time_step;
    _time = initial_time;
    _inner_flux = inner_flux;
    _boundaries_conditions = std::move(boundaries_conditions);

    if (init_dist)
        for(const size_t node : _mesh->container().nodes())
            _temperature[node] = init_dist(_mesh->container().node_coord(node));
    _parameters = evaluate_conductivity(*_mesh, parameters, mesh::utils::nodes_to_qnodes<T>(*_mesh, _temperature));

    static constexpr bool Is_Stationary = false;
    auto settings = init_problem_settings(_mesh->container(), parameters, boundaries_conditions, Is_Stationary);
    const bool is_symmetric = settings.is_symmetric();
    log_problem_settings(settings);
    _conductivity.compute(_parameters, settings);
    convection_condition_2d(_conductivity.matrix(), settings, *_mesh, _boundaries_conditions);
    _boundary_matrix = get_first_kind_matrix(_conductivity.matrix(), settings.is_inner_nodes, is_symmetric);
    remove_first_kind_elements(_conductivity.matrix(), settings.is_inner_nodes);

    _slae_solver = slae::init_iterative_solver(_conductivity.matrix(), is_symmetric);
    
    settings.set_fully_local();
    _capacity.compute(_parameters, settings);
    _capacity.matrix() /= time_step;
    static constexpr bool Set_Diagonal = false;
    remove_first_kind_elements(_capacity.matrix(), settings.is_inner_nodes, Set_Diagonal);

    if (is_symmetric)
        _conductivity.matrix() += _capacity.matrix();
    else
        _conductivity.matrix() += _capacity.matrix().template self_adjoint<metamath::linear::matrix_part::Upper>();
}

template<std::floating_point T>
void nonstationary_heat_equation_solver_2d<T>::calc_step() {
    std::fill(_right_part.begin(), _right_part.end(), T{0});
    // TODO: radiation condition
    boundary_condition_second_kind_2d(_right_part, _conductivity.mesh(), _boundaries_conditions);
    if (_inner_flux)
        integrate_right_part(_right_part, _conductivity.mesh(), _inner_flux);
    using namespace metamath::operators;
    _right_part += _capacity.matrix().template self_adjoint<metamath::linear::matrix_part::Upper>() * temperature();
    _right_part -= _boundary_matrix * calc_first_kind_vector(_mesh->container(), _boundaries_conditions);
    first_kind_fill_2d(_right_part, _mesh->container(), _boundaries_conditions);
    _temperature = _slae_solver->solve(_right_part, temperature());
    _time += time_step();
}

}