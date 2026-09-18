#pragma once

#include "stiffness_matrix_2d.hpp"
#include "mass_matrix.hpp"
#include "mechanical_parameters_2d.hpp"
#include "mechanical_solution_2d.hpp"

#include <solvers/slae/init_solver.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
class motion_equation_solver final {
    std::shared_ptr<mesh::mesh_2d<T>> _mesh;
    std::unique_ptr<slae::iterative_solver_base<T>> _slae_solver;
    stiffness_matrix<T> _stiffness;
    metamath::linear::sparse_matrix<metamath::linear::square_matrix<T, 2>> _boundary_matrix;
    mass_matrix<T> _mass;
    std::vector<std::array<T, 2>> _right_part;
    std::vector<std::array<T, 2>> _displacement_prev;
    std::vector<std::array<T, 2>> _displacement_curr;
    std::vector<std::array<T, 2>> _displacement_next;

    std::function<std::array<T, 2>(const std::array<T, 2>&)> _inner_pressure;
    mechanical_boundaries_conditions_2d<T> _boundaries_conditions;
    evaluated_mechanical_parameters<T> _parameters;
    T _time_step = T{1};
    T _time = T{0};

public:
    explicit motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);

    const std::vector<std::array<T, 2>>& displacement() const noexcept;
    mechanical_solution_2d<T> solution(const bool strain_and_stress = true) const;
    T time_step() const noexcept;
    T time() const noexcept;

    void compute(const raw_mechanical_parameters<T>& parameters,
                 mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                 const T time_step,
                 const std::function<std::array<T, 2>(const std::array<T, 2>&)>& inner_pressure = nullptr,
                 const std::function<std::array<T, 2>(const std::array<T, 2>&)>& init_dist = nullptr,
                 const T initial_time = T{0});

    void calc_step();
};

template<std::floating_point T>
motion_equation_solver<T>::motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _mesh{mesh}
    , _mass{mesh}
    , _stiffness{mesh} 
    , _right_part(mesh->container().nodes_count(), std::array<T, 2>{})
    , _displacement_prev(mesh->container().nodes_count(), std::array<T, 2>{})
    , _displacement_curr(mesh->container().nodes_count(), std::array<T, 2>{})
    , _displacement_next(mesh->container().nodes_count(), std::array<T, 2>{}) {}

template<std::floating_point T>
const std::vector<std::array<T, 2>>& motion_equation_solver<T>::displacement() const noexcept {
    return _displacement_next;
}

template<std::floating_point T>
mechanical_solution_2d<T> motion_equation_solver<T>::solution(const bool strain_and_stress) const {
    mechanical_solution_2d<T> sol{_mass.mesh_ptr(), _parameters, displacement()};
    if (strain_and_stress)
        sol.calc_strain_and_stress(_parameters);
    return sol;
}

template<std::floating_point T>
T motion_equation_solver<T>::time_step() const noexcept {
    return _time_step;
}

template<std::floating_point T>
T motion_equation_solver<T>::time() const noexcept {
    return _time;
}

template<std::floating_point T>
void motion_equation_solver<T>::compute(const raw_mechanical_parameters<T>& parameters,
                                        mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                                        const T time_step,
                                        const std::function<std::array<T, 2>(const std::array<T, 2>&)>& inner_pressure,
                                        const std::function<std::array<T, 2>(const std::array<T, 2>&)>& init_dist,
                                        const T initial_time) {
    _time_step = time_step;
    _time = initial_time;
    _boundaries_conditions = std::move(boundaries_conditions);

    if (init_dist) {
        for(const size_t node : _mesh->container().nodes())
            _displacement_next[node] = (*init_dist)(_mesh->container().node_coord(node));
        _displacement_curr = _displacement_next;
        _displacement_prev = _displacement_curr;
    }

    auto settings = init_problem_settings(_mesh->container(), parameters, _boundaries_conditions);
    const bool is_symmetric = settings.is_symmetric();
    log_problem_settings(settings);
    _parameters = evaluate_mechanical_parameters(_mesh, parameters);
    _stiffness.compute(_parameters, settings);
    _boundary_matrix = get_first_kind_matrix(_stiffness.matrix(), settings.is_inner_nodes, is_symmetric);
    remove_first_kind_elements(_stiffness.matrix(), settings.is_inner_nodes);

    _slae_solver = slae::init_iterative_solver(_stiffness.matrix(), is_symmetric);

    settings.set_fully_local();
    _mass.compute(_parameters, settings);
    _mass.matrix() /= time_step * time_step;
    static constexpr bool Set_Diagonal = false;
    remove_first_kind_elements(_mass.matrix(), settings.is_inner_nodes, Set_Diagonal);

    if (is_symmetric)
        _stiffness.matrix() += _mass.matrix();
    else
        _stiffness.matrix() += _mass.matrix().template self_adjoint<metamath::linear::matrix_part::Upper>();
}

template<std::floating_point T>
void motion_equation_solver<T>::calc_step() {
    std::fill(_right_part.begin(), _right_part.end(), T{0});
    _displacement_prev.swap(_displacement_curr);
    _displacement_curr.swap(_displacement_next);
    boundary_condition_second_kind_2d(_right_part, *_mesh, _boundaries_conditions);
    if (_inner_pressure)
        integrate_right_part(_right_part, *_mesh, _inner_pressure);
    _right_part -= _mass.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * _displacement_prev;
    auto tmp = _mass.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * _displacement_curr;
    tmp *= T{2};
    _right_part += tmp;
    _right_part -= _boundary_matrix * calc_first_kind_vector(_mesh->container(), _boundaries_conditions);
    first_kind_fill_2d(_right_part, _mesh->container(), _boundaries_conditions);
    _displacement_next = _slae_solver->solve(_right_part, _displacement_curr);
    _time += time_step();
}

}