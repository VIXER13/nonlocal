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

template<class T>
class motion_equation_solver final {
    static constexpr size_t DoF = 2;

    std::unique_ptr<slae::iterative_solver_base<T>> slae_solver;
    mass_matrix<T> _mass;
    stiffness_matrix<T> _stiffness;
    mechanical_boundaries_conditions_2d<T> _boundaries_conditions;
    evaluated_mechanical_parameters<T> _parameters;
    std::vector<T> _right_part;
    std::vector<T> _displacement_prev;
    std::vector<T> _displacement_curr;
    std::vector<T> _displacement_next;
    T _time_step = T{1};
    T _time = T{0};

public:
    explicit motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);

    const std::vector<T>& displacement() const noexcept;
    mechanical_solution_2d<T> solution(const bool strain_and_stress = true) const;
    T time_step() const noexcept;
    T time() const noexcept;

    void compute(const raw_mechanical_parameters<T>& parameters,
                 mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                 const T time_step, const T time = T{0},
                 const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& init_dist = std::nullopt);

    void calc_step(const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& right_part = std::nullopt);
};

template<class T>
motion_equation_solver<T>::motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _mass{mesh}
    , _stiffness{mesh} 
    , _right_part(DoF * mesh->container().nodes_count(), T{0})
    , _displacement_prev(DoF * mesh->container().nodes_count(), T{0})
    , _displacement_curr(DoF * mesh->container().nodes_count(), T{0})
    , _displacement_next(DoF * mesh->container().nodes_count(), T{0}) {}

template<class T>
const std::vector<T>& motion_equation_solver<T>::displacement() const noexcept {
    return _displacement_next;
}

template<class T>
mechanical_solution_2d<T> motion_equation_solver<T>::solution(const bool strain_and_stress) const {
    mechanical_solution_2d<T> sol{_mass.mesh_ptr(), _parameters, displacement()};
    if (strain_and_stress)
        sol.calc_strain_and_stress(_parameters);
    return sol;
}

template<class T>
T motion_equation_solver<T>::time_step() const noexcept {
    return _time_step;
}

template<class T>
T motion_equation_solver<T>::time() const noexcept {
    return _time;
}

template<class T>
void motion_equation_solver<T>::compute(const raw_mechanical_parameters<T>& parameters,
                                        mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                                        const T time_step, const T time,
                                        const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& init_dist) {
    _time_step = time_step;
    _time = time;
    _boundaries_conditions = std::move(boundaries_conditions);
    const auto& mesh = _mass.mesh();
    const auto settings = init_problem_settings(mesh.container(), parameters, _boundaries_conditions);
    log_problem_settings(settings);
    _parameters = evaluate_mechanical_parameters(mesh, parameters);
    _mass.compute(_parameters, settings.is_inner_nodes);
    _stiffness.compute(_parameters, settings);

    _mass.matrix().inner() /= time_step * time_step;
    if (settings.is_symmetric())
        _stiffness.matrix().inner() += _mass.matrix().inner();
    else
        _stiffness.matrix().inner() += _mass.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>();
    first_kind_filler(_mass.mesh().process_nodes(), settings.is_inner_nodes, [&matrix = _stiffness.matrix().inner()](const size_t row) {
        matrix.values[matrix.portrait.shifts[row]] = T{1};
    });

    if (init_dist) {
        for(const size_t node : mesh.container().nodes()) {
            const auto displacement = (*init_dist)(mesh.container().node_coord(node));
            _displacement_next[node + X] = displacement[X];
            _displacement_next[node + Y] = displacement[Y];
        }
        _displacement_curr = _displacement_next;
    }

    slae_solver = slae::init_iterative_solver(_stiffness.matrix().inner(), settings.is_symmetric());
}

template<class T>
void motion_equation_solver<T>::calc_step(const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& right_part) {
    std::fill(_right_part.begin(), _right_part.end(), T{0});
    _displacement_prev.swap(_displacement_curr);
    _displacement_curr.swap(_displacement_next);
    const auto& mesh = _mass.mesh();
    boundary_condition_second_kind_2d(_right_part, mesh, _boundaries_conditions);
    if (right_part)
        integrate_right_part<DoF>(_right_part, mesh, *right_part);
    _right_part -= _mass.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * _displacement_prev;
    const std::vector<T> tmp = _mass.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * _displacement_curr;
    _right_part += T{2} * tmp;
    boundary_condition_first_kind_2d(_right_part, mesh, _boundaries_conditions, _stiffness.matrix().bound());
    _displacement_next = slae_solver->solve(_right_part, _displacement_curr);
    _time += time_step();
}

}