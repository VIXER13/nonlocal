#pragma once

#include "stiffness_matrix_2d.hpp"
#include "mass_matrix.hpp"
#include "mechanical_parameters_2d.hpp"
#include "mechanical_solution_2d.hpp"

#include <solvers/slae/init_solver_method.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<class T, class I, class Matrix_Index>
class motion_equation_solver final {
    static constexpr size_t DoF = 2;

    std::unique_ptr<slae::iterative_solver_base<T, Matrix_Index>> slae_solver;
    mass_matrix<T, I, Matrix_Index> _mass;
    stiffness_matrix<T, I, Matrix_Index> _stiffness;
    mechanical_boundaries_conditions_2d<T> _boundaries_conditions;
    evaluated_mechanical_parameters<T> _parameters;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _displacement_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _displacement_curr;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _displacement_next;
    T _time_step = T{1};
    T _time = T{0};

public:
    explicit motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& displacement() const noexcept;
    mechanical_solution_2d<T> solution(const bool strain_and_stress = true) const;
    T time_step() const noexcept;
    T time() const noexcept;

    void compute(const raw_mechanical_parameters<T>& parameters,
                 mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                 const T time_step, const T time = T{0},
                 const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& init_dist = std::nullopt);

    void calc_step(const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& right_part = std::nullopt);
};

template<class T, class I, class Matrix_Index>
motion_equation_solver<T, I, Matrix_Index>::motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _mass{mesh}
    , _stiffness{mesh} 
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _displacement_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _displacement_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _displacement_next{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())} {}

template<class T, class I, class Matrix_Index>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& motion_equation_solver<T, I, Matrix_Index>::displacement() const noexcept {
    return _displacement_next;
}

template<class T, class I, class Matrix_Index>
mechanical_solution_2d<T> motion_equation_solver<T, I, Matrix_Index>::solution(const bool strain_and_stress) const {
    mechanical_solution_2d<T> sol{_mass.mesh_ptr(), _parameters, displacement()};
    if (strain_and_stress)
        sol.calc_strain_and_stress(_parameters);
    return sol;
}

template<class T, class I, class Matrix_Index>
T motion_equation_solver<T, I, Matrix_Index>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I, class Matrix_Index>
T motion_equation_solver<T, I, Matrix_Index>::time() const noexcept {
    return _time;
}

template<class T, class I, class Matrix_Index>
void motion_equation_solver<T, I, Matrix_Index>::compute(const raw_mechanical_parameters<T>& parameters,
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
        _stiffness.matrix().inner() += _mass.matrix().inner().template selfadjointView<Eigen::Upper>();
    first_kind_filler(_mass.mesh().process_nodes(), settings.is_inner_nodes, [&matrix = _stiffness.matrix().inner()](const size_t row) {
        matrix.valuePtr()[matrix.outerIndexPtr()[row]] = T{1};
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

template<class T, class I, class Matrix_Index>
void motion_equation_solver<T, I, Matrix_Index>::calc_step(const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& right_part) {
    _right_part.setZero();
    _displacement_prev.swap(_displacement_curr);
    _displacement_curr.swap(_displacement_next);
    const auto& mesh = _mass.mesh();
    boundary_condition_second_kind_2d(_right_part, mesh, _boundaries_conditions);
    if (right_part)
        integrate_right_part<DoF>(_right_part, mesh, *right_part);
    _right_part -= _mass.matrix().inner().template selfadjointView<Eigen::Upper>() * _displacement_prev;
    const Eigen::Matrix<T, Eigen::Dynamic, 1> tmp = _mass.matrix().inner().template selfadjointView<Eigen::Upper>() * _displacement_curr;
    _right_part += T{2} * tmp;
    boundary_condition_first_kind_2d(_right_part, mesh, _boundaries_conditions, _stiffness.matrix().bound());
    _displacement_next = slae_solver->solve(_right_part, _displacement_curr);
    _time += time_step();
}

}