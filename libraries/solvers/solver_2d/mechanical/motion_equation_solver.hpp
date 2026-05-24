#pragma once

#include "stiffness_matrix_2d.hpp"
#include "mass_matrix.hpp"
#include "mechanical_parameters_2d.hpp"

#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<class T, class I, class Matrix_Index>
class motion_equation_solver final {
    static constexpr size_t DoF = 2;

    std::unique_ptr<slae::conjugate_gradient<T, Matrix_Index>> slae_solver;
    mass_matrix<T, I, Matrix_Index> _mass;
    stiffness_matrix<T, I, Matrix_Index> _stiffness;
    mechanical_boundaries_conditions_2d<T> _boundaries_conditions
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _displacement_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _displacement_curr;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _displacement_next;
    const T _time_step = 1;

public:
    explicit motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& displacement() const noexcept;
    constexpr T time_step() const noexcept;

    void compute(const raw_mechanical_parameters<T>& parameters,
                 mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                 const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>>& init_dist);

    void calc_step(const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>>& right_part);
};

template<class T, class I, class Matrix_Index>
motion_equation_solver<T, I, Matrix_Index>::motion_equation_solver(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step)
    : _conductivity{mesh}
    , _mass{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _displacement_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _displacement_next{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DoF * mesh->container().nodes_count())}
    , _time_step{time_step} {}

template<class T, class I, class Matrix_Index>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& motion_equation_solver<T, I, Matrix_Index>::temperature() const noexcept {
    return _displacement_next;
}

template<class T, class I, class Matrix_Index>
constexpr T motion_equation_solver<T, I, Matrix_Index>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I, class Matrix_Index>
template<class Init_Dist>
void motion_equation_solver<T, I, Matrix_Index>::compute(const raw_mechanical_parameters<T>& parameters,
                                                         mechanical_boundaries_conditions_2d<T>&& boundaries_conditions,
                                                         const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>>& init_dist) {
    _boundaries_conditions = std::move(boundaries_conditions);
    const auto& mesh = _mass.mesh();
    const auto settings = init_problem_settings(mesh.container(), parameters, _boundaries_conditions);
    log_problem_settings(settings);
    const auto evaluated_parameters = evaluate_mechanical_parameters(mesh, parameters);
    _mass.compute(settings.is_inner_nodes);
    _stiffness.compute(evaluated_parameters, settings);

    _stiffness.matrix().inner() *= time_step() * time_step();
    _stiffness.matrix().bound() *= time_step() * time_step();
    _stiffness.matrix().inner() -= _mass.matrix().inner();

    if (init_dist) {
        for(const size_t node : mesh.container().nodes())
            _displacement_next[node] = (*init_dist)(mesh.container().node_coord(node));
        _displacement_curr = _displacement_next;
    }

    slae_solver = std::make_unique<slae::conjugate_gradient<T, Matrix_Index>>(_stiffness.matrix().inner());
}

template<class T, class I, class Matrix_Index>
void motion_equation_solver<T, I, Matrix_Index>::calc_step(const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>>& right_part) {
    _right_part.setZero();
    _displacement_prev.swap(_displacement_curr);
    _displacement_curr.swap(_displacement_next)
    const auto& mesh = _mass.mesh();
    boundary_condition_second_kind_2d(_right_part, mesh, _boundaries_conditions);
    if (right_part)
        integrate_right_part<DoF>(_right_part, mesh, *right_part);
    _right_part *= -time_step() * time_step();
    _right_part += _mass.matrix().inner().template selfadjointView<Eigen::Upper>() * _displacement_prev;
    _right_part -= 2 * _mass.matrix().inner().template selfadjointView<Eigen::Upper>() * _displacement_curr;
    boundary_condition_first_kind_2d(_right_part, mesh, _boundaries_conditions, _stiffness.matrix().bound());
    _displacement_next = slae_solver->solve(_right_part, _displacement_curr);
}

}