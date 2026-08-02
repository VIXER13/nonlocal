#pragma once

#include "conductivity_matrix_2d.hpp"
#include "heat_capacity_matrix_2d.hpp"
#include "convection_condition_2d.hpp"
#include "radiation_condition_2d.hpp"
#include "thermal_parameters_2d.hpp"

#include <solvers/slae/conjugate_gradient.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
class nonstationary_heat_equation_solver_2d final {
    static constexpr size_t DoF = 1;

    std::shared_ptr<mesh::mesh_2d<T>> _mesh;
    std::unique_ptr<slae::conjugate_gradient<T>> slae_solver;
    heat_capacity_matrix_2d<T> _capacity;
    conductivity_matrix_2d<T> _conductivity;
    metamath::linear::sparse_matrix<T> _boundary_matrix;
    std::vector<T> _right_part;
    std::vector<T> _temperature_curr;
    std::vector<T> _temperature_next;

    std::function<T(const std::array<T, 2>&)> _inner_flux;
    thermal_boundaries_conditions_2d<T> _boundaries_conditions;
    evaluated_conductivity_2d<T> _parameters;
    T _time_step = T{1};
    T _time = T{0};

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
                 const T time = T{0});

    void calc_step();
};

template<std::floating_point T>
nonstationary_heat_equation_solver_2d<T>::nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _mesh{mesh}
    , _conductivity{*mesh}
    , _capacity{*mesh}
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
                                                       const T time) {
    // TODO: Rework nonstationary solver

    // std::vector<T> solution(_conductivity.mesh().quad_shift(_conductivity.mesh().container().elements_2d_count()), T{0});
    // const auto conductivity_parameters = evaluate_conductivity(_conductivity.mesh(), parameters, solution);
    // solution = {};

    // const std::vector<bool> is_inner = utils::inner_nodes(_conductivity.mesh().container(), boundaries_conditions);
    // _conductivity.compute(conductivity_parameters, is_inner);
    // convection_condition_2d(_conductivity.matrix().inner(), _conductivity.mesh(), boundaries_conditions, is_inner);
    // _capacity.calc_matrix(parameters, is_inner);

    // _conductivity.matrix().inner() *= time_step();
    // _conductivity.matrix().bound() *= time_step();
    // _conductivity.matrix().inner() += _capacity.matrix().inner();
    // first_kind_filler(_conductivity.mesh().process_nodes(), is_inner, [&matrix = _conductivity.matrix().inner()](const size_t row) {
    //     matrix.values[matrix.portrait.shifts[row]] = T{1};
    // });

    // for(const size_t node : _conductivity.mesh().container().nodes())
    //     _temperature_next[node] = init_dist(_conductivity.mesh().container().node_coord(node));

    // slae_solver = std::make_unique<slae::conjugate_gradient<T>>(_conductivity.matrix().inner());
}

template<std::floating_point T>
void nonstationary_heat_equation_solver_2d<T>::calc_step() {
    // std::fill(_right_part.begin(), _right_part.end(), T{0});
    // _temperature_curr.swap(_temperature_next);
    // radiation_condition_2d(_conductivity.matrix().inner(), _right_part, _conductivity.mesh(), boundaries_conditions, 
    //                        _temperature_curr, time_step());

    // boundary_condition_second_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions);
    // integrate_right_part<DoF>(_right_part, _conductivity.mesh(), right_part);

    // using namespace metamath::operators;
    // _right_part *= time_step();
    // _right_part += _capacity.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * _temperature_curr;
    // boundary_condition_first_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions, _conductivity.matrix().bound());
    // _temperature_next = slae_solver->solve(_right_part, _temperature_curr);
}

}