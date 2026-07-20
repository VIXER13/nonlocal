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

template<class T, std::integral I>
class nonstationary_heat_equation_solver_2d final {
    static constexpr size_t DoF = 1;

    std::unique_ptr<slae::conjugate_gradient<T>> slae_solver;
    heat_capacity_matrix_2d<T, I> _capacity;
    conductivity_matrix_2d<T, I> _conductivity;
    std::vector<T> _right_part;
    std::vector<T> _temperature_prev;
    std::vector<T> _temperature_curr;
    const T _time_step = 1;

public:
    explicit nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step);

    const std::vector<T>& temperature() const noexcept;
    constexpr T time_step() const noexcept;

    template<class Init_Dist>
    void compute(const parameters_2d<T>& parameters,
                 const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                 const Init_Dist& init_dist);

    template<class Right_Part>
    void calc_step(const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                   const Right_Part& right_part);
};

template<class T, std::integral I>
nonstationary_heat_equation_solver_2d<T, I>::nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step)
    : _conductivity{mesh}
    , _capacity{mesh}
    , _right_part(mesh->container().nodes_count(), T{0})
    , _temperature_prev(mesh->container().nodes_count(), T{0})
    , _temperature_curr(mesh->container().nodes_count(), T{0})
    , _time_step{time_step} {}

template<class T, std::integral I>
const std::vector<T>& nonstationary_heat_equation_solver_2d<T, I>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, std::integral I>
constexpr T nonstationary_heat_equation_solver_2d<T, I>::time_step() const noexcept {
    return _time_step;
}

template<class T, std::integral I>
template<class Init_Dist>
void nonstationary_heat_equation_solver_2d<T, I>::compute(const parameters_2d<T>& parameters,
                                                          const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                          const Init_Dist& init_dist) {
    std::vector<T> solution(_conductivity.mesh().quad_shift(_conductivity.mesh().container().elements_2d_count()), T{0});
    const auto conductivity_parameters = evaluate_conductivity(_conductivity.mesh(), parameters, solution);
    solution = {};

    const std::vector<bool> is_inner = utils::inner_nodes(_conductivity.mesh().container(), boundaries_conditions);
    _conductivity.compute(conductivity_parameters, is_inner);
    convection_condition_2d(_conductivity.matrix().inner(), _conductivity.mesh(), boundaries_conditions, is_inner);
    _capacity.calc_matrix(parameters, is_inner);

    _conductivity.matrix().inner() *= time_step();
    _conductivity.matrix().bound() *= time_step();
    _conductivity.matrix().inner() += _capacity.matrix().inner();
    first_kind_filler(_conductivity.mesh().process_nodes(), is_inner, [&matrix = _conductivity.matrix().inner()](const size_t row) {
        matrix.values[matrix.portrait.shifts[row]] = T{1};
    });

    for(const size_t node : _conductivity.mesh().container().nodes())
        _temperature_curr[node] = init_dist(_conductivity.mesh().container().node_coord(node));

    slae_solver = std::make_unique<slae::conjugate_gradient<T>>(_conductivity.matrix().inner());
}

template<class T, std::integral I>
template<class Right_Part>
void nonstationary_heat_equation_solver_2d<T, I>::calc_step(const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                            const Right_Part& right_part) {
    std::fill(_right_part.begin(), _right_part.end(), T{0});
    _temperature_prev.swap(_temperature_curr);
    radiation_condition_2d(_conductivity.matrix().inner(), _right_part, _conductivity.mesh(), boundaries_conditions, 
                           _temperature_prev, time_step());

    boundary_condition_second_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions);
    integrate_right_part<DoF>(_right_part, _conductivity.mesh(), right_part);

    using namespace metamath::operators;
    _right_part *= time_step();
    _right_part += _capacity.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * _temperature_prev;
    boundary_condition_first_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions, _conductivity.matrix().bound());
    _temperature_curr = slae_solver->solve(_right_part, _temperature_prev);
}

}