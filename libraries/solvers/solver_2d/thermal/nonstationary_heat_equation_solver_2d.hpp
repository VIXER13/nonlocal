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

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class nonstationary_heat_equation_solver_2d final {
    static constexpr size_t DoF = 1;

    std::unique_ptr<slae::conjugate_gradient<T, Matrix_Index>> slae_solver;
    heat_capacity_matrix_2d<T, I, Matrix_Index> _capacity;
    conductivity_matrix_2d<T, I, Matrix_Index> _conductivity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    const T _time_step = 1;

public:
    explicit nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    constexpr T time_step() const noexcept;

    template<class Init_Dist>
    void compute(const parameters_2d<T>& parameters,
                 const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                 const Init_Dist& init_dist);

    template<class Right_Part>
    void calc_step(const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                   const Right_Part& right_part);
};

template<class T, class I, class Matrix_Index>
nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step)
    : _conductivity{mesh}
    , _capacity{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count())}
    , _time_step{time_step} {}

template<class T, class I, class Matrix_Index>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I, class Matrix_Index>
constexpr T nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I, class Matrix_Index>
template<class Init_Dist>
void nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::compute(const parameters_2d<T>& parameters,
                                                                        const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                        const Init_Dist& init_dist) {
    const std::vector<bool> is_inner = utils::inner_nodes(_conductivity.mesh().container(), boundaries_conditions);
    _conductivity.compute(parameters, is_inner);
    convection_condition_2d(_conductivity.matrix().inner(), _conductivity.mesh(), boundaries_conditions);
    _capacity.calc_matrix(parameters, is_inner);

    _conductivity.matrix().inner() *= time_step();
    _conductivity.matrix().bound() *= time_step();
    _conductivity.matrix().inner() += _capacity.matrix().inner();
    first_kind_filler(_conductivity.mesh().process_nodes(), is_inner, [&matrix = _conductivity.matrix().inner()](const size_t row) {
        matrix.valuePtr()[matrix.outerIndexPtr()[row]] = T{1};
    });

    for(const size_t node : _conductivity.mesh().container().nodes())
        _temperature_curr[node] = init_dist(_conductivity.mesh().container().node_coord(node));

    slae_solver = std::make_unique<slae::conjugate_gradient<T, Matrix_Index>>(_conductivity.matrix().inner());
}

template<class T, class I, class Matrix_Index>
template<class Right_Part>
void nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::calc_step(const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                          const Right_Part& right_part) {
    _right_part.setZero();
    _temperature_prev.swap(_temperature_curr);
    radiation_condition_2d(_conductivity.matrix().inner(), _right_part, _conductivity.mesh(), boundaries_conditions, 
                           _temperature_prev, time_step());

    boundary_condition_second_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions);
    integrate_right_part<DoF>(_right_part, _conductivity.mesh(), right_part);

    _right_part *= time_step();
    _right_part += _capacity.matrix().inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    boundary_condition_first_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions, _conductivity.matrix().bound());
    _temperature_curr = slae_solver->solve(_right_part, _temperature_prev);
}

}