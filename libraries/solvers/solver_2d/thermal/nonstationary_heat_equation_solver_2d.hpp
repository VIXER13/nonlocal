#ifndef NONSTATIONARY_HEAT_EQUATION_SOLVER_2D_HPP
#define NONSTATIONARY_HEAT_EQUATION_SOLVER_2D_HPP

#include "thermal_conductivity_matrix_2d.hpp"
#include "heat_capacity_matrix_2d.hpp"
#include "right_part_2d.hpp"
#include "boundary_condition_first_kind_2d.hpp"
#include "boundary_condition_second_kind_2d.hpp"
#include "convection_condition_2d.hpp"
#include "thermal_parameters_2d.hpp"

#include "conjugate_gradient.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
class nonstationary_heat_equation_solver_2d final {
    static constexpr size_t DoF = 1;

    std::unique_ptr<slae::conjugate_gradient<T, Matrix_Index>> slae_solver;
    heat_capacity_matrix_2d<T, I, Matrix_Index> _capacity;
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> _conductivity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    const T _time_step = 1;

public:
    explicit nonstationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const T time_step);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    constexpr T time_step() const noexcept;

    template<class Init_Dist, class Influence_Function>
    void compute(const parameter_2d<T>& parameters,
                 const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                 const Init_Dist& init_dist,
                 const T p1, const Influence_Function& influence_function);

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
template<class Init_Dist, class Influence_Function>
void nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::compute(const parameter_2d<T>& parameters,
                                                                        const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                        const Init_Dist& init_dist,
                                                                        const T p1, const Influence_Function& influence_function) {
    const std::vector<bool> is_inner = utils::inner_nodes(_conductivity.mesh().container(), boundaries_conditions);
    _conductivity.template compute(parameters.conductivity, parameters.material, is_inner, p1, influence_function);
    convection_condition_2d(_conductivity.matrix_inner(), _conductivity.mesh(), boundaries_conditions);
    _capacity.calc_matrix(parameters.capacity, parameters.density, is_inner);

    _conductivity.matrix_inner() *= time_step();
    _conductivity.matrix_bound() *= time_step();
    _conductivity.matrix_inner() += _capacity.matrix_inner();
    for(const size_t node : std::ranges::iota_view{0u, is_inner.size()})
        if(!is_inner[node])
            _conductivity.matrix_inner().coeffRef(node, node) = T{1};

    for(const size_t node : _conductivity.mesh().container().nodes())
        _temperature_prev[node] = init_dist(_conductivity.mesh().container().node_coord(node));

    slae_solver = std::make_unique<slae::conjugate_gradient<T, Matrix_Index>>(_conductivity.matrix_inner());
}

template<class T, class I, class Matrix_Index>
template<class Right_Part>
void nonstationary_heat_equation_solver_2d<T, I, Matrix_Index>::calc_step(const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                          const Right_Part& right_part) {
    _right_part.setZero();
    _temperature_prev.swap(_temperature_curr);
    boundary_condition_second_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions);
    integrate_right_part<DoF>(_right_part, _conductivity.mesh(), right_part);
    _right_part *= time_step();
    _right_part += _capacity.matrix_inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    boundary_condition_first_kind_2d(_right_part, _conductivity.mesh(), boundaries_conditions, _conductivity.matrix_bound());
    _temperature_curr = slae_solver->solve(_right_part, _temperature_prev);
}

}

#endif