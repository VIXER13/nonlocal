#ifndef NONSTATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define NONSTATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "heat_capacity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class nonstationary_heat_equation_solver_1d final {
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>> _slae_solver;
    heat_capacity_matrix_1d<T, I> _capacity;
    thermal_conductivity_matrix_1d<T, I> _conductivity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    const T _time_step = T{1};

public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step);

    const T time_step() const noexcept;
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;

    template<class Init_Dist>
    void compute(const std::vector<equation_parameters<1, T, parameters_1d>>& parameters,
                 const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition,
                 const Init_Dist& init_dist);

    template<class Right_Part>
    void calc_step(const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition,
                   const Right_Part& right_part);
};

template<class T, class I>
nonstationary_heat_equation_solver_1d<T, I>::nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step)
    : _conductivity{mesh}
    , _capacity{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _time_step{time_step} {}

template<class T, class I>
const T nonstationary_heat_equation_solver_1d<T, I>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_heat_equation_solver_1d<T, I>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I>
template<class Init_Dist>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const std::vector<equation_parameters<1, T, parameters_1d>>& parameters,
                                                          const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition,
                                                          const Init_Dist& init_dist) {
    const std::array<bool, 2> is_first_kind = {
        bool(dynamic_cast<stationary_temperature_1d<T>*>(boundary_condition.front().get())),
        bool(dynamic_cast<stationary_temperature_1d<T>*>(boundary_condition.back ().get()))
    };
    _capacity.calc_matrix(parameters, is_first_kind);
    _conductivity.template calc_matrix(parameters, is_first_kind);

    const std::array<size_t, 2> indices = {0, _capacity.mesh().nodes_count() - 1};
    for(const size_t b : std::ranges::iota_view{0u, boundary_condition.size()})
        convection_condition_1d(_conductivity.matrix_inner(), *boundary_condition[b], indices[b]);

    _conductivity.matrix_inner() *= time_step();
    _conductivity.matrix_inner() += _capacity.matrix_inner();
    if (is_first_kind.front())
        _conductivity.matrix_inner().coeffRef(0, 0) = T{1};
    if (is_first_kind.back())
        _conductivity.matrix_inner().coeffRef(_conductivity.matrix_inner().rows() - 1, _conductivity.matrix_inner().cols() - 1) = T{1};
    _slae_solver.compute(_conductivity.matrix_inner());
    
    for(std::unordered_map<size_t, T>& matrix_part : _conductivity.matrix_bound())
        for(auto& [_, val] : matrix_part)
            val *= time_step();

    for(const size_t i : std::ranges::iota_view(0u, _conductivity.mesh().nodes_count()))
        _temperature_prev[i] = init_dist(_conductivity.mesh().node_coord(i));
}

template<class T, class I>
template<class Right_Part>
void nonstationary_heat_equation_solver_1d<T, I>::calc_step(const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition,
                                                            const Right_Part& right_part) {
    _right_part.setZero();
    _temperature_prev.swap(_temperature_curr);
    const std::array<size_t, 2> indices = {0, _capacity.mesh().nodes_count() - 1};
    for(const size_t b : std::ranges::iota_view{0u, boundary_condition.size()})
        boundary_condition_second_kind_1d<T>(_right_part, *boundary_condition[b], indices[b]);
    integrate_right_part(_right_part, _conductivity.mesh(), right_part);
    _right_part *= time_step();
    _right_part += _capacity.matrix_inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    for(const size_t b : std::ranges::iota_view{0u, boundary_condition.size()})
        boundary_condition_first_kind_1d<T>(_right_part, _conductivity.matrix_bound()[b], *boundary_condition[b], indices[b]);
    _temperature_curr = _slae_solver.solve(_right_part);
}

}

#endif