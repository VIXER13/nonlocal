#ifndef NONSTATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define NONSTATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "heat_capacity_matrix_1d.hpp"
#include "convection_condition_1d.hpp"
#include "right_part_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class nonstationary_heat_equation_solver_1d final {
    thermal_conductivity_matrix_1d<T, I> _conductivity;
    heat_capacity_matrix_1d<T, I> _capacity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    T _time = 0;
    const T _tau = 1;

public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T tau, const T time = 0);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;

    template<class Init_Dist, class Influence_Function>
    void compute(const heat_equation_parameters_1d<T>& equation_param,
                 const std::array<boundary_condition_t, 2> boundary_condition,
                 const Init_Dist& init_dist,
                 const T p1,
                 const Influence_Function& influence_function);

    template<class Right_Part>
    void calc_step(const std::array<nonstatinary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                   const Right_Part& right_part);
};

template<class T, class I>
nonstationary_heat_equation_solver_1d<T, I>::nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T tau, const T time)
    : _conductivity{mesh}
    , _capacity{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _time{time}
    , _tau{tau} {}

template<class T, class I>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_heat_equation_solver_1d<T, I>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I>
template<class Init_Dist, class Influence_Function>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const heat_equation_parameters_1d<T>& equation_param,
                                                          const std::array<boundary_condition_t, 2> boundary_condition,
                                                          const Init_Dist& init_dist,
                                                          const T p1,
                                                          const Influence_Function& influence_function) {
    _conductivity.template calc_matrix(equation_param.lambda, p1, influence_function, boundary_condition);
    convection_condition_1d(_conductivity.matrix_inner(), boundary_condition, equation_param.alpha);
    _capacity.calc_matrix(equation_param.c, equation_param.rho, boundary_condition);

    _conductivity.matrix_inner() *= _tau;
    _conductivity.matrix_inner() += _capacity.matrix_inner();
    if (boundary_condition.front() == boundary_condition_t::TEMPERATURE)
        _conductivity.matrix_inner().coeffRef(0, 0) = T{1};
    if (boundary_condition.back() == boundary_condition_t::TEMPERATURE)
        _conductivity.matrix_inner().coeffRef(_conductivity.matrix_inner().rows()-1, _conductivity.matrix_inner().cols()-1) = T{1};
    for(std::unordered_map<size_t, T>& matrix_part : _conductivity.matrix_bound())
        for(auto& [_, val] : matrix_part)
            val *= _tau;

    for(const size_t i : std::views::iota(size_t{0}, _conductivity.mesh()->nodes_count()))
        _temperature_prev[i] = init_dist(_conductivity.mesh()->node_coord(i));
}

template<class T, class I>
template<class Right_Part>
void nonstationary_heat_equation_solver_1d<T, I>::calc_step(const std::array<nonstatinary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                                            const Right_Part& right_part) {
    _time += _tau;
    _right_part.setZero();
    const auto stationary_bound = to_stationary(boundary_condition, _time);
    boundary_condition_second_kind_1d(_right_part, stationary_bound, std::array{size_t{0}, size_t(_right_part.size() - 1)});
    integrate_right_part(_right_part, *_conductivity.mesh(), [&right_part, t = _time](const T x) { return right_part(t, x); });
    _right_part *= _tau;
    _right_part += _capacity.matrix_inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    boundary_condition_first_kind_1d(_right_part, stationary_bound, _conductivity.matrix_bound());
    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{_conductivity.matrix_inner()};
    _temperature_curr = solver.template solveWithGuess(_right_part, _temperature_prev);
    _temperature_prev.swap(_temperature_curr);
}

}

#endif