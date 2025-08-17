#pragma once

#include "thermal_conductivity_assembler.hpp"
#include "heat_capacity_assembler.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"
#include "init_problem_settings.hpp"

#include <solvers/solver_1d/base/assemble_matrix_portrait.hpp>
#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_second_kind_1d.hpp>
#include <solvers/solver_1d/base/right_part_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T, class I>
class nonstationary_heat_equation_solver_1d final {
    static constexpr bool Is_Stationary = false;

    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    finite_element_matrix_1d<T, I> _capacity;
    finite_element_matrix_1d<T, I> _conductivity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    std::array<T, 2> _capacity_initial_values = {T(0), T(0)};
    std::array<T, 2> _conductivity_initial_values = {T(0), T(0)};
    const T _time_step = T{1};

    static std::array<T, 2> get_initial_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);
    static void reset_initial_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                     const std::array<T, 2>& values);
 
public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step);

    const T time_step() const noexcept;
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;

    void compute(const parameters_1d<T>& parameters,
                 const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                 const std::optional<std::function<T(const T)>>& initial_temperature = std::nullopt);

    void calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                   const std::optional<std::function<T(const T)>>& right_part = std::nullopt);
};

template<class T, class I>
nonstationary_heat_equation_solver_1d<T, I>::nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step)
    : _mesh{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _time_step{time_step} {}

template<class T, class I>
std::array<T, 2> nonstationary_heat_equation_solver_1d<T, I>::get_initial_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
    return {matrix.coeff(0, 0), matrix.coeff(matrix.rows() - 1, matrix.cols() - 1)};
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::reset_initial_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                                                       const std::array<T, 2>& values) {
    matrix.coeffRef(0, 0) = values.front();
    matrix.coeffRef(matrix.rows() - 1, matrix.cols() - 1) = values.back();
}

template<class T, class I>
const T nonstationary_heat_equation_solver_1d<T, I>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_heat_equation_solver_1d<T, I>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I>
const mesh::mesh_1d<T>& nonstationary_heat_equation_solver_1d<T, I>::mesh() const {
    return *_mesh;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& nonstationary_heat_equation_solver_1d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const parameters_1d<T>& parameters,
                                                          const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                          const std::optional<std::function<T(const T)>>& initial_temperature) {
    problem_settings settings = init_problem_settings(parameters, boundaries_conditions, Is_Stationary);

    init_matrix_portrait(_conductivity.inner, mesh(), settings);
    thermal_conductivity_assembler_1d<T, I> conductivity_assembler{_conductivity, _mesh};
    _conductivity.set_zero();
    conductivity_assembler.calc_matrix(parameters, settings);
    convection_condition_1d(_conductivity.inner, boundaries_conditions);

    settings.theories = std::vector<theory_t>(parameters.size(), theory_t::LOCAL);
    init_matrix_portrait(_capacity.inner, mesh(), settings);
    _capacity.set_zero();
    heat_capacity_assembler_1d<T, I> capacity_assembler{_capacity, _mesh};
    capacity_assembler.calc_matrix(parameters, settings.is_first_kind);
    _capacity_initial_values = get_initial_values(_capacity.inner);

    _conductivity.inner *= time_step();
    _conductivity.inner += _capacity.inner;
    _conductivity_initial_values = get_initial_values(_conductivity.inner);
    reset_initial_values(_conductivity.inner, {
        settings.is_first_kind.front() ? T{1} : _conductivity_initial_values.front(),
        settings.is_first_kind.back()  ? T{1} : _conductivity_initial_values.back()
    });
    _conductivity_initial_values = get_initial_values(_conductivity.inner);
    
    for(std::unordered_map<size_t, T>& matrix_part : _conductivity.bound)
        for(auto& val : matrix_part | std::views::values)
            val *= time_step();

    if (initial_temperature)
        for(const size_t i : std::ranges::iota_view{0u, mesh().nodes_count()})
            _temperature_curr[i] = (*initial_temperature)(mesh().node_coord(i));
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                            const std::optional<std::function<T(const T)>>& right_part) {
    _right_part.setZero();
    _temperature_prev.swap(_temperature_curr);
    reset_initial_values(_conductivity.inner, _conductivity_initial_values);
    reset_initial_values(_capacity.inner, _capacity_initial_values);

    boundary_condition_second_kind_1d<T>(_right_part, boundaries_conditions);
    if (right_part)
        integrate_right_part(_right_part, mesh(), *right_part);
    _right_part *= time_step();
    _right_part += _capacity.inner.template selfadjointView<Eigen::Upper>() * _temperature_prev;
    
    radiation_condition_1d<Is_Stationary>(_conductivity.inner, _right_part, boundaries_conditions, _temperature_prev, time_step());
    boundary_condition_first_kind_1d(_right_part, _conductivity.bound, boundaries_conditions);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{_conductivity.inner};
    _temperature_curr = solver.solveWithGuess(_right_part, _temperature_prev);
}

}