#pragma once

#include "thermal_conductivity_assembler.hpp"
#include "heat_capacity_assembler.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_second_kind_1d.hpp>
#include <solvers/solver_1d/base/right_part_1d.hpp>

namespace nonlocal::thermal {

template<class T, class I>
class nonstationary_heat_equation_solver_1d final {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    finite_element_matrix_1d<T, I> _capacity;
    finite_element_matrix_1d<T, I> _conductivity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    const T _time_step = T{1};
    std::array<T, 2> _capacity_initial_values = {T(0), T(0)};
    std::array<T, 2> _conductivity_initial_values = {T(0), T(0)};

    static std::array<T, 2> get_init_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);
    static void reset_to_init_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                     const std::array<T, 2>& values);
 
public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step);

    const T time_step() const noexcept;
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;

    template<class Init_Dist>
    void compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                 const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                 const Init_Dist& init_dist);

    template<class Right_Part>
    void calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                   const Right_Part& right_part);
};

template<class T, class I>
nonstationary_heat_equation_solver_1d<T, I>::nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step)
    : _mesh{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _time_step{time_step} {}

template<class T, class I>
std::array<T, 2> nonstationary_heat_equation_solver_1d<T, I>::get_init_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
    return {matrix.coeff(0, 0), matrix.coeff(matrix.rows() - 1, matrix.cols() - 1)};
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::reset_to_init_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
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
template<class Init_Dist>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                                                          const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                          const Init_Dist& init_dist) {
    const std::array<bool, 2> is_first_kind = {
        bool(dynamic_cast<const temperature_1d<T>*>(boundaries_conditions.front().get())),
        bool(dynamic_cast<const temperature_1d<T>*>(boundaries_conditions.back ().get()))
    };
    heat_capacity_assembler_1d<T, I> capacity_assembler{_capacity, _mesh};
    thermal_conductivity_assembler_1d<T, I> conductivity_assembler{_conductivity, _mesh};
    capacity_assembler.calc_matrix(parameters, is_first_kind);
    conductivity_assembler.template calc_matrix(parameters, is_first_kind);
    convection_condition_1d(_conductivity.inner, boundaries_conditions);

    auto& conductivity = _conductivity.inner;
    conductivity *= time_step();
    conductivity += _capacity.inner;
    reset_to_init_values(conductivity, {
        is_first_kind.front() ? T{1} : conductivity.coeffRef(0, 0),
        is_first_kind.back()  ? T{1} : conductivity.coeffRef(conductivity.rows() - 1, conductivity.cols() - 1)
    });
    
    for(std::unordered_map<size_t, T>& matrix_part : _conductivity.bound)
        for(auto& val : matrix_part | std::views::values)
            val *= time_step();

    _conductivity_initial_values = get_init_values(_conductivity.inner);
    _capacity_initial_values = get_init_values(_capacity.inner);

    if constexpr (!std::is_same_v<Init_Dist, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        for(const size_t i : std::ranges::iota_view{0u, mesh().nodes_count()})
            _temperature_curr[i] = init_dist(mesh().node_coord(i));
}

template<class T, class I>
template<class Right_Part>
void nonstationary_heat_equation_solver_1d<T, I>::calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                            const Right_Part& right_part) {
    _right_part.setZero();
    _temperature_prev.swap(_temperature_curr);
    reset_to_init_values(_conductivity.inner, _conductivity_initial_values);
    reset_to_init_values(_capacity.inner, _capacity_initial_values);

    boundary_condition_second_kind_1d<T>(_right_part, boundaries_conditions);
    if constexpr (!std::is_same_v<Right_Part, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        integrate_right_part(_right_part, mesh(), right_part);
    _right_part *= time_step();
    _right_part += _capacity.inner.template selfadjointView<Eigen::Upper>() * _temperature_prev;
    
    radiation_condition_1d(_conductivity.inner, _right_part, boundaries_conditions, _temperature_prev, time_step());
    boundary_condition_first_kind_1d(_right_part, _conductivity.bound, boundaries_conditions);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{_conductivity.inner};
    _temperature_curr = solver.solveWithGuess(_right_part, _temperature_prev);
}

}