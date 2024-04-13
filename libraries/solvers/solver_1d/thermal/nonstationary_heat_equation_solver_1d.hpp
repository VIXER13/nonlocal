#pragma once

#include "matrix_portrait_assembler_1d.hpp"
#include "thermal_conductivity_matrix_assembler_1d.hpp"
#include "heat_capacity_matrix_assembler_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

#include "nonlocal_constants.hpp"

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

    static std::array<T, 2> get_initial_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix) {
        return {matrix.coeff(0, 0), matrix.coeff(matrix.rows() - 1, matrix.cols() - 1)};
    }

    static void reset_initial_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                     const std::array<T, 2>& values) {
        matrix.coeffRef(0, 0) = values.front();
        matrix.coeffRef(matrix.rows() - 1, matrix.cols() - 1) = values.back();
    }
 
public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const T time_step);

    finite_element_matrix_1d<T, I>& capacity() noexcept;
    finite_element_matrix_1d<T, I>& conductivity() noexcept;
    const finite_element_matrix_1d<T, I>& capacity() const noexcept;
    const finite_element_matrix_1d<T, I>& conductivity() const noexcept;
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    const T time_step() const noexcept;

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
finite_element_matrix_1d<T, I>& nonstationary_heat_equation_solver_1d<T, I>::capacity() noexcept {
    return _capacity;
}

template<class T, class I>
finite_element_matrix_1d<T, I>& nonstationary_heat_equation_solver_1d<T, I>::conductivity() noexcept {
    return _conductivity;
}

template<class T, class I>
const finite_element_matrix_1d<T, I>& nonstationary_heat_equation_solver_1d<T, I>::capacity() const noexcept {
    return _capacity;
}

template<class T, class I>
const finite_element_matrix_1d<T, I>& nonstationary_heat_equation_solver_1d<T, I>::conductivity() const noexcept {
    return _conductivity;
}

template<class T, class I>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& nonstationary_heat_equation_solver_1d<T, I>::temperature() const noexcept {
    return _temperature_curr;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& nonstationary_heat_equation_solver_1d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
const mesh::mesh_1d<T>& nonstationary_heat_equation_solver_1d<T, I>::mesh() const {
    return *mesh_ptr();
}

template<class T, class I>
const T nonstationary_heat_equation_solver_1d<T, I>::time_step() const noexcept {
    return _time_step;
}

template<class T, class I>
template<class Init_Dist>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const nonlocal::thermal::parameters_1d<T>& parameters,
                                                          const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                          const Init_Dist& init_dist) {
    const std::array<bool, 2> is_first = is_first_kind(boundaries_conditions);
    matrix_portrait_assembler_1d<T, I> capacity_portrait_assembler{_capacity, mesh_ptr()};
    matrix_portrait_assembler_1d<T, I> conductivity_portrait_assembler{_conductivity, mesh_ptr()};
    capacity_portrait_assembler.compute(std::vector<theory_t>(parameters.size(), theory_t::LOCAL), is_first);
    conductivity_portrait_assembler.compute(theories_types(parameters), is_first);

    heat_capacity_matrix_assembler_1d<T, I> capacity_assembler{_capacity, mesh_ptr()};
    thermal_conductivity_matrix_assembler_1d<T, I> conductivity_assembler{_conductivity, mesh_ptr()};
    capacity_assembler.compute(parameters, is_first);
    conductivity_assembler.compute(parameters, is_first);
    convection_condition_1d(_conductivity.inner(), boundaries_conditions);

    _conductivity.inner() *= time_step();
    _conductivity.inner() += _capacity.inner();
    reset_initial_values(_conductivity.inner(), {
        is_first.front() ? T{1} : _conductivity.inner().coeffRef(0, 0),
        is_first.back()  ? T{1} : _conductivity.inner().coeffRef(_conductivity.inner().rows() - 1, _conductivity.inner().cols() - 1)
    });
    
    for(std::unordered_map<size_t, T>& matrix_part : _conductivity.bound())
        for(auto& val : matrix_part | std::views::values)
            val *= time_step();

    _conductivity_initial_values = get_initial_values(_conductivity.inner());
    _capacity_initial_values = get_initial_values(_capacity.inner());

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
    reset_initial_values(_conductivity.inner(), _conductivity_initial_values);
    reset_initial_values(_capacity.inner(), _capacity_initial_values);

    radiation_condition_1d(_capacity.inner(), boundaries_conditions, time_step());
    boundary_condition_second_kind_1d<T>(_right_part, boundaries_conditions);
    if constexpr (!std::is_same_v<Right_Part, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        integrate_right_part(_right_part, mesh(), right_part);
    _right_part *= time_step();
    _right_part += _capacity.inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    
    radiation_condition_1d(_conductivity.inner(), boundaries_conditions, time_step());
    boundary_condition_first_kind_1d(_right_part, _conductivity.bound(), boundaries_conditions);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{_conductivity.inner()};
    _temperature_curr = solver.solveWithGuess(_right_part, _temperature_prev);
}

}