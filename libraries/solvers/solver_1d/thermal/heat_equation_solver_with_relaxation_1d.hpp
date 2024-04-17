#pragma once

#include "nonlocal_constants.hpp"

#include "matrix_portrait_assembler_1d.hpp"
#include "thermal_conductivity_matrix_assembler_1d.hpp"
#include "heat_capacity_matrix_assembler_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class nonstationary_heat_equation_solver_1d final {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    finite_element_matrix_1d<T, I> _capacity;
    finite_element_matrix_1d<T, I> _conductivity;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_prev;
    Eigen::Matrix<T, Eigen::Dynamic, 1> _temperature_curr;
    parameters_1d<T> _parameters;
    const T _time_step = T{1}; // The solver is designed for a constant time step.
    T _time = T{0}; // We assume that initial time is zero, but if you want to change it, 
                    // please update compute function according to this changes.
    // Relaxation problems parameters
    std::array<finite_element_matrix_1d<T, I>, 2> _relaxation_matrix; // The first half stores only values obtained from even segments, 
                                                                      // the second from odd segments.
    std::array<Eigen::Matrix<T, Eigen::Dynamic, 1>, 2> _relaxation; // The vector accumulating the relaxation term on the right part. 
                                                                    // The structure is similar, the first half is the values ​​obtained from even segments,
                                                                    // the second from odd segments.

    void assemble_portraits(const std::array<bool, 2> is_temperature);
    void assemble_initial_matrices(const std::array<bool, 2> is_temperature);
    void fix_boundaries_values(const std::array<bool, 2> is_temperature);

    template<class Right_Part>
    void calculate_right_part(const Right_Part& right_part);
    void accumulate_relaxation_vector();
    void add_relaxation_vector();
 
public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, 
                                                   const parameters_1d<T>& parameters,
                                                   const T time_step);

    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    const T time_step() const noexcept;
    const T time() const noexcept;

    void compute(const thermal_boundaries_conditions_1d<T>& boundaries_conditions);
    template<class Initial_Distribution>
    void initialize_temperature(const Initial_Distribution& initial_distribution);
    template<class Right_Part>
    void calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                   const Right_Part& right_part);
};

template<class T, class I>
nonstationary_heat_equation_solver_1d<T, I>::nonstationary_heat_equation_solver_1d(
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh, 
    const parameters_1d<T>& parameters,
    const T time_step)
    : _mesh{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _parameters{parameters}
    , _time_step{time_step} {
    if (_parameters.size() != mesh->segments_count())
        throw std::domain_error{"The parameters number does not match the segments number in the mesh."};
    for(const size_t segment : mesh->segments()) {
        {   // Copy physical parameter, because parameters are shared_ptr. 
            // TODO: need rework parameters structure
            auto& physical = parameter_cast<coefficients_t::CONSTANTS>(*_parameters[segment].physical);
            _parameters[segment].physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(physical);
        }
        if (auto& physical = parameter_cast<coefficients_t::CONSTANTS>(*_parameters[segment].physical);
            physical.relaxation_time > T{0}) {
            if (const size_t parity = segment % 2; _relaxation[parity].size() == 0)
                _relaxation[parity] = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count());
            physical.conductivity *= T{1} - std::exp(-time_step / physical.relaxation_time);
        }
    }
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
const T nonstationary_heat_equation_solver_1d<T, I>::time() const noexcept {
    return _time;
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::assemble_portraits(const std::array<bool, 2> is_temperature) {
    const std::vector<theory_t> theories = theories_types(_parameters);
    matrix_portrait_assembler_1d<T, I> capacity_portrait_assembler{_capacity, mesh_ptr()};
    matrix_portrait_assembler_1d<T, I> conductivity_portrait_assembler{_conductivity, mesh_ptr()};
    capacity_portrait_assembler.compute(std::vector<theory_t>(_parameters.size(), theory_t::LOCAL), is_temperature);
    conductivity_portrait_assembler.compute(theories, is_temperature);

    std::array<std::vector<size_t>, 2> initialize_nodes;
    for(const size_t segment : mesh().segments())
        if (const size_t parity = segment % 2; _parameters[segment].physical->relaxation_time > T{0})
            for(const size_t node : mesh().nodes(segment))
                initialize_nodes[parity].push_back(node);
    for(const size_t parity : std::ranges::iota_view{0u, 2u}) {
        matrix_portrait_assembler_1d<T, I> relaxation_portrait_assembler{_relaxation_matrix[parity], mesh_ptr(), initialize_nodes[parity]};
        relaxation_portrait_assembler.compute(theories, {false, false});
    }
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::assemble_initial_matrices(const std::array<bool, 2> is_temperature) {
    heat_capacity_matrix_assembler_1d<T, I> capacity_assembler{_capacity, mesh_ptr()};
    thermal_conductivity_matrix_assembler_1d<T, I> conductivity_assembler{_conductivity, mesh_ptr()};
    capacity_assembler.compute(_parameters, is_temperature);
    conductivity_assembler.compute(_parameters, is_temperature);

    std::array<std::vector<size_t>, 2> initialize_nodes;
    for(const size_t segment : mesh().segments())
        if (const size_t parity = segment % 2; _parameters[segment].physical->relaxation_time > T{0})
            for(const size_t node : mesh().nodes(segment))
                initialize_nodes[parity].push_back(node);
    for(const size_t parity : std::ranges::iota_view{0u, 2u}) {
        auto parameters = _parameters;
        for(const size_t segment : mesh().segments()) {
            {   // Copy physical parameter, because parameters are shared_ptr. 
                // TODO: need rework parameters structure
                auto& physical = parameter_cast<coefficients_t::CONSTANTS>(*parameters[segment].physical);
                parameters[segment].physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(physical);
            }
            auto& physical = parameter_cast<coefficients_t::CONSTANTS>(*parameters[segment].physical);
            physical.conductivity *= physical.relaxation_time > T{0} && (parity % 2 == segment % 2) ? std::exp(time() / physical.relaxation_time) : T{0};
        }
        thermal_conductivity_matrix_assembler_1d<T, I> relaxation_assembler{_relaxation_matrix[parity], mesh_ptr(), initialize_nodes[parity]};
        relaxation_assembler.compute(parameters, {false, false});
    }
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::fix_boundaries_values(const std::array<bool, 2> is_temperature) {
    if (is_temperature.front())
        _conductivity.inner().valuePtr()[0] = T{1};
    if (is_temperature.back())
        _conductivity.inner().valuePtr()[_conductivity.inner().nonZeros() - 1] = T{1};
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const thermal_boundaries_conditions_1d<T>& boundaries_conditions) {
    const std::array<bool, 2> is_temperature = is_first_kind(boundaries_conditions);
    assemble_portraits(is_temperature);
    assemble_initial_matrices(is_temperature);
    convection_condition_1d(_conductivity.inner(), boundaries_conditions);

    _conductivity.inner() *= time_step();
    for(std::unordered_map<size_t, T>& matrix_part : _conductivity.bound())
        for(auto& val : matrix_part | std::views::values)
            val *= time_step();

    _conductivity.inner() += _capacity.inner();
    fix_boundaries_values(is_temperature);
}

template<class T, class I>
template<class Initial_Distribution>
void nonstationary_heat_equation_solver_1d<T, I>::initialize_temperature(const Initial_Distribution& initial_distribution) {
    if constexpr (std::is_same_v<Initial_Distribution, std::vector<T>> ||
                  std::is_same_v<Initial_Distribution, Eigen::Matrix<T, Eigen::Dynamic, 1>>) {
        if (initial_distribution.size() != temperature().size())
            throw std::domain_error{"The temperature initializing vector size is not equal to the temperature vector size."};
        for(const size_t i : std::ranges::iota_view{0u, initial_distribution.size()})
            _temperature_curr[i] = initial_distribution[i];
    } else if constexpr (!std::is_same_v<Initial_Distribution, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        for(const size_t i : std::ranges::iota_view{0u, mesh().nodes_count()})
            _temperature_curr[i] = initial_distribution(mesh().node_coord(i));
}

template<class T, class I>
template<class Right_Part>
void nonstationary_heat_equation_solver_1d<T, I>::calculate_right_part(const Right_Part& right_part) {
    if constexpr (std::is_same_v<Right_Part, std::vector<T>> ||
                  std::is_same_v<Right_Part, Eigen::Matrix<T, Eigen::Dynamic, 1>>) {
        if (right_part.size() != _right_part().size())
            throw std::domain_error{"The right part vector size is not equal to the solver's right part vector size."};
        for(const size_t i : std::ranges::iota_view{0u, right_part.size()})
            _right_part[i] = right_part[i];
    }
    else if constexpr (!std::is_same_v<Right_Part, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        integrate_right_part(_right_part, mesh(), right_part);
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::accumulate_relaxation_vector() {
    for(const size_t segment : mesh().segments())
        if (const size_t parity = segment % 2; _parameters[segment].physical->relaxation_time > T{0}) {
            const T value = std::exp(time() / _parameters[segment].physical->relaxation_time);
            auto& matrix = _relaxation_matrix[parity].inner();
            for(const size_t row : mesh().nodes(segment)) {
                const size_t ind = matrix.outerIndexPtr()[row];
                _relaxation[parity][row] += value * matrix.valuePtr()[ind] * _temperature_prev[matrix.innerIndexPtr()[ind]];
                for(const size_t i : std::ranges::iota_view{ind + 1, size_t(matrix.outerIndexPtr()[row + 1])}) {
                    _relaxation[parity][row] += value * matrix.valuePtr()[i] * _temperature_prev[matrix.innerIndexPtr()[i]];
                    _relaxation[parity][matrix.innerIndexPtr()[i]] += value * matrix.valuePtr()[i] * _temperature_prev[row];
                }
            }
        }
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::add_relaxation_vector() {
    for(const size_t segment : mesh().segments())
        if (const size_t parity = segment % 2; _parameters[segment].physical->relaxation_time > T{0}) {
            const T value = std::exp(-time() / _parameters[segment].physical->relaxation_time);
            for(const size_t node : mesh().nodes(segment))
                _right_part[node] -= value * _relaxation[parity][node];
        }
}

template<class T, class I>
template<class Right_Part>
void nonstationary_heat_equation_solver_1d<T, I>::calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                            const Right_Part& right_part) {
    _temperature_prev.swap(_temperature_curr);
    _right_part.setZero();
    boundary_condition_second_kind_1d<T>(_right_part, boundaries_conditions);
    calculate_right_part(right_part);
    accumulate_relaxation_vector();
    _time += time_step();
    add_relaxation_vector();
    _right_part *= time_step();
    _right_part += _capacity.inner().template selfadjointView<Eigen::Upper>() * _temperature_prev;
    boundary_condition_first_kind_1d(_right_part, _conductivity.bound(), boundaries_conditions);
    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{_conductivity.inner()};
    _temperature_curr = solver.solveWithGuess(_right_part, _temperature_prev);
}

}