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

template<class T>
using temperature_function = std::variant<
    std::function<T(const T)>,
    std::vector<T>,
    Eigen::Matrix<T, Eigen::Dynamic, 1>
>;

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
    parameters_1d<T> _parameters;
    const T _time_step = T{1};

    // relaxation model
    T _time = T{0};
    // The first half stores only values obtained from even segments, the second from odd segments.
    std::array<finite_element_matrix_1d<T, I>, 2> _relaxation_matrix;
    // The vector accumulating the relaxation term on the right part. 
    // The structure is similar, the first half is the values ​​obtained from even segments, the second from odd segments.
    std::array<Eigen::Matrix<T, Eigen::Dynamic, 1>, 2> _relaxation;

    static coefficient_t<T, 1u> update_conductivity(const coefficient_t<T, 1u>& conductivity, const T relaxation_factor);
    static std::array<T, 2> get_initial_values(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix);
    static void reset_initial_values(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                     const std::array<T, 2>& values);

    void accumulate_relaxation_vector();
    void add_relaxation_vector();
 
public:
    explicit nonstationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const parameters_1d<T>& parameters,
                                                   const T time_step, const T initial_time = T{0});

    T time() const noexcept;
    T time_step() const noexcept;
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;

    void initialize_temperature(const temperature_function<T>& temperature);

    void compute(const thermal_boundaries_conditions_1d<T>& boundaries_conditions);

    void calc_step(const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                   const std::optional<std::function<T(const T)>>& right_part = std::nullopt);
};

template<class T, class I>
nonstationary_heat_equation_solver_1d<T, I>::nonstationary_heat_equation_solver_1d(
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const parameters_1d<T>& parameters, const T time_step, const T initial_time)
    : _mesh{mesh}
    , _right_part{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_prev{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _temperature_curr{Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count())}
    , _parameters{parameters}
    , _time_step{time_step}
    , _time{initial_time} {
        for(const size_t segment : mesh->segments()) {
            auto& phys = _parameters[segment].physical;
            if (phys.relaxation_time > T{0}) {
                if (const size_t parity = segment % 2; _relaxation[parity].size() == 0)
                    _relaxation[parity] = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count());
                const T relaxation_factor = T{1} - std::exp(-time_step / phys.relaxation_time);
                phys.conductivity = update_conductivity(phys.conductivity, relaxation_factor);
            }
        }
    }

template<class T, class I>
coefficient_t<T, 1u> nonstationary_heat_equation_solver_1d<T, I>::update_conductivity(const coefficient_t<T, 1u>& conductivity, const T relaxation_factor) {
    return std::visit(metamath::visitor{
        [relaxation_factor](const T value) -> coefficient_t<T, 1u> { return relaxation_factor * value; },
        [relaxation_factor](const spatial_dependency<T, 1u>& value) -> coefficient_t<T, 1u> { 
            return [value, relaxation_factor](const point<T, 1u>& x) { return relaxation_factor * value(x); };
        },
        [relaxation_factor](const solution_dependency<T, 1u>& value) -> coefficient_t<T, 1u> {
            return [value, relaxation_factor](const point<T, 1u>& x, const T temperature) { return relaxation_factor * value(x, temperature); };
        }
    }, conductivity);
}

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
T nonstationary_heat_equation_solver_1d<T, I>::time() const noexcept {
    return _time;
}

template<class T, class I>
T nonstationary_heat_equation_solver_1d<T, I>::time_step() const noexcept {
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
void nonstationary_heat_equation_solver_1d<T, I>::initialize_temperature(const temperature_function<T>& temperature) {
    std::visit(metamath::visitor{
        [this](const std::function<T(const T)>& temperature) {
            for(const size_t i : std::ranges::iota_view{0u, mesh().nodes_count()})
                _temperature_curr[i] = temperature(mesh().node_coord(i));
        },
        [this](const auto& temperature) {
            if (temperature.size() != mesh().nodes_count())
                throw std::domain_error{"The initialization vector size does not match the number of nodes in the mesh."};
            for (const size_t i : std::ranges::iota_view{0u, mesh().nodes_count()})
                _temperature_curr[i] = temperature[i];
        }
    }, temperature);
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::compute(const thermal_boundaries_conditions_1d<T>& boundaries_conditions) {
    problem_settings settings = init_problem_settings(_parameters, boundaries_conditions, Is_Stationary);
    log_problem_settings(settings);

    init_matrix_portrait(_conductivity.inner, mesh(), settings);
    thermal_conductivity_assembler_1d<T, I> conductivity_assembler{_conductivity, _mesh};
    _conductivity.set_zero();
    conductivity_assembler.calc_matrix(_parameters, settings);
    convection_condition_1d(_conductivity.inner, boundaries_conditions);

    auto save_theories = settings.theories;
    settings.theories = std::vector<theory_t>(_parameters.size(), theory_t::LOCAL);
    init_matrix_portrait(_capacity.inner, mesh(), settings);
    _capacity.set_zero();
    heat_capacity_assembler_1d<T, I> capacity_assembler{_capacity, _mesh};
    capacity_assembler.calc_matrix(_parameters, settings.is_first_kind);
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

    // Initialize relaxation matrices
    settings.is_first_kind = {false, false};
    settings.theories = std::move(save_theories);
    for (const size_t parity : std::ranges::iota_view{0u, 2u}) {
        auto parameters = _parameters;
        bool need_initialization = false;
        for(const size_t segment : mesh().segments()) {
            auto& phys = parameters[segment].physical;
            const T relaxation_factor = phys.relaxation_time > T{0} && parity == segment % 2 ? std::exp(time() / phys.relaxation_time) : T{0};
            need_initialization = need_initialization || relaxation_factor != T{0};
            phys.conductivity = update_conductivity(phys.conductivity, relaxation_factor);
        }
        if (need_initialization) {
            init_matrix_portrait(_relaxation_matrix[parity].inner, mesh(), settings);
            thermal_conductivity_assembler_1d<T, I> assembler{_relaxation_matrix[parity], _mesh};
            _relaxation_matrix[parity].set_zero();
            assembler.calc_matrix(parameters, settings);
        }
    }
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::accumulate_relaxation_vector() {
    for(const size_t segment : mesh().segments())
        if (auto& phys = _parameters[segment].physical; phys.relaxation_time > T{0}) {
            const size_t parity = segment % 2;
            const T value = std::exp(time() / phys.relaxation_time);
            const auto& matrix = _relaxation_matrix[parity].inner;
            for(const size_t row : mesh().nodes(segment)) {
                const size_t ind = matrix.outerIndexPtr()[row];
                _relaxation[parity][row] += value * matrix.valuePtr()[ind] * _temperature_prev[matrix.innerIndexPtr()[ind]];
                for(const size_t i : std::ranges::iota_view{ind + 1, size_t(matrix.outerIndexPtr()[row + 1])}) {
                    const T factor = value * matrix.valuePtr()[i];
                    _relaxation[parity][row] += factor * _temperature_prev[matrix.innerIndexPtr()[i]];
                    _relaxation[parity][matrix.innerIndexPtr()[i]] += factor * _temperature_prev[row];
                }
            }
        }
}

template<class T, class I>
void nonstationary_heat_equation_solver_1d<T, I>::add_relaxation_vector() {
    for(const size_t segment : mesh().segments())
        if (auto& phys = _parameters[segment].physical; phys.relaxation_time > T{0}) {
            const size_t parity = segment % 2;
            const T value = std::exp(-time() / phys.relaxation_time );
            for(const size_t node : mesh().nodes(segment))
                _right_part[node] -= value * _relaxation[parity][node];
        }
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
    accumulate_relaxation_vector();
    _time += time_step();
    add_relaxation_vector();
    _right_part *= time_step();
    _right_part += _capacity.inner.template selfadjointView<Eigen::Upper>() * _temperature_prev;
    
    radiation_condition_1d<Is_Stationary>(_conductivity.inner, _right_part, boundaries_conditions, _temperature_prev, time_step());
    boundary_condition_first_kind_1d(_right_part, _conductivity.bound, boundaries_conditions);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{_conductivity.inner};
    _temperature_curr = solver.solveWithGuess(_right_part, _temperature_prev);
}

}