#pragma once

#include "thermal_parameters_1d.hpp"

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/base/equation_parameters.hpp>
#include <solvers/solver_1d/base/solution_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T>
class heat_equation_solution_1d : public solution_1d<T> {
    using _base = solution_1d<T>;

    std::vector<T> _temperature;
    std::vector<parameter_1d<T>> _parameters;
    std::vector<T> _flux;
    std::vector<metamath::types::vector_with_shifted_index<T>> _relaxation;
    T _time = T{ 0 };
    T _time_step = T{ 0 };

  public:
    using _base::mesh;

    template<class Parameter, class Vector>
    explicit heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const std::vector<Parameter>& parameters,
                                       const Vector& solution, const T time_step = T{ 0 });
    ~heat_equation_solution_1d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::vector<T>& flux() const;
    const parameter_1d<T>& parameter(const size_t segment) const noexcept;

    template<class Vector>
    void temperature(const Vector& solution);
    void time_step(const T step);
    void time(const T current_time);

    bool is_flux_calculated() const noexcept;
    const std::vector<T>& calc_flux();
};

template<class T>
coefficient_t<T, 1u> update_conductivity(const coefficient_t<T, 1u>& conductivity, const T relaxation_factor) {
    // clang-format off
    return std::visit(metamath::types::visitor{
        [relaxation_factor](const T value) -> coefficient_t<T, 1u> {
            return relaxation_factor * value;
        },
        [relaxation_factor](const spatial_dependency<T, 1u>& value) -> coefficient_t<T, 1u> { 
            return [value, relaxation_factor](const point<T, 1u>& x) { return relaxation_factor * value(x); };
        },
        [relaxation_factor](const solution_dependency<T, 1u>& value) -> coefficient_t<T, 1u> {
            return [value, relaxation_factor](const point<T, 1u>& x, const T temperature) { return relaxation_factor * value(x, temperature); };
        }
    }, conductivity);
    // clang-format on
}

template<class T>
template<class Parameter, class Vector>
heat_equation_solution_1d<T>::heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                        const std::vector<Parameter>& parameters, const Vector& solution,
                                                        const T time_step)
    : _base{ mesh, get_models(parameters) }, _temperature(solution.cbegin(), std::next(solution.cbegin(), mesh->nodes_count())),
      _parameters{ get_physical_parameters(parameters) }, _time_step{ time_step }, _relaxation(mesh->segments_count()) {
    for (const size_t segment : _base::mesh().segments())
        if (auto& param = _parameters[segment]; param.relaxation_time > T{ 0 }) {
            _relaxation[segment].container.resize(_base::mesh().qnodes_count(segment));
            _relaxation[segment].shift = _base::mesh().qnodes(segment).front();
            const T relaxation_factor = T{ 1 } - std::exp(-time_step / param.relaxation_time);
            param.conductivity = update_conductivity(param.conductivity, relaxation_factor);
        }
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::temperature() const noexcept {
    return _temperature;
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::flux() const {
    return _flux;
}

template<class T>
const parameter_1d<T>& heat_equation_solution_1d<T>::parameter(const size_t segment) const noexcept {
    return _parameters[segment];
}

template<class T>
template<class Vector>
void heat_equation_solution_1d<T>::temperature(const Vector& solution) {
    _temperature = std::vector<T>(solution.cbegin(), std::next(solution.cbegin(), _base::mesh().nodes_count()));
}

template<class T>
void heat_equation_solution_1d<T>::time(const T current_time) {
    _time = current_time;
}

template<class T>
bool heat_equation_solution_1d<T>::is_flux_calculated() const noexcept {
    return !_flux.empty();
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::calc_flux() {
    auto flux = this->weighted_gradient(temperature(), _parameters, &parameter_1d<T>::conductivity, T{ -1 });
    for (const size_t segment : mesh().segments())
        if (const T relax = _parameters[segment].relaxation_time; relax > T{ 0 }) {
            const T decay = std::exp(-_time / relax);
            for (const size_t q : mesh().qnodes(segment)) {
                const T current_flux = flux[q];
                flux[q] += decay * _relaxation[segment][q];
                _relaxation[segment][q] += decay * current_flux;
            }
        }
    _flux = mesh::utils::from_qnodes_to_nodes(mesh(), flux);
    return _flux;
}

} // namespace nonlocal::solver_1d::thermal
