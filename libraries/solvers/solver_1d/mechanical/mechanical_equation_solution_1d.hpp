#pragma once

#include "mechanical_parameters_1d.hpp"

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/base/solution_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
class mechanical_equation_solution_1d : public solution_1d<T> {
    using _base = solution_1d<T>;

    const std::vector<T> _displacement;
    const std::vector<parameter_1d<T>> _parameters;
    std::vector<T> _stress;

  public:
    using _base::mesh;

    template<class Parameter, class Vector>
    explicit mechanical_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                             const std::vector<Parameter>& parameters, const Vector& solution);
    ~mechanical_equation_solution_1d() noexcept override = default;

    const std::vector<T>& displacement() const noexcept;
    const std::vector<T>& stress() const;
    const parameter_1d<T>& parameter(const size_t segment) const noexcept;

    bool is_stress_calculated() const noexcept;
    const std::vector<T>& calc_stress();
};

template<std::floating_point T>
template<class Parameter, class Vector>
mechanical_equation_solution_1d<T>::mechanical_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                    const std::vector<Parameter>& parameters,
                                                                    const Vector& solution)
    : _base{ mesh, get_models(parameters) }, _displacement(solution.cbegin(), std::next(solution.cbegin(), mesh->nodes_count())),
      _parameters{ get_physical_parameters(parameters) } {}

template<std::floating_point T>
const std::vector<T>& mechanical_equation_solution_1d<T>::displacement() const noexcept {
    return _displacement;
}

template<std::floating_point T>
const std::vector<T>& mechanical_equation_solution_1d<T>::stress() const {
    return _stress;
}

template<std::floating_point T>
const parameter_1d<T>& mechanical_equation_solution_1d<T>::parameter(const size_t segment) const noexcept {
    return _parameters[segment];
}

template<std::floating_point T>
bool mechanical_equation_solution_1d<T>::is_stress_calculated() const noexcept {
    return !_stress.empty();
}

template<std::floating_point T>
const std::vector<T>& mechanical_equation_solution_1d<T>::calc_stress() {
    if (!is_stress_calculated()) {
        const auto stress = this->weighted_gradient(displacement(), _parameters, &parameter_1d<T>::youngs_modulus, T{ 1 });
        _stress = mesh::utils::from_qnodes_to_nodes(mesh(), stress);
    }
    return _stress;
}

} // namespace nonlocal::solver_1d::mechanical
