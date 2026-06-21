#pragma once

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/base/solution_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
class mechanical_equation_solution_1d : public solution_1d<T> {
    using _base = solution_1d<T>;

    const std::vector<T> _displacement;
    const std::vector<parameter_1d<T>> _parameters;
    std::vector<T> _stress;

    
    T evaluate(const coefficient_t<T, 1>& youngs_modulus, const std::vector<T>& solution, const size_t e, const size_t q) const;

    void calc_local_stress();
    void calc_nonlocal_stress();

public:
    using _base::mesh;

    template<class Parameter, class Vector>
    explicit mechanical_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                             const std::vector<Parameter>& parameters,
                                             const Vector& solution);
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
    : _base{mesh, get_models(parameters)}
    , _displacement(solution.cbegin(), std::next(solution.cbegin(), mesh->nodes_count()))
    , _parameters{get_physical_parameters(parameters)} {}

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
T mechanical_equation_solution_1d<T>::evaluate(const coefficient_t<T, 1>& youngs_modulus, const std::vector<T>& solution, const size_t e, const size_t q) const {
    return std::visit(metamath::visitor{
        [](const T value) noexcept { return value; },
        [this, e, q](const spatial_dependency<T, 1u>& value) { return value(_base::mesh().qnode_coord(e, q)); },
        [this, &solution, e, q](const solution_dependency<T, 1u>& value) { 
            const size_t qshift = _base::mesh().qnode_number(e, q);
            return value(_base::mesh().qnode_coord(e, q), solution[qshift]); 
        }
    }, youngs_modulus);
}

template<std::floating_point T>
void mechanical_equation_solution_1d<T>::calc_local_stress() {
    const auto& el = mesh().element();
    _stress = mesh::utils::gradient_in_qnodes(mesh(), displacement());
    const bool is_any_nonlinear = 
        std::any_of(_parameters.begin(), _parameters.end(), [](const auto& parameter) constexpr noexcept {
            return std::holds_alternative<solution_dependency<T, 1>>(parameter.youngs_modulus);
        });
    const auto displacement_in_qnodes = is_any_nonlinear ? mesh::utils::from_nodes_to_qnodes(mesh(), displacement()) : std::vector<T>{};
    for(const size_t segment : mesh().segments()) {
        const auto& param = parameter(segment);
        for(const size_t e : mesh().elements(segment)) {
            size_t qshift = e * el.qnodes_count();
            for(const size_t q : el.qnodes())
                _stress[qshift++] *= -evaluate(parameter(segment).youngs_modulus, displacement_in_qnodes, e, q);
        }
    }
}

template<std::floating_point T>
void mechanical_equation_solution_1d<T>::calc_nonlocal_stress() {
    const auto& el = mesh().element();
    const std::vector<T> stress = _stress;
    for(const size_t segment : mesh().segments())
        if (theory_type(_base::model(segment).local_weight) == theory_t::NONLOCAL) {
            for(const size_t eL : mesh().elements(segment)) {
                size_t qshiftL = eL * el.qnodes_count();
                for(const size_t qL : el.qnodes()) {
                    T nonlocal_stress = T{0};
                    const T qcoordL = mesh().qnode_coord(eL,  qL);
                    for(const size_t eNL : mesh().neighbours(eL)) {
                        size_t qshiftNL = eNL * el.qnodes_count();
                        for(const size_t qNL : el.qnodes()) {
                            const T qcoordNL = mesh().qnode_coord(eNL, qNL);
                            const T influence_weight = _base::model(segment).influence(qcoordL, qcoordNL);
                            nonlocal_stress += el.weight(qNL) * influence_weight * stress[qshiftNL++];
                        }
                    }
                    nonlocal_stress *= nonlocal_weight(_base::model(segment).local_weight) * mesh().jacobian(segment);
                    _stress[qshiftL] *= _base::model(segment).local_weight;
                    _stress[qshiftL++] += nonlocal_stress;
                }
            }
        }
}

template<std::floating_point T>
const std::vector<T>& mechanical_equation_solution_1d<T>::calc_stress() {
    if (!is_stress_calculated()) {
        calc_local_stress();
        calc_nonlocal_stress();
        _stress = mesh::utils::from_qnodes_to_nodes(mesh(), _stress);
        for (auto& x: _stress) x *= static_cast<T>(-1);
    }
    return _stress;
}


}