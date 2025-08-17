#pragma once

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/base/solution_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T>
class heat_equation_solution_1d : public solution_1d<T> {
    using _base = solution_1d<T>;

    const std::vector<T> _temperature;
    const std::vector<parameter_1d<T>> _parameters;
    std::vector<T> _flux;

    T evaluate(const coefficient_t<T, 1>& conductivity, const std::vector<T>& solution, const size_t e, const size_t q) const;

    void calc_local_flux();
    void calc_nonlocal_flux();
    
public:
    using _base::mesh;

    template<class Parameter, class Vector>
    explicit heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                       const std::vector<Parameter>& parameters,
                                       const Vector& solution);
    ~heat_equation_solution_1d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::vector<T>& flux() const;
    const parameter_1d<T>& parameter(const size_t segment) const noexcept;

    bool is_flux_calculated() const noexcept;
    const std::vector<T>& calc_flux();
    const std::vector<T>& calc_relaxation_flux(const std::vector<T>& relaxation_integral, const T time, const T relaxation_time);
};

template<class T>
template<class Parameter, class Vector>
heat_equation_solution_1d<T>::heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                        const std::vector<Parameter>& parameters,
                                                        const Vector& solution)
    : _base{mesh, get_models(parameters)}
    , _temperature(solution.cbegin(), std::next(solution.cbegin(), mesh->nodes_count()))
    , _parameters{get_physical_parameters(parameters)} {}

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
bool heat_equation_solution_1d<T>::is_flux_calculated() const noexcept {
    return !_flux.empty();
}

template<class T>
T heat_equation_solution_1d<T>::evaluate(const coefficient_t<T, 1>& conductivity, const std::vector<T>& solution, const size_t e, const size_t q) const {
    return std::visit(metamath::visitor{
        [](const T value) noexcept { return value; },
        [this, e, q](const spatial_dependency<T, 1u>& value) { return value(_base::mesh().qnode_coord(e, q)); },
        [this, &solution, e, q](const solution_dependency<T, 1u>& value) { 
            const size_t qshift = _base::mesh().qnode_number(e, q);
            return value(_base::mesh().qnode_coord(e, q), solution[qshift]); 
        }
    }, conductivity);
}

template<class T>
void heat_equation_solution_1d<T>::calc_local_flux() {
    const auto& el = mesh().element();
    _flux = mesh::utils::gradient_in_qnodes(mesh(), temperature());
    const bool is_any_nonlinear = 
        std::any_of(_parameters.begin(), _parameters.end(), [](const auto& parameter) constexpr noexcept {
            return std::holds_alternative<solution_dependency<T, 1>>(parameter.conductivity);
        });
    const auto temperature_in_qnodes = is_any_nonlinear ? mesh::utils::from_nodes_to_qnodes(mesh(), temperature()) : std::vector<T>{};
    for(const size_t segment : mesh().segments()) {
        const auto& param = parameter(segment);
        for(const size_t e : mesh().elements(segment)) {
            size_t qshift = e * el.qnodes_count();
            for(const size_t q : el.qnodes())
                _flux[qshift++] *= -evaluate(parameter(segment).conductivity, temperature_in_qnodes, e, q);
        }
    }
}

template<class T>
void heat_equation_solution_1d<T>::calc_nonlocal_flux() {
    const auto& el = mesh().element();
    const std::vector<T> flux = _flux;
    for(const size_t segment : mesh().segments())
        if (theory_type(_base::model(segment).local_weight) == theory_t::NONLOCAL) {
            for(const size_t eL : mesh().elements(segment)) {
                size_t qshiftL = eL * el.qnodes_count();
                for(const size_t qL : el.qnodes()) {
                    T nonlocal_flux = T{0};
                    const T qcoordL = mesh().qnode_coord(eL,  qL);
                    for(const size_t eNL : mesh().neighbours(eL)) {
                        size_t qshiftNL = eNL * el.qnodes_count();
                        for(const size_t qNL : el.qnodes()) {
                            const T qcoordNL = mesh().qnode_coord(eNL, qNL);
                            const T influence_weight = _base::model(segment).influence(qcoordL, qcoordNL);
                            nonlocal_flux += el.weight(qNL) * influence_weight * flux[qshiftNL++];
                        }
                    }
                    nonlocal_flux *= nonlocal_weight(_base::model(segment).local_weight) * mesh().jacobian(segment);
                    _flux[qshiftL] *= _base::model(segment).local_weight;
                    _flux[qshiftL++] += nonlocal_flux;
                }
            }
        }
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::calc_flux() {
    if (!is_flux_calculated()) {
        calc_local_flux();
        calc_nonlocal_flux();
        _flux = mesh::utils::from_qnodes_to_nodes(mesh(), _flux);
    }
    return _flux;
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::calc_relaxation_flux(
    const std::vector<T>& relaxation_integral, const T time, const T relaxation_time) {
    using namespace metamath::functions;
    calc_flux();
    _flux *= std::exp(-time / relaxation_time);
    _flux += relaxation_integral;
    return _flux;
}

}