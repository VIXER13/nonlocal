#pragma once

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
    std::array<std::vector<T>, 2> _relaxation;
    T _time = T{0};
    T _time_step = T{0};

    T evaluate(const coefficient_t<T, 1>& conductivity, const std::vector<T>& solution, const size_t e, const size_t q) const;

    std::array<std::vector<T>, 2> calc_local_flux();
    std::vector<T> calc_nonlocal_flux(const std::vector<T>& local_flux);
    
public:
    using _base::mesh;

    template<class Parameter, class Vector>
    explicit heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                       const std::vector<Parameter>& parameters,
                                       const Vector& solution,
                                       const T time_step = T{0});
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
    return std::visit(metamath::visitor{
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
}

template<class T>
template<class Parameter, class Vector>
heat_equation_solution_1d<T>::heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                        const std::vector<Parameter>& parameters,
                                                        const Vector& solution,
                                                        const T time_step)
    : _base{mesh, get_models(parameters)}
    , _temperature(solution.cbegin(), std::next(solution.cbegin(), mesh->nodes_count()))
    , _parameters{get_physical_parameters(parameters)}
    , _time_step{time_step} {
        _relaxation.front().resize(mesh->qnodes_count(), T{0});
        _relaxation.back().resize(mesh->qnodes_count(), T{0});
        for(auto& param : _parameters)
            if (param.relaxation_time > T{0}) {
                const T relaxation_factor = T{1} - std::exp(-time_step / param.relaxation_time);
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
std::array<std::vector<T>, 2> heat_equation_solution_1d<T>::calc_local_flux() {
    const auto& el = mesh().element();
    std::vector<T> gradient = mesh::utils::gradient_in_qnodes(mesh(), temperature());
    const bool is_any_nonlinear = 
        std::any_of(_parameters.begin(), _parameters.end(), [](const auto& parameter) constexpr noexcept {
            return std::holds_alternative<solution_dependency<T, 1>>(parameter.conductivity);
        });
    const auto temperature_in_qnodes = is_any_nonlinear ? mesh::utils::from_nodes_to_qnodes(mesh(), temperature()) : std::vector<T>{};
    std::array<std::vector<T>, 2> flux = {std::vector<T>(gradient.size(), T{0}), std::vector<T>(gradient.size(), T{0})};
    for(const size_t segment : mesh().segments()) {
        const size_t parity = segment % 2;
        const auto& param = parameter(segment);
        for(const size_t e : mesh().elements(segment)) {
            size_t qshift = e * el.qnodes_count();
            for(const size_t q : el.qnodes()) {
                flux[parity][qshift] -= gradient[qshift] * evaluate(parameter(segment).conductivity, temperature_in_qnodes, e, q);
                ++qshift;
            }
        }
    }
    return flux;
}

template<class T>
std::vector<T> heat_equation_solution_1d<T>::calc_nonlocal_flux(const std::vector<T>& local_flux) {
    const auto& el = mesh().element();
    std::vector<T> nonlocal_flux = local_flux;
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
                            nonlocal_flux += el.weight(qNL) * influence_weight * local_flux[qshiftNL++];
                        }
                    }
                    nonlocal_flux *= nonlocal_weight(_base::model(segment).local_weight) * mesh().jacobian(segment);
                    nonlocal_flux[qshiftL] *= _base::model(segment).local_weight;
                    nonlocal_flux[qshiftL++] += nonlocal_flux;
                }
            }
        }
    return nonlocal_flux;
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::calc_flux() {
    //if (!is_flux_calculated()) {
    using namespace metamath::functions;
    std::vector<T> flux(_base::mesh().qnodes_count(), T{0});
    std::array<std::vector<T>, 2> local_flux = calc_local_flux();
    for(const size_t q : std::ranges::iota_view{0zu, flux.size()})
        flux[q] = local_flux[0][q] + local_flux[1][q];

    for(const size_t segment : _base::mesh().segments())
        if (const T relax = _parameters[segment].relaxation_time; relax > T{0}) {
            const size_t parity = segment % 2;
            for(const size_t q : _base::mesh().qnodes(segment))
                flux[q] += std::exp(-_time / relax) * _relaxation[parity][q];
        }

    _flux = mesh::utils::from_qnodes_to_nodes(mesh(), flux);

    for(const size_t segment : _base::mesh().segments())
        if (const T relax = _parameters[segment].relaxation_time; relax > T{0}) {
            const size_t parity = segment % 2;
            for(const size_t q : _base::mesh().qnodes(segment))
                _relaxation[parity][q] += std::exp(_time / relax) * local_flux[parity][q];
        }

    //}
    return _flux;
}

}