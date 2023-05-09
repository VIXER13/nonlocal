#ifndef NONLOCAL_HEAT_EQUATION_SOLUTION_1D_HPP
#define NONLOCAL_HEAT_EQUATION_SOLUTION_1D_HPP

#include "solution_1d.hpp"
#include "mesh_1d_utils.hpp"

namespace nonlocal::thermal {

template<class T>
class heat_equation_solution_1d : public solution_1d<T> {
    using _base = solution_1d<T>;

    const std::vector<T> _temperature;
    const std::vector<parameter_1d_sptr<T>> _parameters;
    std::optional<std::vector<T>> _flux;
    
public:
    using _base::mesh;

    template<class Parameter, class Vector>
    explicit heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                       const std::vector<Parameter>& parameters,
                                       const Vector& solution);
    ~heat_equation_solution_1d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::vector<T>& flux() const;
    const parameter_1d_sptr<T>& parameter(const size_t segment) const noexcept;

    bool is_flux_calculated() const noexcept;
    const std::vector<T>& calc_flux();
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
    return *_flux;
}

template<class T>
const parameter_1d_sptr<T>& heat_equation_solution_1d<T>::parameter(const size_t segment) const noexcept {
    return _parameters[segment];
}

template<class T>
bool heat_equation_solution_1d<T>::is_flux_calculated() const noexcept {
    return _flux.has_value();
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::calc_flux() {
    if (is_flux_calculated())
        return *_flux;

    const auto& el = mesh().element();
    const bool is_nonlinear = std::any_of(_parameters.begin(), _parameters.end(), [](const auto& parameter) noexcept {
        return parameter->type == coefficients_t::SOLUTION_DEPENDENT;
    });
    const auto temperature_in_qnodes = is_nonlinear ? mesh::utils::from_nodes_to_qnodes(mesh(), temperature()) : std::vector<T>{};
    std::vector<T> flux = mesh::utils::gradient_in_qnodes(mesh(), temperature());
    std::vector<T> flux_nonlocal;

    for(const size_t segment : mesh().segments()) {
        const auto& param = *parameter(segment);
        const auto segment_elements = mesh().elements(segment);
        for(const size_t e : segment_elements) {
            size_t qshift = e * el.qnodes_count();
            for(const size_t q : el.qnodes()) {
                using enum coefficients_t;
                const T conductivity =
                    param.type == CONSTANTS ?
                    parameter_cast<CONSTANTS>(param).conductivity :
                    param.type == SPACE_DEPENDENT ?
                    parameter_cast<SPACE_DEPENDENT>(param).conductivity(mesh().qnode_coord(e, q)) :
                    param.type == SOLUTION_DEPENDENT ?
                    parameter_cast<SOLUTION_DEPENDENT>(param).conductivity(mesh().qnode_coord(e, q), temperature_in_qnodes[qshift]) :
                    throw std::domain_error{"Unknown parameter type"};
                flux[qshift] *= -conductivity;
                ++qshift;
            }
        }

        if (theory_type(_base::model(segment).local_weight) == theory_t::NONLOCAL) {
            flux_nonlocal.clear();
            flux_nonlocal.resize(segment_elements.size() * el.qnodes_count(), T{0});
            for(const size_t eL : segment_elements) {
                size_t qshiftL = eL * el.qnodes_count();
                for(const size_t qL : el.qnodes()) {
                    const T qcoordL = mesh().qnode_coord(eL,  qL);
                    for(const size_t eNL : mesh().neighbours(eL)) {
                        size_t qshiftNL = eNL * el.qnodes_count();
                        for(const size_t qNL : el.qnodes()) {
                            const T qcoordNL = mesh().qnode_coord(eNL, qNL);
                            const T influence_weight = _base::model(segment).influence(qcoordL, qcoordNL);
                            flux_nonlocal[qshiftL] += el.weight(qNL) * influence_weight * flux[qshiftNL++];
                        }
                    }
                    ++qshiftL;
                }
            }
            using namespace metamath::functions;
            flux *= _base::model(segment).local_weight;
            flux_nonlocal *= nonlocal_weight(_base::model(segment).local_weight) * mesh().jacobian(segment);
            size_t qshift = segment_elements.front() * el.qnodes_count();
            for(const T value : flux_nonlocal)
                flux[qshift++] += value;
        }
    }

    _flux = mesh::utils::from_qnodes_to_nodes(mesh(), flux);
    return *_flux;
}

}

#endif