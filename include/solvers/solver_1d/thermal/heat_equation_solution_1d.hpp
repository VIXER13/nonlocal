#ifndef NONLOCAL_HEAT_EQUATION_SOLUTION_1D_HPP
#define NONLOCAL_HEAT_EQUATION_SOLUTION_1D_HPP

#include "../solution_1d.hpp"
#include "mesh_1d_utils.hpp"

namespace nonlocal::thermal {

template<class T>
class heat_equation_solution_1d : public solution_1d<T> {
    using _base = solution_1d<T>;

    const std::vector<T> _temperature;
    const std::vector<T> _conductivity;
    std::optional<std::vector<T>> _flux;

    template<class Parameter>
    static std::vector<T> get_conductivity(const std::vector<Parameter>& parameters);
    
public:
    template<class Parameter, class Vector>
    explicit heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                       const std::vector<Parameter>& parameters,
                                       const Vector& solution);
    ~heat_equation_solution_1d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::vector<T>& flux() const;

    bool is_flux_calculated() const noexcept;
    const std::vector<T>& calc_flux();
};

template<class T>
template<class Parameter>
std::vector<T> heat_equation_solution_1d<T>::get_conductivity(const std::vector<Parameter>& parameters) {
    std::vector<T> conductivity(parameters.size());
    std::transform(parameters.cbegin(), parameters.cend(), conductivity.begin(),
                   [](const Parameter& parameter) { return parameter.physical.conductivity; } );
    return conductivity;
}

template<class T>
template<class Parameter, class Vector>
heat_equation_solution_1d<T>::heat_equation_solution_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                        const std::vector<Parameter>& parameters,
                                                        const Vector& solution)
    : _base{mesh, get_models(parameters)}
    , _temperature(solution.cbegin(), std::next(solution.cbegin(), mesh->nodes_count()))
    , _conductivity{get_conductivity(parameters)} {}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::temperature() const noexcept {
    return _temperature;
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::flux() const {
    return *_flux;
}

template<class T>
bool heat_equation_solution_1d<T>::is_flux_calculated() const noexcept {
    return _flux;
}

template<class T>
const std::vector<T>& heat_equation_solution_1d<T>::calc_flux() {
    const std::vector<T> gradient = mesh::utils::gradient_in_qnodes(_base::mesh(), temperature());
    _flux = mesh::utils::from_qnodes_to_nodes(_base::mesh(), gradient);
    std::vector<T> flux(_flux->size(), T{0});
    //auto& flux = *_flux;
    for(const size_t segment : _base::mesh().segments()) {
        const auto segment_nodes = _base::mesh().nodes(segment);
        for(const size_t node : segment_nodes)
           flux[node] -= _base::model(segment).local_weight * _conductivity[segment] * (*_flux)[node];

        if (theory_type(_base::model(segment).local_weight) == theory_t::NONLOCAL) {
            const auto& el = _base::mesh().element();
            const T jacobian = _base::mesh().jacobian(segment);
            std::vector<T> nonlocal_flux(segment_nodes.size(), T{0});
            for(const size_t eL : _base::mesh().elements(segment))
                for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()}) {
                    const size_t node = _base::mesh().node_number(eL, i);
                    const T node_coord = _base::mesh().node_coord(node);
                    for(const size_t eNL : _base::mesh().neighbours(eL)) {
                        const size_t qshift = eNL * el.qnodes_count();
                        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
                            const T influence_weight = _base::model(segment).influence(node_coord, _base::mesh().qnode_coord(eNL, q));
                            nonlocal_flux[node - segment_nodes.front()] += el.weight(q) * influence_weight * gradient[qshift + q] * jacobian;
                        }
                    }
                }
            const T nonloc_weight = nonlocal_weight(_base::model(segment).local_weight) * _conductivity[segment];
            for(const size_t node : segment_nodes) {
                flux[node] -= nonloc_weight * nonlocal_flux[node - segment_nodes.front()] / _base::mesh().node_elements(node).count();
            }
        }
    }
    _flux = flux;
    return *_flux;
}

}

#endif