#ifndef NONLOCAL_HEAT_EQUATION_SOLUTION_1D_HPP
#define NONLOCAL_HEAT_EQUATION_SOLUTION_1D_HPP

#include "../solution_1d.hpp"

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
    , _conductivity{} {}

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
    _flux = std::vector<T>{};
    return *_flux;
}

}

#endif