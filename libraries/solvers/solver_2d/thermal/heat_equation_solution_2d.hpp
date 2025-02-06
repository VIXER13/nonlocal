#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/base/solution_2d.hpp>

namespace nonlocal::thermal {

template<class T, class I = uint32_t>
class heat_equation_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;

    std::vector<T> _temperature;
    std::array<std::vector<T>, 2> _flux;
    std::unordered_map<std::string, parameter_2d_sptr<T>> _parameters;

    std::array<std::vector<T>, 2> local_flux_in_qnodes() const;

public:
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    template<class Vector>
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                       const parameters_2d<T>& parameters, const Vector& temperature);
    ~heat_equation_solution_2d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::array<std::vector<T>, 2>& flux() const;
    const parameter_2d<T>& parameter(const std::string& group) const;

    T calc_energy() const;
    bool is_flux_calculated() const noexcept;
    const std::array<std::vector<T>, 2>& calc_flux();
};

template<class T, class I>
heat_equation_solution_2d<T, I>::heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh}
    , _temperature(mesh->container().nodes_count(), T{0}) {}

template<class T, class I>
template<class Vector>
heat_equation_solution_2d<T, I>::heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                           const parameters_2d<T>& parameters, const Vector& temperature)
    : _base{mesh, get_models(parameters)}
    , _temperature{temperature.cbegin(), std::next(temperature.cbegin(), mesh->container().nodes_count())}
    , _parameters{get_physical_parameters(parameters)} {}

template<class T, class I>
const std::vector<T>& heat_equation_solution_2d<T, I>::temperature() const noexcept {
    return _temperature;
}

template<class T, class I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::flux() const {
    if (!is_flux_calculated())
        throw std::runtime_error{"Flux wasn't calculated"};
    return _flux;
}

template<class T, class I>
const parameter_2d<T>& heat_equation_solution_2d<T, I>::parameter(const std::string& group) const {
    return _parameters.at(group);
}

template<class T, class I>
T heat_equation_solution_2d<T, I>::calc_energy() const {
    //return mesh::integrate(*_base::mesh_proxy(), temperature());
    return 0;
}

template<class T, class I>
bool heat_equation_solution_2d<T, I>::is_flux_calculated() const noexcept {
    return !_flux[X].empty() && !_flux[Y].empty();
}

template<class T, class I>
std::array<std::vector<T>, 2> heat_equation_solution_2d<T, I>::local_flux_in_qnodes() const {
    auto flux = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
    for(const auto& [group, parameter] : _parameters)
        for(const size_t e : _base::mesh().container().elements(group))
            for(const size_t qshift : _base::mesh().quad_shifts_count(e)) {
                const std::array<T, 2>& qnode = _base::mesh().quad_coord(qshift);
                using enum coefficients_t;
                switch (parameter->material) {
                case material_t::ISOTROPIC: {
                    const T conductivity =
                        parameter->type == CONSTANTS ?
                        parameter_cast<CONSTANTS>(*parameter).conductivity[X][X] :
                        parameter->type == SPACE_DEPENDENT ?
                        parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[X][X](qnode) :
                        parameter->type == SOLUTION_DEPENDENT ?
                        parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[X][X](qnode, 0.0) :
                        throw std::domain_error{"Unknown parameter type"};
                    flux[X][qshift] *= -conductivity;
                    flux[Y][qshift] *= -conductivity;
                } break;

                case material_t::ORTHOTROPIC: {
                    using U = std::array<T, 2>;
                    const U conductivity =
                        parameter->type == CONSTANTS ?
                        U{parameter_cast<CONSTANTS>(*parameter).conductivity[X][X], 
                          parameter_cast<CONSTANTS>(*parameter).conductivity[Y][Y]} :
                        parameter->type == SPACE_DEPENDENT ?
                        U{parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[X][X](qnode),
                          parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[Y][Y](qnode)} :
                        parameter->type == SOLUTION_DEPENDENT ?
                        U{parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[X][X](qnode, 0.0),
                          parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[Y][Y](qnode, 0.0)} :
                        throw std::domain_error{"Unknown parameter type"};
                    flux[X][qshift] *= -conductivity[X];
                    flux[Y][qshift] *= -conductivity[X];
                } break;

                case material_t::ANISOTROPIC: {
                    using U = metamath::types::square_matrix<T, 2>;
                    const U conductivity =
                        parameter->type == CONSTANTS ?
                        U{parameter_cast<CONSTANTS>(*parameter).conductivity[X][X],
                          parameter_cast<CONSTANTS>(*parameter).conductivity[X][Y],
                          parameter_cast<CONSTANTS>(*parameter).conductivity[Y][X],
                          parameter_cast<CONSTANTS>(*parameter).conductivity[Y][Y]} :
                        parameter->type == SPACE_DEPENDENT ?
                        U{parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[X][X](qnode),
                          parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[X][Y](qnode),
                          parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[Y][X](qnode),
                          parameter_cast<SPACE_DEPENDENT>(*parameter).conductivity[Y][Y](qnode)} :
                        parameter->type == SOLUTION_DEPENDENT ?
                        U{parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[X][X](qnode, 0.0),
                          parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[X][Y](qnode, 0.0),
                          parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[Y][X](qnode, 0.0),
                          parameter_cast<SOLUTION_DEPENDENT>(*parameter).conductivity[Y][Y](qnode, 0.0)} :
                        throw std::domain_error{"Unknown parameter type"};
                    const std::array<T, 2> temp = {
                        -(conductivity[X][X] * flux[X][qshift] + conductivity[X][Y] * flux[Y][qshift]),
                        -(conductivity[Y][X] * flux[X][qshift] + conductivity[Y][Y] * flux[Y][qshift])
                    };
                    flux[X][qshift] = temp[X];
                    flux[Y][qshift] = temp[Y];
                } break;

                default:
                    throw std::runtime_error{"Unknown material type."};
                }
            }
    return flux;
}

template<class T, class I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::calc_flux() {
    if (is_flux_calculated())
        return _flux;

    _flux = local_flux_in_qnodes();
    std::array<std::vector<T>, 2> flux = _flux;
    for(const auto& [group, parameter] : _parameters)
        if (const model_parameters<2, T>& model = _base::model(group); theory_type(model.local_weight) == theory_t::NONLOCAL) {
            const T nonlocal_weight = nonlocal::nonlocal_weight(model.local_weight);
            for(const size_t eL : _base::mesh().container().elements(group))
                for(const size_t qshiftL : _base::mesh().quad_shifts_count(eL)) {
                    std::array<T, 2> nonlocal_gradient = {};
                    const auto& qcoordL = _base::mesh().quad_coord(qshiftL);
                    for(const size_t eNL : _base::mesh().neighbours(eL)) {
                        size_t qshiftNL = _base::mesh().quad_shift(eNL);
                        const auto& elNL = _base::mesh().container().element_2d(eNL);
                        for(const size_t qNL : elNL.qnodes()) {
                            const T influence_weight = elNL.weight(qNL) * mesh::jacobian(_base::mesh().jacobi_matrix(qshiftNL)) *
                                                       model.influence(qcoordL, _base::mesh().quad_coord(qshiftNL));
                            nonlocal_gradient[X] += influence_weight * _flux[X][qshiftNL];
                            nonlocal_gradient[Y] += influence_weight * _flux[Y][qshiftNL];
                            ++qshiftNL;
                        }
                    }
                    using namespace metamath::functions;
                    nonlocal_gradient *= nonlocal_weight;
                    flux[X][qshiftL] *= model.local_weight;
                    flux[Y][qshiftL] *= model.local_weight;
                    flux[X][qshiftL] += nonlocal_gradient[X];
                    flux[Y][qshiftL] += nonlocal_gradient[Y];
                }
        }
    _flux[X] = mesh::utils::qnodes_to_nodes(_base::mesh(), flux[X]);
    _flux[Y] = mesh::utils::qnodes_to_nodes(_base::mesh(), flux[Y]);
    _flux[X] = parallel::all_to_all(_flux[X], _base::mesh().MPI_ranges());
    _flux[Y] = parallel::all_to_all(_flux[Y], _base::mesh().MPI_ranges());
    return _flux;
}

}