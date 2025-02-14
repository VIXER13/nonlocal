#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/base/solution_2d.hpp>

namespace nonlocal::thermal {

template<std::floating_point T, std::integral I = uint32_t>
class heat_equation_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;

    std::vector<T> _solution; // stub for nonlinear problems

    std::vector<T> _temperature;
    std::array<std::vector<T>, 2> _flux;
    std::unordered_map<std::string, parameter_2d<T>> _parameters;

    std::array<std::vector<T>, 2> local_flux_in_qnodes() const;

    T evaluate(const isotropic_conductivity_t<T>& conductivity, const size_t qshift) const;
    std::array<T, 2> evaluate(const orthotropic_conductivity_t<T>& conductivity, const size_t qshift) const;
    metamath::types::square_matrix<T, 2> evaluate(const anisotropic_conductivity_t<T>& conductivity, const size_t qshift) const;

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

template<std::floating_point T, std::integral I>
heat_equation_solution_2d<T, I>::heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh}
    , _temperature(mesh->container().nodes_count(), T{0}) {}

template<std::floating_point T, std::integral I>
template<class Vector>
heat_equation_solution_2d<T, I>::heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                           const parameters_2d<T>& parameters, const Vector& temperature)
    : _base{mesh, get_models(parameters)}
    , _temperature{temperature.cbegin(), std::next(temperature.cbegin(), mesh->container().nodes_count())}
    , _parameters{get_physical_parameters(parameters)} {}

template<std::floating_point T, std::integral I>
const std::vector<T>& heat_equation_solution_2d<T, I>::temperature() const noexcept {
    return _temperature;
}

template<std::floating_point T, std::integral I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::flux() const {
    if (!is_flux_calculated())
        throw std::runtime_error{"Flux wasn't calculated"};
    return _flux;
}

template<std::floating_point T, std::integral I>
const parameter_2d<T>& heat_equation_solution_2d<T, I>::parameter(const std::string& group) const {
    return _parameters.at(group);
}

template<std::floating_point T, std::integral I>
T heat_equation_solution_2d<T, I>::calc_energy() const {
    //return mesh::integrate(*_base::mesh_proxy(), temperature());
    return 0;
}

template<std::floating_point T, std::integral I>
bool heat_equation_solution_2d<T, I>::is_flux_calculated() const noexcept {
    return !_flux[X].empty() && !_flux[Y].empty();
}

template<std::floating_point T, std::integral I>
T heat_equation_solution_2d<T, I>::evaluate(const coefficient_t<T>& conductivity, const size_t qshift) const {
    return std::visit(visitor{
        [](const T value) noexcept { return value; },
        [this, qshift](const spatial_dependency<T>& value) { return value(_base::mesh().quad_coord(qshift)); },
        [this, qshift](const solution_dependency<T>& value) { return value(_solution[qshift], _base::mesh().quad_coord(qshift)); }
    }, conductivity);
}

template<std::floating_point T, std::integral I>
std::array<T, 2> heat_equation_solution_2d<T, I>::evaluate(const orthotropic_conductivity_t<T>& conductivity, const size_t qshift) const {
    return { evaluate(conductivity[X], qshift), evaluate(conductivity[Y], qshift) };
}

template<std::floating_point T, std::integral I>
metamath::types::square_matrix<T, 2> heat_equation_solution_2d<T, I>::evaluate(const anisotropic_conductivity_t<T>& conductivity, const size_t qshift) const {
    return {
        evaluate(conductivity[X][X], qshift), evaluate(conductivity[X][Y], qshift),
        evaluate(conductivity[Y][X], qshift), evaluate(conductivity[Y][Y], qshift)
    };
}

template<std::floating_point T, std::integral I>
std::array<std::vector<T>, 2> heat_equation_solution_2d<T, I>::local_flux_in_qnodes() const {
    auto flux = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
    for (const auto& [group, parameter] : _parameters)
        for(const size_t e : _base::mesh().container().elements(group))
            for(const size_t qshift : _base::mesh().quad_shifts_count(e)) {
                std::visit(visitor{
                    [&](const isotropic_conductivity_t<T>& conductivity) { 
                        const T cond = evaluate(conductivity, qshift);
                        flux[X][qshift] *= -cond;
                        flux[Y][qshift] *= -cond;
                    },
                    [&](const orthotropic_conductivity_t<T>& conductivity) {
                        const auto cond = evaluate(conductivity, qshift);
                        flux[X][qshift] *= -cond[X];
                        flux[Y][qshift] *= -cond[Y];
                    },
                    [&](const anisotropic_conductivity_t<T>& conductivity) {
                        const auto cond = evaluate(conductivity, qshift);
                        std::tie(flux[X][qshift], flux[Y][qshift]) = std::make_tuple(
                            -cond[X][X] * flux[X][qshift] - cond[X][Y] * flux[Y][qshift],
                            -cond[Y][X] * flux[X][qshift] - cond[Y][Y] * flux[Y][qshift]
                        );
                    }
                }, parameter.conductivity);
            }
    return flux;
}

template<std::floating_point T, std::integral I>
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