#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/base/solution_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T, std::integral I = uint32_t>
class heat_equation_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;

    std::vector<T> _temperature;
    std::array<std::vector<T>, 2> _flux;
    std::unordered_map<std::string, evaluated_conductivity_t<T>> _conductivity;

    std::array<std::vector<T>, 2> local_flux_in_qnodes() const;

public:
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    template<class Vector>
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                       const evaluated_conductivity_2d<T>& conductivity,
                                       const Vector& temperature);
    ~heat_equation_solution_2d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::array<std::vector<T>, 2>& flux() const;
    const evaluated_conductivity_t<T>& conductivity(const std::string& group) const;

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
                                                           const evaluated_conductivity_2d<T>& conductivity,
                                                           const Vector& temperature)
    : _base{mesh, get_models(conductivity)}
    , _temperature{temperature.cbegin(), std::next(temperature.cbegin(), mesh->container().nodes_count())}
    , _conductivity{get_physical_parameters(conductivity)} {}

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
const evaluated_conductivity_t<T>& heat_equation_solution_2d<T, I>::conductivity(const std::string& group) const {
    return _conductivity.at(group);
}

template<std::floating_point T, std::integral I>
bool heat_equation_solution_2d<T, I>::is_flux_calculated() const noexcept {
    return !_flux[X].empty() && !_flux[Y].empty();
}

template<std::floating_point T, std::integral I>
std::array<std::vector<T>, 2> heat_equation_solution_2d<T, I>::local_flux_in_qnodes() const {
    auto flux = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
    for (const auto& [group, conductivity] : _conductivity)
        for(const size_t e : _base::mesh().container().elements(group))
            for(const size_t qshift : _base::mesh().quad_shifts_count(e)) {
                std::visit(metamath::types::visitor{
                    [&flux, qshift](const evaluated_isotropic_conductivity_t<T>& conductivity) { 
                        const T conduct = conductivity.index() ? std::get<Nonconstant>(conductivity)[qshift] :
                                                                 std::get<Constant>(conductivity);
                        flux[X][qshift] *= -conduct;
                        flux[Y][qshift] *= -conduct;
                    },
                    [&flux, qshift](const evaluated_orthotropic_conductivity_t<T>& conductivity) {
                        const auto& conduct = conductivity.index() ? std::get<Nonconstant>(conductivity)[qshift] :
                                                                     std::get<Constant>(conductivity);
                        flux[X][qshift] *= -conduct[X];
                        flux[Y][qshift] *= -conduct[Y];
                    },
                    [&flux, qshift](const evaluated_anisotropic_conductivity_t<T>& conductivity) {
                        const auto& conduct = conductivity.index() ? std::get<Nonconstant>(conductivity)[qshift] :
                                                                     std::get<Constant>(conductivity);
                        std::tie(flux[X][qshift], flux[Y][qshift]) = std::make_tuple(
                            -conduct[XX] * flux[X][qshift] - conduct[XY] * flux[Y][qshift],
                            -conduct[XY] * flux[X][qshift] - conduct[YY] * flux[Y][qshift]
                        );
                    }
                }, conductivity);
            }
    return flux;
}

template<std::floating_point T, std::integral I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::calc_flux() {
    if (is_flux_calculated())
        return _flux;

    _flux = local_flux_in_qnodes();
    std::array<std::vector<T>, 2> flux = _flux;
    for(const auto& [group, parameter] : _conductivity)
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
                            const T influence_weight = elNL.weight(qNL) * _base::mesh().jacobian(qshiftNL) *
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