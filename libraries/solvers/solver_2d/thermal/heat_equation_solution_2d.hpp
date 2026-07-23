#pragma once

#include "thermal_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/base/solution_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
class heat_equation_solution_2d : public solution_2d<T> {
    using _base = solution_2d<T>;

    std::vector<T> _temperature;
    std::vector<std::array<T, 2>> _flux;
    std::unordered_map<std::string, evaluated_conductivity_t<T>> _conductivity;

    std::vector<std::array<T, 2>> local_flux_in_qnodes() const;

public:
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                       const evaluated_conductivity_2d<T>& conductivity,
                                       const std::vector<T>& temperature);
    ~heat_equation_solution_2d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::vector<std::array<T, 2>>& flux() const;
    const evaluated_conductivity_t<T>& conductivity(const std::string& group) const;

    bool is_flux_calculated() const noexcept;
    const std::vector<std::array<T, 2>>& calc_flux();
};

template<std::floating_point T>
heat_equation_solution_2d<T>::heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _base{mesh}
    , _temperature(mesh->container().nodes_count(), T{0}) {}

template<std::floating_point T>
heat_equation_solution_2d<T>::heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                                        const evaluated_conductivity_2d<T>& conductivity,
                                                        const std::vector<T>& temperature)
    : _base{mesh, get_models(conductivity)}
    , _temperature{temperature.cbegin(), std::next(temperature.cbegin(), mesh->container().nodes_count())}
    , _conductivity{get_physical_parameters(conductivity)} {}

template<std::floating_point T>
const std::vector<T>& heat_equation_solution_2d<T>::temperature() const noexcept {
    return _temperature;
}

template<std::floating_point T>
const std::vector<std::array<T, 2>>& heat_equation_solution_2d<T>::flux() const {
    if (!is_flux_calculated())
        throw std::runtime_error{"Flux wasn't calculated"};
    return _flux;
}

template<std::floating_point T>
const evaluated_conductivity_t<T>& heat_equation_solution_2d<T>::conductivity(const std::string& group) const {
    return _conductivity.at(group);
}

template<std::floating_point T>
bool heat_equation_solution_2d<T>::is_flux_calculated() const noexcept {
    return !_flux.empty();
}

template<std::floating_point T>
std::vector<std::array<T, 2>> heat_equation_solution_2d<T>::local_flux_in_qnodes() const {
    auto flux = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
    for (const auto& [group, conductivity] : _conductivity)
        for(const size_t e : _base::mesh().container().elements(group))
            for(const size_t qshift : _base::mesh().quad_shifts_count(e)) {
                std::visit([&flux, qshift](const auto& conductivity) {
                    const auto& conduct = conductivity.index() ? std::get<Variable>(conductivity)[qshift] :
                                                                 std::get<Constant>(conductivity);
                    using namespace metamath::operators;
                    using conductivity_t = std::remove_cvref_t<decltype(conductivity)>;
                    if constexpr (std::is_same_v<conductivity_t, evaluated_isotropic_conductivity_t<T>>)
                        flux[qshift] *= -conduct;
                    else if constexpr (std::is_same_v<conductivity_t, evaluated_orthotropic_conductivity_t<T>>)
                        flux[qshift] = {-conduct[X] * flux[qshift][X], -conduct[Y] * flux[qshift][Y]};
                    else if constexpr (std::is_same_v<conductivity_t, evaluated_anisotropic_conductivity_t<T>>)
                        flux[qshift] = {-conduct[XX] * flux[X][qshift] - conduct[XY] * flux[Y][qshift],
                                        -conduct[XY] * flux[X][qshift] - conduct[YY] * flux[Y][qshift]};
                    else
                        static_assert(false, "Unknown conductivity coefficients type.");
                }, conductivity);
            }
    return flux;
}

template<std::floating_point T>
const std::vector<std::array<T, 2>>& heat_equation_solution_2d<T>::calc_flux() {
    if (is_flux_calculated())
        return _flux;

    using namespace metamath::operators;
    _flux = local_flux_in_qnodes();
    auto flux = _flux;
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
                            nonlocal_gradient += influence_weight * _flux[qshiftNL];
                            ++qshiftNL;
                        }
                    }
                    nonlocal_gradient *= nonlocal_weight;
                    flux[qshiftL] *= model.local_weight;
                    flux[qshiftL] += nonlocal_gradient;
                }
        }
    _flux = mesh::utils::qnodes_to_nodes(_base::mesh(), flux);
    _flux = parallel::all_to_all(_flux, _base::mesh().MPI_ranges());
    return _flux;
}

}