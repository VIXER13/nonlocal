#ifndef NONLOCAL_HEAT_EQUATION_SOLUTION_2D_HPP
#define NONLOCAL_HEAT_EQUATION_SOLUTION_2D_HPP

#include "solution_2d.hpp"
#include "thermal_parameters_2d.hpp"

#include "mesh_2d_utils.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class heat_equation_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;

    std::vector<T> _temperature;
    std::array<std::vector<T>, 2> _flux;
    std::unordered_map<std::string, parameter_2d<T>> _parameters;

    void add_flux(const material_t material, const metamath::types::square_matrix<T, 2>& factor, 
                  const std::array<T, 2>& gradient, const size_t qshift);
    template<class Influence>
    std::array<T, 2> calc_nonlocal_gradient(const size_t eL,
                                            const std::array<std::vector<T>, 2>& gradient,
                                            const Influence& influence);

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
    void calc_flux();

    void save_as_vtk(std::ofstream& output) const override;
    void save_as_vtk(const std::filesystem::path& path) const;
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
void heat_equation_solution_2d<T, I>::add_flux(const material_t material, const metamath::types::square_matrix<T, 2>& factor,
                                               const std::array<T, 2>& gradient, const size_t qshift) {
    switch(material) {
        case material_t::ISOTROPIC:
            _flux[X][qshift] += factor[X][X] * gradient[X];
            _flux[Y][qshift] += factor[X][X] * gradient[Y];
        break;

        case material_t::ORTHOTROPIC:
            _flux[X][qshift] += factor[X][X] * gradient[X];
            _flux[Y][qshift] += factor[Y][Y] * gradient[Y];
        break;

        case material_t::ANISOTROPIC:
            _flux[X][qshift] += factor[X][X] * gradient[X] + factor[X][Y] * gradient[Y];
            _flux[Y][qshift] += factor[Y][X] * gradient[X] + factor[Y][Y] * gradient[Y];
        break;
    
        default:
            throw std::domain_error{"Unknown material type."};
    }
}

template<class T, class I>
template<class Influence>
std::array<T, 2> heat_equation_solution_2d<T, I>::calc_nonlocal_gradient(const size_t eL,
                                                                         const std::array<std::vector<T>, 2>& gradient,
                                                                         const Influence& influence) {
    std::array<T, 2> nonlocal_gradient = {};
    for(const size_t eNL : _base::mesh().neighbours(eL)) {
        const auto& elNL = _base::mesh().container().element_2d(eNL);
        const size_t qshiftNL = _base::mesh().quad_shift(eNL);
        for(const size_t qNL : elNL.qnodes()) {
            const T influence_weight = elNL.weight(qNL) * mesh::jacobian(_base::mesh().jacobi_matrix(qshiftNL + qNL)) *
                                       influence(_base::mesh().quad_coord(qshiftNL + qNL));
            nonlocal_gradient[X] += influence_weight * gradient[X][qshiftNL + qNL];
            nonlocal_gradient[Y] += influence_weight * gradient[Y][qshiftNL + qNL];
        }
    }
    return nonlocal_gradient;
}

template<class T, class I>
void heat_equation_solution_2d<T, I>::calc_flux() {
    const auto gradient = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
    _flux[X].resize(gradient[X].size(), T{0});
    _flux[Y].resize(gradient[Y].size(), T{0});
    for(const auto& [group, parameter] : _parameters) {
        using namespace metamath::functions;
        const model_parameters<2, T>& model = _base::model(group);
        const metamath::types::square_matrix<T, 2> local_factor = -model.local_weight * parameter.conductivity;
        const metamath::types::square_matrix<T, 2> nonlocal_factor = -nonlocal_weight(model.local_weight) * parameter.conductivity;
        const auto elements = _base::mesh().container().elements(group);
#pragma omp parallel for default(none) shared(gradient, model, parameter, local_factor, nonlocal_factor, elements) schedule(dynamic)
        for(size_t eL = elements.front(); eL < *elements.end(); ++eL)
            for(const size_t qshiftL : std::ranges::iota_view{_base::mesh().quad_shift(eL), _base::mesh().quad_shift(eL + 1)}) {
                if (theory_type(model.local_weight) == theory_t::NONLOCAL) {
                    const auto influence = [&influence = model.influence, &coordL = _base::mesh().quad_coord(qshiftL)]
                                           (const std::array<T, 2>& coordNL) { return influence(coordL, coordNL); };
                    add_flux(parameter.material, nonlocal_factor, calc_nonlocal_gradient(eL, gradient, influence), qshiftL);
                }
                add_flux(parameter.material, local_factor, {gradient[X][qshiftL], gradient[Y][qshiftL]}, qshiftL);
            }
    }
    _flux[X] = mesh::utils::qnodes_to_nodes(_base::mesh(), _flux[X]);
    _flux[Y] = mesh::utils::qnodes_to_nodes(_base::mesh(), _flux[Y]);
}

template<class T, class I>
void heat_equation_solution_2d<T, I>::save_as_vtk(std::ofstream& output) const {
    _base::save_as_vtk(output);
    _base::save_scalars(output, temperature(), "Temperature");
    if (is_flux_calculated())
        _base::save_vectors(output, flux(), "Flux");
}

template<class T, class I>
void heat_equation_solution_2d<T, I>::save_as_vtk(const std::filesystem::path& path) const {
    std::ofstream output{path};
    save_as_vtk(output);
}

}

#endif