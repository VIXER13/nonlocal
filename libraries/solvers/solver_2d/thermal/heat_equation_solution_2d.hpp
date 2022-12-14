#ifndef NONLOCAL_HEAT_EQUATION_SOLUTION_2D_HPP
#define NONLOCAL_HEAT_EQUATION_SOLUTION_2D_HPP

#include "solution_2d.hpp"
#include "thermal_parameters_2d.hpp"

#include "mesh_2d_utils.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class heat_equation_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;
    using typename solution_2d<T, I>::influence_function_t;

    std::vector<T> _temperature;
    std::array<std::vector<T>, 2> _flux;
    parameter_2d<T> _parameters;

    void flux_from_gradient(std::array<std::vector<T>, 2>& gradient, const parameter_2d<T>& parameters, const T model_weight) const;
    void calc_local_flux(const std::array<std::vector<T>, 2>& gradient);
    void calc_nonlocal_flux(const std::array<std::vector<T>, 2>& gradient);

public:
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    template<class Vector>
    explicit heat_equation_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                       const T local_weight, const influence_function_t& influence_function,
                                       const parameter_2d<T>& parameters, const Vector& temperature);
    ~heat_equation_solution_2d() noexcept override = default;

    const std::vector<T>& temperature() const noexcept;
    const std::array<std::vector<T>, 2>& flux() const;
    const parameter_2d<T>& parameters() const;

    T calc_energy() const;
    bool is_flux_calculated() const noexcept;
    const std::array<std::vector<T>, 2>& calc_flux();

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
                                                           const T local_weight, const influence_function_t& influence_function,
                                                           const parameter_2d<T>& parameters, const Vector& temperature)
    : _base{mesh, local_weight, influence_function}
    , _temperature{temperature.cbegin(), std::next(temperature.cbegin(), mesh->container().nodes_count())}
    , _parameters{parameters} {}

template<class T, class I>
const std::vector<T>& heat_equation_solution_2d<T, I>::temperature() const noexcept {
    return _temperature;
}

template<class T, class I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::flux() const {
    return _flux;
}

template<class T, class I>
const parameter_2d<T>& heat_equation_solution_2d<T, I>::parameters() const {
    return _parameters;
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
void heat_equation_solution_2d<T, I>::flux_from_gradient(std::array<std::vector<T>, 2>& gradient, const parameter_2d<T>& parameters, const T model_weight) const {
    using namespace metamath::functions;
    const metamath::types::square_matrix<T, 2> factor = { 
        -model_weight * parameters.conductivity[X], 
        -model_weight * parameters.conductivity[Y]
    };
    switch (parameters.material) {
        case material_t::ISOTROPIC:
            gradient[X] *= factor[X][X];
            gradient[Y] *= factor[X][X];
        break;

        case material_t::ORTHOTROPIC:
            gradient[X] *= factor[X][X];
            gradient[Y] *= factor[Y][Y];
        break;

        case material_t::ANISOTROPIC:
            for(const size_t i : _base::mesh().container().nodes())
                for(const size_t comp : std::ranges::iota_view{0u, 2u})
                    gradient[comp] = factor[comp][X] * gradient[X] + factor[comp][Y] * gradient[Y];
        break;
    
        default:
            throw std::domain_error{"Unknown material type."};
    }
}

template<class T, class I>
void heat_equation_solution_2d<T, I>::calc_local_flux(const std::array<std::vector<T>, 2>& gradient) {
    _flux[X] = mesh::utils::qnodes_to_nodes(_base::mesh(), gradient[X]);
    _flux[Y] = mesh::utils::qnodes_to_nodes(_base::mesh(), gradient[Y]);
    flux_from_gradient(_flux, parameters(), _base::local_weight());
}

template<class T, class I>
void heat_equation_solution_2d<T, I>::calc_nonlocal_flux(const std::array<std::vector<T>, 2>& gradient) {
    std::array<std::vector<T>, 2> nonlocal_gradient = {
        std::vector<T>(gradient[X].size(), T{0}),
        std::vector<T>(gradient[Y].size(), T{0})
    };

    for(const size_t eL : _base::mesh().container().elements_2d()) {
        const auto& elL = _base::mesh().container().element_2d(eL);
        const size_t qshiftL = _base::mesh().quad_shift(eL);
        for(const size_t eNL : _base::mesh().neighbours(eL)) {
            const auto& elNL = _base::mesh().container().element_2d(eNL);
            const size_t qshiftNL = _base::mesh().quad_shift(eNL);
            for(const size_t qL : std::ranges::iota_view{0u, elL.qnodes_count()})
                for(const size_t qNL : std::ranges::iota_view{0u, elNL.qnodes_count()}) {
                    const T influence_weight = elNL.weight(qNL) * mesh::jacobian(_base::mesh().jacobi_matrix(qshiftNL + qNL)) *
                                               _base::influence_function()(_base::mesh().quad_coord(qshiftL  + qL ), 
                                                                           _base::mesh().quad_coord(qshiftNL + qNL));
                    nonlocal_gradient[X][qshiftL + qL] += influence_weight * gradient[X][qshiftNL + qNL];
                    nonlocal_gradient[Y][qshiftL + qL] += influence_weight * gradient[Y][qshiftNL + qNL];
                }
        }
    }

    nonlocal_gradient[X] = mesh::utils::qnodes_to_nodes(_base::mesh(), nonlocal_gradient[X]);
    nonlocal_gradient[Y] = mesh::utils::qnodes_to_nodes(_base::mesh(), nonlocal_gradient[Y]);
    flux_from_gradient(nonlocal_gradient, parameters(), nonlocal_weight(_base::local_weight()));
    using namespace metamath::functions;
    _flux[X] += nonlocal_gradient[X];
    _flux[Y] += nonlocal_gradient[Y];
}

template<class T, class I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::calc_flux() {
    if (!is_flux_calculated()) {
        const auto gradient = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
        calc_local_flux(gradient);
        if (theory_type(_base::local_weight()) == theory_t::NONLOCAL)
            calc_nonlocal_flux(gradient);
    }
    return _flux;
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