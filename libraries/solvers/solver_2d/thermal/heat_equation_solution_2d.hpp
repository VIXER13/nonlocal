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
void heat_equation_solution_2d<T, I>::calc_local_flux(const std::array<std::vector<T>, 2>& gradient) {
    _flux[X] = mesh::utils::qnodes_to_nodes(_base::mesh(), gradient[X]);
    _flux[Y] = mesh::utils::qnodes_to_nodes(_base::mesh(), gradient[Y]);
    using namespace metamath::functions;
    const metamath::types::square_matrix<T, 2> factor = { 
        _base::local_weight() * parameters().conductivity[X], 
        _base::local_weight() * parameters().conductivity[Y]
    };
    switch (parameters().material) {
        case material_t::ISOTROPIC:
            _flux[X] *= factor[X][X];
            _flux[Y] *= factor[X][X];
        break;

        case material_t::ORTHOTROPIC:
            _flux[X] *= factor[X][X];
            _flux[Y] *= factor[Y][Y];
        break;

        case material_t::ANISOTROPIC:
            for(const size_t i : _base::mesh().container().nodes())
                for(const size_t comp : std::ranges::iota_view{0u, 2u})
                    _flux[comp] = factor[comp][X] * _flux[X] + factor[comp][Y] * _flux[Y];
        break;
    
        default:
            throw std::domain_error{"Unknown material type."};
    }
}

template<class T, class I>
void heat_equation_solution_2d<T, I>::calc_nonlocal_flux(const std::array<std::vector<T>, 2>& gradient_in_quads) {
    // using namespace metamath::functions;
    // const std::array<T, 2> nonlocal_factor = (T{1} - _base::local_weight()) * _thermal_conductivity;
    // _base::template calc_nonlocal(
    //     [this, &gradient_in_quads, &nonlocal_factor](const size_t e, const size_t node) {
    //         std::array<T, 2> integral = {};
    //         auto J = _base::mesh()->jacobi_matrix(e);
    //         auto qshift = _base::mesh()->quad_shift(e);
    //         auto qcoord = _base::mesh()->quad_coord(e);
    //         const auto& eNL = _base::mesh()->mesh().element_2d(e);
    //         for(size_t q = 0; q < eNL->qnodes_count(); ++q, ++qshift, ++qcoord, ++J) {
    //             const T influence_weight = eNL->weight(q) * mesh::jacobian(*J) *
    //                                        _base::influence_function()(*qcoord, _base::mesh_proxy()->mesh().node(node));
    //             integral[X] += influence_weight * gradient_in_quads[X][qshift];
    //             integral[Y] += influence_weight * gradient_in_quads[Y][qshift];
    //         }
    //         _flux[X][node] += nonlocal_factor[X] * integral[X];
    //         _flux[Y][node] += nonlocal_factor[Y] * integral[Y];
    //     }
    // );
}

template<class T, class I>
const std::array<std::vector<T>, 2>& heat_equation_solution_2d<T, I>::calc_flux() {
    // if (!is_flux_calculated()) {
    //     const std::array<std::vector<T>, 2> gradient_in_quads = mesh::approximate_gradient_in_qnodes(*_base::mesh_proxy(), temperature());
    //     calc_local_flux(gradient_in_quads);
    //     if (_base::local_weight() < MAX_NONLOCAL_WEIGHT<T>)
    //         calc_nonlocal_flux(gradient_in_quads);
    // }
    if (!is_flux_calculated()) {
        const auto gradient = mesh::utils::gradient_in_qnodes(_base::mesh(), _temperature);
        calc_local_flux(gradient);
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