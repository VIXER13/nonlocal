#ifndef NONLOCAL_HEAT_EQUATION_SOLUTION_HPP
#define NONLOCAL_HEAT_EQUATION_SOLUTION_HPP

#include "mesh_2d.hpp"
#include "heat_equation_parameters_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class solution final {
    std::shared_ptr<mesh::mesh_proxy<T, I>> _mesh_proxy;
    std::vector<T> _temperature;
    std::array<std::vector<T>, 2> _flux;
    std::array<T, 2> _lambda = {1, 1};
    T _p1 = T{1};
    std::function<T(const std::array<T, 2>& x, const std::array<T, 2>& y)>
        _influence_fun = [](const std::array<T, 2>&, const std::array<T, 2>&) constexpr noexcept { return T{}; };

    void calc_local_flux();
    void calc_nonlocal_flux();

public:
    explicit solution(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);

    template<material_t Material, class Influence_Function, class Vector>
    explicit solution(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                      const equation_parameters_2d<T, Material>& parameters,
                      const T p1, const Influence_Function& influence_fun,
                      const Vector& temperature);

    const std::vector<T>& temperature() const noexcept;
    const std::array<std::vector<T>, 2>& flux() const noexcept;
    const std::shared_ptr<mesh::mesh_proxy<T, I>> mesh_proxy() const noexcept;

    T calc_energy() const;

    void calc_flux();

    void save_as_vtk(const std::string& path) const;
};

template<class T, class I>
solution<T, I>::solution(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy)
    : _mesh_proxy{mesh_proxy} {
    _temperature.resize(_mesh_proxy->mesh().nodes_count(), 0);
}

template<class T, class I>
template<material_t Material, class Influence_Function, class Vector>
solution<T, I>::solution(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                         const equation_parameters_2d<T, Material>& parameters,
                         const T p1, const Influence_Function& influence_fun,
                         const Vector& temperature)
    : _mesh_proxy{mesh_proxy}
    , _p1{p1}
    , _influence_fun{influence_fun} {
    if constexpr (Material == material_t::ISOTROPIC)
        _lambda = {parameters.lambda, parameters.lambda};
    else
        _lambda = parameters.lambda;
    _temperature.resize(_mesh_proxy->mesh().nodes_count());
    for(size_t i = 0; i < _mesh_proxy->mesh().nodes_count(); ++i)
        _temperature[i] = temperature[i];
}

template<class T, class I>
const std::vector<T>& solution<T, I>::temperature() const noexcept {
    return _temperature;
}

template<class T, class I>
const std::array<std::vector<T>, 2>& solution<T, I>::flux() const noexcept {
    return _flux;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_proxy<T, I>> solution<T, I>::mesh_proxy() const noexcept {
    return _mesh_proxy;
}

template<class T, class I>
T solution<T, I>::calc_energy() const {
    return _mesh_proxy->integrate_solution(_temperature);
}

template<class T, class I>
void solution<T, I>::calc_local_flux() {
    using namespace metamath::function;
    const std::array<T, 2> local_factor = _p1 * _lambda;
    _flux[X] *= local_factor[X];
    _flux[Y] *= local_factor[Y];
}

template<class T, class I>
void solution<T, I>::calc_nonlocal_flux() {
    using namespace metamath::function;
    const std::array<T, 2> nonlocal_factor = (T{1} - _p1) * _lambda;
    const std::array<std::vector<T>, 2> gradient_in_quads{
        _mesh_proxy->approx_in_quad(_flux[0]),
        _mesh_proxy->approx_in_quad(_flux[1])
    };
    calc_local_flux();
    std::vector<bool> neighbors(_mesh_proxy->mesh().elements_count(), false);
#pragma omp parallel for default(none) shared(nonlocal_factor, gradient_in_quads) firstprivate(neighbors)
    for(size_t node = _mesh_proxy->first_node(); node < _mesh_proxy->last_node(); ++node) {
        std::fill(neighbors.begin(), neighbors.end(), false);
        for(const I eL : _mesh_proxy->nodes_elements_map(node))
            for(const I eNL : _mesh_proxy->neighbors(eL))
                neighbors[eNL] = true;
        std::array<T, 2> nonlocal_integral = {};
        for(size_t e = 0; e < _mesh_proxy->mesh().elements_count(); ++e)
            if (neighbors[e]) {
                const auto& eNL = _mesh_proxy->mesh().element_2d(e);
                auto qshift     = _mesh_proxy->quad_shift(e);
                auto qcoord     = _mesh_proxy->quad_coord(e);
                auto J          = _mesh_proxy->jacobi_matrix(e);
                for(size_t q = 0; q < eNL->qnodes_count(); ++q, ++qshift, ++qcoord, ++J) {
                    const T influence_weight = eNL->weight(q) * _mesh_proxy->jacobian(*J) *
                                               _influence_fun(*qcoord, _mesh_proxy->mesh().node(node));
                    nonlocal_integral[X] += influence_weight * gradient_in_quads[X][qshift];
                    nonlocal_integral[Y] += influence_weight * gradient_in_quads[Y][qshift];
                }
            }
        _flux[X][node] += nonlocal_factor[X] * nonlocal_integral[X];
        _flux[Y][node] += nonlocal_factor[Y] * nonlocal_integral[Y];
    }
}

template<class T, class I>
void solution<T, I>::calc_flux() {
    _flux = _mesh_proxy->template gradient(_temperature);
    if (_p1 < MAX_NONLOCAL_WEIGHT<T>)
        calc_nonlocal_flux();
    else
        calc_local_flux();
}

template<class T, class I>
void solution<T, I>::save_as_vtk(const std::string& path) const {
    static constexpr std::string_view data_type = std::is_same_v<T, float> ? "float" : "double";
    std::ofstream fout{path};
    fout.precision(std::numeric_limits<T>::max_digits10);

    _mesh_proxy->mesh().save_as_vtk(fout);

    fout << "POINT_DATA " << _mesh_proxy->mesh().nodes_count() << '\n';
    fout << "SCALARS Temperature " << data_type << " 1\n"
         << "LOOKUP_TABLE default\n";
    for(const T val : _temperature)
        fout << val << '\n';

    if (!_flux[X].empty() && !_flux[Y].empty()) {
        fout << "VECTORS Flux " << data_type << '\n';
        for (size_t i = 0; i < _mesh_proxy->mesh().nodes_count(); ++i)
            fout << _flux[X][i] << ' ' << _flux[Y][i] << " 0\n";
    }
}

}

#endif