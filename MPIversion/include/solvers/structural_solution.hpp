#ifndef STRUCTURAL_SOLUTION_HPP
#define STRUCTURAL_SOLUTION_HPP

#include "mesh_info.hpp"

namespace nonlocal::structural {

template<class T>
struct parameters final {
    T nu = 0, // Коэффициент Пуассона
      E  = 0; // Модуль Юнга
};

// Матрица Гука, которая имеет следующий портрет:
// arr[0] arr[1]   0
// arr[1] arr[0]   0
//   0      0    arr[2]
template<class T>
std::array<T, 3> hooke_matrix(const parameters<T>& params) noexcept {
    return {             params.E / (1 - params.nu*params.nu),
             params.nu * params.E / (1 - params.nu*params.nu),
                   0.5 * params.E / (1 + params.nu) };
}

template<class T, class I>
struct solution final {
    std::shared_ptr<mesh::mesh_info<T, I>> _mesh_info;
    std::array<T, 3> _D;
    T _p1 = 1;
    std::function<T(const std::array<T, 2>& x, const std::array<T, 2>& y)> _influence_fun;

    std::array<std::vector<T>, 2> _u;
    std::array<std::vector<T>, 3> _strain, _stress;

    void strain_and_stress_loc();
    void stress_nonloc();

public:
    template<class Vector>
    explicit solution(const std::shared_ptr<mesh::mesh_info<T, I>>& mesh_info, const std::array<T, 3>& D, const Vector& u,
                      const T p1, const std::function<T(const std::array<T, 2>& x, const std::array<T, 2>& y)>& influence_fun)
            : _mesh_info{mesh_info}
            , _D{D}
            , _p1{p1}
            , _influence_fun{influence_fun}
            , _u{std::vector<T>(mesh_info->mesh().nodes_count()), std::vector<T>(mesh_info->mesh().nodes_count())} {
        for(size_t comp = 0; comp < _u.size(); ++comp)
            for(size_t i = 0; i < _mesh_info->mesh().nodes_count(); ++i)
                _u[comp][i] = u[2*i+comp];
    }

    void calc_strain_and_stress();
    T calc_energy() const;
    void save_as_vtk(const std::string& path) const;
};

template<class T, class I>
void solution<T, I>::strain_and_stress_loc() {
    for(size_t el = 0; el < _mesh_info->mesh().elements_count(); ++el) {
        const auto& e = _mesh_info->mesh().element_2d(el);
        for(size_t i = 0; i < e->nodes_count(); ++i) {
            std::array<T, 4> jacobi_matrix = {};
            for(size_t j = 0; j < e->nodes_count(); ++j) {
                const std::array<T, 2>& node = _mesh_info->mesh().node(_mesh_info->mesh().node_number(el, j));
                jacobi_matrix[0] += node[0] * e->Nxi (j, e->node(i));
                jacobi_matrix[1] += node[0] * e->Neta(j, e->node(i));
                jacobi_matrix[2] += node[1] * e->Nxi (j, e->node(i));
                jacobi_matrix[3] += node[1] * e->Neta(j, e->node(i));
            }

            std::array<T, 3> strain_loc = {};
            for(size_t j = 0; j < e->nodes_count(); ++j) {
                const T jac = _mesh_info->jacobian(jacobi_matrix),
                        dx1 = ( jacobi_matrix[3] * e->Nxi(j, e->node(i)) - jacobi_matrix[2] * e->Neta(j, e->node(i))) / jac,
                        dx2 = (-jacobi_matrix[1] * e->Nxi(j, e->node(i)) + jacobi_matrix[0] * e->Neta(j, e->node(i))) / jac;
                strain_loc[0] += dx1 * _u[0][_mesh_info->mesh().node_number(el, j)];
                strain_loc[1] += dx2 * _u[1][_mesh_info->mesh().node_number(el, j)];
                strain_loc[2] += dx1 * _u[1][_mesh_info->mesh().node_number(el, j)] +
                                 dx2 * _u[0][_mesh_info->mesh().node_number(el, j)];
            }

            _stress[0][_mesh_info->mesh().node_number(el, i)] += _D[0] * strain_loc[0] + _D[1] * strain_loc[1];
            _stress[1][_mesh_info->mesh().node_number(el, i)] += _D[1] * strain_loc[0] + _D[0] * strain_loc[1];
            _stress[2][_mesh_info->mesh().node_number(el, i)] += _D[2] * strain_loc[2];
            for(size_t j = 0; j < 3; ++j)
                _strain[j][_mesh_info->mesh().node_number(el, i)] += strain_loc[j];
        }
    }

    for(size_t i = 0; i < _mesh_info->mesh().nodes_count(); ++i) {
        const size_t repeating_count = _mesh_info->nodes_elements_map(i).size();
        _strain[0][i] /=     repeating_count;
        _strain[1][i] /=     repeating_count;
        _strain[2][i] /= 2 * repeating_count;
        _stress[0][i] /=     repeating_count;
        _stress[1][i] /=     repeating_count;
        _stress[2][i] /= 2 * repeating_count;
    }
}

template<class T, class I>
void solution<T, I>::stress_nonloc() {
    const T p2 = 1 - _p1;
    const std::array<std::vector<T>, 3> strains_in_quads{
        _mesh_info->approx_in_quad(_strain[0]),
        _mesh_info->approx_in_quad(_strain[1]),
        _mesh_info->approx_in_quad(_strain[2])
    };

#pragma omp parallel for default(none)
    for(size_t node = 0; node < _mesh_info->mesh().nodes_count(); ++node) {
        std::unordered_set<I> neighbors;
        for(const I elL : _mesh_info->nodes_elements_map(node))
            for(const I elNL : _mesh_info->neighbors(elL))
                neighbors.insert(elNL);
        for(const I elNL : neighbors) {
            const auto& eNL = _mesh_info->mesh().element_2d(_mesh_info->mesh().element_2d_type(elNL));
            for(size_t q = 0, shift = _mesh_info->quad_shift(elNL); q < eNL->qnodes_count(); ++q, ++shift) {
                const T influence_weight = p2 * eNL->weight(q) * _mesh_info->jacobian(shift) * _influence_fun(_mesh_info->quad_coord(shift), _mesh_info->mesh().node(node));
                _stress[0][node] += influence_weight * (_D[0] * strains_in_quads[0][shift] + _D[1] * strains_in_quads[1][shift]);
                _stress[1][node] += influence_weight * (_D[1] * strains_in_quads[0][shift] + _D[0] * strains_in_quads[1][shift]);
                _stress[2][node] += influence_weight *  _D[2] * strains_in_quads[2][shift];
            }
        }
    }
}

template<class T, class I>
void solution<T, I>::calc_strain_and_stress() {
    for(size_t comp = 0; comp < _strain.size(); ++comp) {
        _strain[comp].resize(_mesh_info->mesh().nodes_count());
        _stress[comp].resize(_mesh_info->mesh().nodes_count());
    }

    strain_and_stress_loc();
    if(_p1 < 0.999) { // Нелокальная задача
        for(size_t comp = 0; comp < _stress.size(); ++comp)
            for(size_t node = 0; node < _mesh_info->mesh().nodes_count(); ++node)
                _stress[comp][node] *= _p1;
        stress_nonloc();
    }
}

template<class T, class I>
T solution<T, I>::calc_energy() const {
    T integral = 0;
    if(!_strain[0].empty()) {
        for(size_t el = 0; el < _mesh_info->mesh().elements_count(); ++el) {
            const auto& e = _mesh_info->mesh().element_2d(el);
            for(size_t i = 0; i < e->nodes_count(); ++i)
                for(size_t q = 0, shift = _mesh_info->quad_shift(el); q < e->qnodes_count(); ++q, ++shift) {
                    const size_t node = _mesh_info->mesh().node_number(el, i);
                    integral += e->weight(q) * e->qN(i, q) * _mesh_info->jacobian(shift) *
                                (    _strain[0][node] * _stress[0][node] +
                                     _strain[1][node] * _stress[1][node] +
                                 2 * _strain[2][node] * _stress[2][node]);
                }
        }
    }
    return 0.5 * integral;
}

template<class T, class I>
void solution<T, I>::save_as_vtk(const std::string& path) const {
    static constexpr std::string_view data_type = std::is_same_v<T, float> ? "float" : "double";
    static constexpr std::array<std::string_view, 3> strain_number = {"strain11", "strain22", "strain12"},
                                                     stress_number = {"stress11", "stress22", "stress12"};

    std::ofstream fout{path};
    fout.precision(std::numeric_limits<T>::max_digits10);

    _mesh_info->mesh().save_as_vtk(fout);

    fout << "POINT_DATA " << _mesh_info->mesh().nodes_count() << '\n';
    fout << "VECTORS Displacement " << data_type << '\n';
    for(size_t i = 0; i < _mesh_info->mesh().nodes_count(); ++i)
        fout << _u[0][i] << ' ' << _u[1][i] << " 0\n";

    if(!_strain[0].empty()) {
        for(size_t comp = 0; comp < 3; ++comp) {
            fout << "SCALARS " << strain_number[comp] << ' ' << data_type << " 1\n"
                 << "LOOKUP_TABLE default\n";
            for(size_t i = 0; i < _mesh_info->mesh().nodes_count(); ++i)
                fout << _strain[comp][i] << '\n';
        }

        for(size_t comp = 0; comp < 3; ++comp) {
            fout << "SCALARS " << stress_number[comp] << ' ' << data_type << " 1\n"
                 << "LOOKUP_TABLE default\n";
            for(size_t i = 0; i < _mesh_info->mesh().nodes_count(); ++i)
                fout << _stress[comp][i] << '\n';
        }

        fout << "SCALARS mises " << data_type << " 1\n"
             << "LOOKUP_TABLE default\n";
        for(size_t i = 0; i < _mesh_info->mesh().nodes_count(); ++i)
            fout << std::sqrt(_stress[0][i] * _stress[0][i] +
                              _stress[1][i] * _stress[1][i] -
                              _stress[0][i] * _stress[1][i] +
                          3 * _stress[2][i] * _stress[2][i]) << '\n';
    }
}

}

#endif