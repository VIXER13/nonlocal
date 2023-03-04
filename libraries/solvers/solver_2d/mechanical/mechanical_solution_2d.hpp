#ifndef NONLOCAL_STRUCTURAL_SOLUTION_HPP
#define NONLOCAL_STRUCTURAL_SOLUTION_HPP

#include "solution_2d.hpp"
#include "mechanical_parameters_2d.hpp"

#include "mesh_2d_utils.hpp"

namespace nonlocal::mechanical {

template<class T, class I>
class mechanical_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;

    enum : size_t { _11, _22, _12 };

    std::array<std::vector<T>, 2> _displacement;
    std::array<std::vector<T>, 3> _strain, _stress;
    std::unordered_map<std::string, parameter_2d<T>> _parameters;

    std::array<std::vector<T>, 3> strains_in_quadratures() const;
    template<class Influence>
    std::array<T, 3> calc_nonlocal_strain(const size_t eL, const std::array<std::vector<T>, 3>& strains, const Influence& influence) const;
    void add_stress(const hooke_matrix<T>& hooke, const std::array<T, 3>& strain, const size_t qshift);

public:
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    template<class Vector>
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                    const parameters_2d<T>& parameters, const Vector& displacement);
    ~mechanical_solution_2d() noexcept override = default;

    const std::array<std::vector<T>, 2>& displacement() const noexcept;
    const std::array<std::vector<T>, 3>& strain() const noexcept;
    const std::array<std::vector<T>, 3>& stress() const noexcept;

    T calc_energy() const;
    bool is_strain_and_stress_calculated() const noexcept;
    void calc_strain_and_stress();

    void save_as_vtk(std::ofstream& output) const override;
    void save_as_vtk(const std::filesystem::path& path) const;
};

template<class T, class I>
mechanical_solution_2d<T, I>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {
    for(std::vector<T>& displacement : _displacement)
        displacement.resize(mesh->container().nodes_count(), T{0});
}

template<class T, class I>
template<class Vector>
mechanical_solution_2d<T, I>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                     const parameters_2d<T>& parameters, const Vector& displacement)
    : _base{mesh, get_models(parameters)}
    , _parameters{get_physical_parameters(parameters)} {
    for(std::vector<T>& displacement : _displacement)
        displacement.resize(_base::mesh().container().nodes_count(), 0);
    for(const size_t i : std::views::iota(0u, _base::mesh().container().nodes_count())) {
        _displacement[X][i] = displacement[2 * i + X];
        _displacement[Y][i] = displacement[2 * i + Y];
    }
}

template<class T, class I>
const std::array<std::vector<T>, 2>& mechanical_solution_2d<T, I>::displacement() const noexcept {
    return _displacement;
}

template<class T, class I>
const std::array<std::vector<T>, 3>& mechanical_solution_2d<T, I>::strain() const noexcept {
    return _strain;
}

template<class T, class I>
const std::array<std::vector<T>, 3>& mechanical_solution_2d<T, I>::stress() const noexcept {
    return _stress;
}

template<class T, class I>
T mechanical_solution_2d<T, I>::calc_energy() const {
    T integral = 0;
    // if(is_strain_and_stress_calculated()) {
    //     for(size_t e = 0; e < _base::mesh_proxy()->mesh().elements_count(); ++e) {
    //         const auto& el     = _base::mesh_proxy()->mesh().element_2d(e);
    //         const auto J_start = _base::mesh_proxy()->jacobi_matrix(e);
    //         for(size_t i = 0; i < el->nodes_count(); ++i) {
    //             auto J = J_start;
    //             for(size_t q = 0; q < el->qnodes_count(); ++q, ++J) {
    //                 const size_t node = _base::mesh_proxy()->mesh().node_number(e, i);
    //                 integral += el->weight(q) * el->qN(i, q) * mesh::jacobian(*J) *
    //                             (    strain()[_11][node] * stress()[_11][node] +
    //                                  strain()[_22][node] * stress()[_22][node] +
    //                              2 * strain()[_12][node] * stress()[_12][node]);
    //             }
    //         }
    //     }
    // }
    return 0.5 * integral;
}

template<class T, class I>
bool mechanical_solution_2d<T, I>::is_strain_and_stress_calculated() const noexcept {
    return !strain()[_11].empty() && !strain()[_22].empty() && !strain()[_12].empty() &&
           !stress()[_11].empty() && !stress()[_22].empty() && !stress()[_12].empty();
}

template<class T, class I>
std::array<std::vector<T>, 3> mechanical_solution_2d<T, I>::strains_in_quadratures() const {
    auto [strain11_in_quad, strain12_in_quad] = mesh::utils::gradient_in_qnodes(_base::mesh(), displacement()[X]);
    auto [strain21_in_quad, strain22_in_quad] = mesh::utils::gradient_in_qnodes(_base::mesh(), displacement()[Y]);
    for(const size_t q : std::ranges::iota_view{size_t{0}, strain12_in_quad.size()})
        strain12_in_quad[q] = 0.5 * (strain12_in_quad[q] + strain21_in_quad[q]);
    return {std::move(strain11_in_quad), std::move(strain22_in_quad), std::move(strain12_in_quad)};
}

template<class T, class I>
template<class Influence>
std::array<T, 3> mechanical_solution_2d<T, I>::calc_nonlocal_strain(const size_t eL,
                                                                    const std::array<std::vector<T>, 3>& strains,
                                                                    const Influence& influence) const {
    std::array<T, 3> nonlocal_stress = {};
    for(const size_t eNL : _base::mesh().neighbours(eL)) {
        const auto& elNL = _base::mesh().container().element_2d(eNL);
        const size_t qshiftNL = _base::mesh().quad_shift(eNL);
        for(const size_t qNL : elNL.qnodes()) {
            const size_t qshift = qshiftNL + qNL;
            const T influence_weight = elNL.weight(qNL) * mesh::jacobian(_base::mesh().jacobi_matrix(qshift)) *
                                       influence(_base::mesh().quad_coord(qshift));
            for(const size_t i : std::ranges::iota_view{0u, nonlocal_stress.size()})
                nonlocal_stress[i] += influence_weight * strains[i][qshift];
        }
    }
    return nonlocal_stress;
}

template<class T, class I>
void mechanical_solution_2d<T, I>::add_stress(const hooke_matrix<T>& hooke, const std::array<T, 3>& strain, const size_t qshift) {
    _stress[_11][qshift] += hooke[_11] * strain[_11] + hooke[_22] * strain[_22];
    _stress[_22][qshift] += hooke[_22] * strain[_11] + hooke[_11] * strain[_22];
    _stress[_12][qshift] += hooke[_12] * strain[_12];
}

template<class T, class I>
void mechanical_solution_2d<T, I>::calc_strain_and_stress() {
    auto strains = strains_in_quadratures();
    for(const size_t i : std::ranges::iota_view{0u, 3u}) {
        _strain[i] = mesh::utils::qnodes_to_nodes(_base::mesh(), strains[i]);
        _stress[i].resize(strains[i].size(), T{0});
    }
    // strains[_11] -= temperature_strains; // TODO: thermoelasticity problem
    // strains[_22] -= temperature_strains; // TODO: thermoelasticity problem
    for(const auto& [group, parameter] : _parameters) {
        using namespace metamath::functions;
        const model_parameters<2, T>& model = _base::model(group);
        const hooke_matrix local_hooke = model.local_weight * parameter.hooke();
        const hooke_matrix nonlocal_hooke = nonlocal_weight(model.local_weight) * parameter.hooke();
        const auto elements = _base::mesh().container().elements(group);
#pragma omp parallel for default(none) shared(strains, model, local_hooke, nonlocal_hooke, elements) schedule(dynamic)
        for(size_t eL = elements.front(); eL < *elements.end(); ++eL)
            for(const size_t qshiftL : std::ranges::iota_view{_base::mesh().quad_shift(eL), _base::mesh().quad_shift(eL + 1)}) {
                if (theory_type(model.local_weight) == theory_t::NONLOCAL) {
                    const auto influence = [&influence = model.influence, &qnodeL = _base::mesh().quad_coord(qshiftL)]
                                           (const std::array<T, 2>& qnodeNL) { return influence(qnodeL, qnodeNL); };
                    add_stress(nonlocal_hooke, calc_nonlocal_strain(eL, strains, influence), qshiftL);
                }
                add_stress(local_hooke, {strains[_11][qshiftL], strains[_22][qshiftL], strains[_12][qshiftL]}, qshiftL);
            }
    }
    for(const size_t i : std::ranges::iota_view{0u, 3u})
        _stress[i] = mesh::utils::qnodes_to_nodes(_base::mesh(), _stress[i]);
}

template<class T, class I>
void mechanical_solution_2d<T, I>::save_as_vtk(std::ofstream& output) const {
    _base::save_as_vtk(output);
    _base::save_vectors(output, displacement(), "Displacement");
    if (is_strain_and_stress_calculated()) {
        _base::save_scalars(output, strain()[_11], "strain11");
        _base::save_scalars(output, strain()[_22], "strain22");
        _base::save_scalars(output, strain()[_12], "strain12");
        _base::save_scalars(output, stress()[_11], "stress11");
        _base::save_scalars(output, stress()[_22], "stress22");
        _base::save_scalars(output, stress()[_12], "stress12");

        std::vector<T> mises(_base::mesh().container().nodes_count(), T{0});
        for(const size_t i : std::ranges::iota_view{0u, mises.size()})
            mises[i] = std::sqrt(stress()[_11][i] * stress()[_11][i] -
                                 stress()[_11][i] * stress()[_22][i] +
                                 stress()[_22][i] * stress()[_22][i] +
                             3 * stress()[_12][i] * stress()[_12][i]);
        _base::save_scalars(output, mises, "mises");
    }
}

template<class T, class I>
void mechanical_solution_2d<T, I>::save_as_vtk(const std::filesystem::path& path) const {
    std::ofstream output{path};
    save_as_vtk(output);
}

}

#endif