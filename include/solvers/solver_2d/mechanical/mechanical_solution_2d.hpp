#ifndef NONLOCAL_STRUCTURAL_SOLUTION_HPP
#define NONLOCAL_STRUCTURAL_SOLUTION_HPP

#include "mechanical_parameters_2d.hpp"
#include "solution_2d.hpp"
#include <ranges>

namespace nonlocal::mechanical {

template<class T, class I>
class mechanical_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;
    using Influence_Function = typename solution_2d<T, I>::Influence_Function;

    enum : size_t { _11, _22, _12, COMPONENTS_COUNT };
    enum class collect : bool { ONLY_STRAIN, WITH_STRESS };

    std::array<std::vector<T>, 2> _displacement;
    std::array<std::vector<T>, 3> _strain, _stress;
    equation_parameters<T> _parameters;

    template<collect Type>
    void collect_solution();
    void strain_and_stress_loc();
    void stress_nonloc();

public:
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);
    template<class Vector>
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                    const T local_weight, const Influence_Function& influence_function,
                                    const equation_parameters<T>& parameters, const Vector& displacement);
    ~mechanical_solution_2d() noexcept override = default;

    const std::array<std::vector<T>, 2>& displacement() const noexcept;
    const std::array<std::vector<T>, 3>& strain() const noexcept;
    const std::array<std::vector<T>, 3>& stress() const noexcept;

    T calc_energy() const;
    bool is_strain_and_stress_calculated() const noexcept;
    void calc_strain_and_stress();

    void save_as_vtk(std::ofstream& output) const override;
    void save_as_vtk(const std::string& path) const;
};

template<class T, class I>
mechanical_solution_2d<T, I>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy)
    : _base{mesh_proxy} {
    for(std::vector<T>& displacement : _displacement)
        displacement.resize(_base::mesh_proxy()->mesh().nodes_count(), 0);
}

template<class T, class I>
template<class Vector>
mechanical_solution_2d<T, I>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                                     const T local_weight, const Influence_Function& influence_function,
                                                     const equation_parameters<T>& parameters, const Vector& displacement)
    : _base{mesh_proxy, local_weight, influence_function}
    , _parameters{parameters} {
    for(std::vector<T>& displacement : _displacement)
        displacement.resize(_base::mesh_proxy()->mesh().nodes_count(), 0);
    for(const size_t i : std::views::iota(size_t{0}, _base::mesh_proxy()->mesh().nodes_count())) {
        _displacement[X][i] = displacement[2 * i];
        _displacement[Y][i] = displacement[2 * i + 1];
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
    if(is_strain_and_stress_calculated()) {
        for(size_t e = 0; e < _base::mesh_proxy()->mesh().elements_count(); ++e) {
            const auto& el     = _base::mesh_proxy()->mesh().element_2d(e);
            const auto J_start = _base::mesh_proxy()->jacobi_matrix(e);
            for(size_t i = 0; i < el->nodes_count(); ++i) {
                auto J = J_start;
                for(size_t q = 0; q < el->qnodes_count(); ++q, ++J) {
                    const size_t node = _base::mesh_proxy()->mesh().node_number(e, i);
                    integral += el->weight(q) * el->qN(i, q) * _base::mesh_proxy()->jacobian(*J) *
                                (    strain()[_11][node] * stress()[_11][node] +
                                     strain()[_22][node] * stress()[_22][node] +
                                 2 * strain()[_12][node] * stress()[_12][node]);
                }
            }
        }
    }
    return 0.5 * integral;
}

template<class T, class I>
bool mechanical_solution_2d<T, I>::is_strain_and_stress_calculated() const noexcept {
    return !strain()[_11].empty() && !strain()[_22].empty() && !strain()[_12].empty() &&
           !stress()[_11].empty() && !stress()[_22].empty() && !stress()[_12].empty();
}

template<class T, class I>
template<typename mechanical_solution_2d<T, I>::collect Type>
void mechanical_solution_2d<T, I>::collect_solution() {
    for(std::vector<T>& strain : _strain)
        strain = MPI_utils::all_to_all<T>(strain, _base::mesh_proxy()->ranges());
    if constexpr (Type == collect::WITH_STRESS)
        for(std::vector<T>& stress : _stress)
            stress = MPI_utils::all_to_all<T>(stress, _base::mesh_proxy()->ranges());
}

template<class T, class I>
void mechanical_solution_2d<T, I>::strain_and_stress_loc() {
    {
        auto [du1dx1, du1dx2] = _base::mesh_proxy()->gradient(displacement()[X]);
        auto [du2dx1, du2dx2] = _base::mesh_proxy()->gradient(displacement()[Y]);
        _strain[_11] = std::move(du1dx1);
        _strain[_22] = std::move(du2dx2);
        _strain[_12] = std::move(du1dx2);
        for(const size_t node : std::views::iota(_base::mesh_proxy()->first_node(), _base::mesh_proxy()->last_node()))
            _strain[_12][node] = 0.5 * (_strain[_12][node] + du2dx1[node]);
    }

    const std::array<T, 3> D = hooke_matrix(_parameters);
    for(std::vector<T>& stress : _stress)
        stress.resize(_base::mesh_proxy()->mesh().nodes_count());
    for(const size_t node : std::views::iota(_base::mesh_proxy()->first_node(), _base::mesh_proxy()->last_node())) {
        _stress[_11][node] = D[_11] * _strain[_11][node] + D[_22] * _strain[_22][node];
        _stress[_22][node] = D[_22] * _strain[_11][node] + D[_11] * _strain[_22][node];
        _stress[_12][node] = D[_12] * _strain[_12][node];
    }
}

template<class T, class I>
void mechanical_solution_2d<T, I>::stress_nonloc() {
    const std::array<std::vector<T>, 3> strains_in_quads{
        _base::mesh_proxy()->approx_in_quad(_strain[0]),
        _base::mesh_proxy()->approx_in_quad(_strain[1]),
        _base::mesh_proxy()->approx_in_quad(_strain[2])
    };
    const std::array<T, 3> D = hooke_matrix(_parameters);
    _base::template calc_nonlocal(
        [this, &strains_in_quads, &D, nonlocal_weight = 1 - _base::local_weight()](const size_t e, const size_t node) {
            std::array<T, 3> integral = {};
            auto J = _base::mesh_proxy()->jacobi_matrix(e);
            auto qshift = _base::mesh_proxy()->quad_shift(e);
            auto qcoord = _base::mesh_proxy()->quad_coord(e);
            const auto& eNL = _base::mesh_proxy()->mesh().element_2d(e);
            for(size_t q = 0; q < eNL->qnodes_count(); ++q, ++qshift, ++qcoord, ++J) {
                const T influence_weight = eNL->weight(q) * _base::mesh_proxy()->jacobian(*J) *
                                           _base::influence_function()(*qcoord, _base::mesh_proxy()->mesh().node(node));
                integral[_11] += influence_weight * (D[_11] * strains_in_quads[_11][qshift] + D[_22] * strains_in_quads[_22][qshift]);
                integral[_22] += influence_weight * (D[_22] * strains_in_quads[_11][qshift] + D[_11] * strains_in_quads[_22][qshift]);
                integral[_12] += influence_weight *  D[_12] * strains_in_quads[_12][qshift];
            }
            _stress[_11][node] += nonlocal_weight * integral[_11];
            _stress[_22][node] += nonlocal_weight * integral[_22];
            _stress[_12][node] += nonlocal_weight * integral[_12];
        }
    );
}

template<class T, class I>
void mechanical_solution_2d<T, I>::calc_strain_and_stress() {
    strain_and_stress_loc();
    if(_base::local_weight() < MAX_NONLOCAL_WEIGHT<T>) {
        using namespace metamath::function;
        for(std::vector<T>& stress : _stress)
            stress *= _base::local_weight();
        collect_solution<collect::ONLY_STRAIN>();
        stress_nonloc();
    }
    collect_solution<collect::WITH_STRESS>();
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

        std::vector<T> mises(_base::mesh_proxy()->mesh().nodes_count(), 0);
        for(size_t i = 0; i < mises.size(); ++i)
            mises[i] = std::sqrt(stress()[_11][i] * stress()[_11][i] -
                                 stress()[_11][i] * stress()[_22][i] +
                                 stress()[_22][i] * stress()[_22][i] +
                             3 * stress()[_12][i] * stress()[_12][i]);
        _base::save_scalars(output, mises, "mises");
    }
}

template<class T, class I>
void mechanical_solution_2d<T, I>::save_as_vtk(const std::string& path) const {
    std::ofstream output{path};
    save_as_vtk(output);
}

}

#endif