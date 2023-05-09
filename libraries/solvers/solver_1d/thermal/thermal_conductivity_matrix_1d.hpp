#ifndef NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_1D_HPP
#define NONLOCAL_THERMAL_CONDUCTIVITY_MATRIX_1D_HPP

#include "../../equation_parameters.hpp"

#include "finite_element_matrix_1d.hpp"
#include "thermal_parameters_1d.hpp"

#include <optional>

namespace nonlocal::thermal {

template<class T, class I>
class thermal_conductivity_matrix_1d : public finite_element_matrix_1d<T, I> {
    using _base = finite_element_matrix_1d<T, I>;

protected:
    T integrate_basic(const size_t e, const size_t i) const;

    template<class Integrator>
    T integrate_loc(const Integrator& integrator) const;
    T integrate_loc(const T conductivity, const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const std::function<T(const T)>& conductivity, const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const std::function<T(const T, const T)>& conductivity, const std::vector<T>& solution,
                    const size_t e, const size_t i, const size_t j) const;

    template<class Integrator>
    T integrate_nonloc(const size_t eL, const size_t iL, const Integrator& integrator) const;
    template<class Influence_Function>
    T integrate_nonloc(const T conductivity, const Influence_Function& influence,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;
    template<class Influence_Function>
    T integrate_nonloc(const std::function<T(const T)>& conductivity, const Influence_Function& influence,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;
    template<class Influence_Function>
    T integrate_nonloc(const std::function<T(const T, const T)>& conductivity, const Influence_Function& influence,
                       const std::vector<T>& solution,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;


    void neumann_problem_col_fill();
    void create_matrix_portrait(const std::vector<theory_t>& theories,
                                const std::array<bool, 2> is_first_kind, const bool is_neumann);

public:
    explicit thermal_conductivity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~thermal_conductivity_matrix_1d() override = default;

    void calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2> is_first_kind, const bool is_neumann = false,
                     const std::optional<std::vector<T>>& solution = std::nullopt);
};

template<class T, class I>
thermal_conductivity_matrix_1d<T, I>::thermal_conductivity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _base{mesh} {}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
template<class Integrator>
T thermal_conductivity_matrix_1d<T, I>::integrate_loc(const Integrator& integrator) const {
    const auto qnodes = _base::mesh().element().qnodes();
    return std::accumulate(qnodes.begin(), qnodes.end(), T{0}, integrator);
}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_loc(const T conductivity, const size_t e, const size_t i, const size_t j) const {
    const T integral = integrate_loc([&el = _base::mesh().element(), i, j](const T integral, const size_t q) {
        return integral + el.weight(q) * el.qNxi(i, q) * el.qNxi(j, q);
    });
    return conductivity * integral / _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_loc(
    const std::function<T(const T)>& conductivity, const size_t e, const size_t i, const size_t j) const {
    const T integral = integrate_loc([this, &conductivity, e, i, j](const T integral, const size_t q) {
        const auto& el = _base::mesh().element();
        const T qcoord = _base::mesh().qnode_coord(e, q);
        return integral + el.weight(q) * el.qNxi(i, q) * el.qNxi(j, q) * conductivity(qcoord);
    });
    return integral / _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_loc(
    const std::function<T(const T, const T)>& conductivity, const std::vector<T>& solution, 
    const size_t e, const size_t i, const size_t j) const {
    const T integral = integrate_loc([this, &conductivity, &solution, e, i, j](const T integral, const size_t q) {
        const auto& el = _base::mesh().element();
        const T qcoord = _base::mesh().qnode_coord(e, q);
        const size_t qshift = _base::mesh().qnode_number(e, q); 
        return integral + el.weight(q) * el.qNxi(i, q) * el.qNxi(j, q) * conductivity(qcoord, solution[qshift]);
    });
    return integral / _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
template<class Integrator>
T thermal_conductivity_matrix_1d<T, I>::integrate_nonloc(const size_t eL, const size_t iL, const Integrator& integrator) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    const auto qnodes = el.qnodes();
    for(const size_t qL : qnodes) {
        const T inner_integral = std::accumulate(qnodes.begin(), qnodes.end(), T{0},
            [&integrator, qcoordL = _base::mesh().qnode_coord(eL, qL)](const T integral, const size_t qNL) {
                return integral + integrator(qcoordL, qNL);
            });
        integral += el.weight(qL) * el.qNxi(iL, qL) * inner_integral;
    }
    return integral;
}

template<class T, class I>
template<class Influence_Function>
T thermal_conductivity_matrix_1d<T, I>::integrate_nonloc(
    const T conductivity, const Influence_Function& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    return conductivity * integrate_nonloc(eL, iL, [this, &influence, eNL, jNL](const T qcoordL, const size_t qNL) {
        const auto& el = _base::mesh().element();
        const T qcoordNL = _base::mesh().qnode_coord(eNL, qNL);
        return el.weight(qNL) * influence(qcoordL, qcoordNL) * el.qNxi(jNL, qNL);
    });
}

template<class T, class I>
template<class Influence_Function>
T thermal_conductivity_matrix_1d<T, I>::integrate_nonloc(
    const std::function<T(const T)>& conductivity, const Influence_Function& influence,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    return integrate_nonloc(eL, iL, [this, &conductivity, &influence, eNL, jNL](const T qcoordL, const size_t qNL) {
        const auto& el = _base::mesh().element();
        const T qcoordNL = _base::mesh().qnode_coord(eNL, qNL);
        return el.weight(qNL) * influence(qcoordL, qcoordNL) * conductivity(qcoordNL) * el.qNxi(jNL, qNL);
    });
}

template<class T, class I>
template<class Influence_Function>
T thermal_conductivity_matrix_1d<T, I>::integrate_nonloc(
    const std::function<T(const T, const T)>& conductivity, const Influence_Function& influence, const std::vector<T>& solution,
    const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    return integrate_nonloc(eL, iL, [this, &conductivity, &influence, &solution, eNL, jNL](const T qcoordL, const size_t qNL) {
        const auto& el = _base::mesh().element();
        const T qcoordNL = _base::mesh().qnode_coord(eNL, qNL);
        const size_t qshiftNL = _base::mesh().qnode_number(eNL, qNL);
        return el.weight(qNL) * influence(qcoordL, qcoordNL) * conductivity(qcoordNL, solution[qshiftNL]) * el.qNxi(jNL, qNL);
    });
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::neumann_problem_col_fill() {
#pragma omp parallel for default(none)
    for(size_t node = 0; node < _base::mesh().nodes_count(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node, _base::mesh().nodes_count());
        for(const auto& [e, i] : _base::mesh().node_elements(node).to_array())
            if (e) val += integrate_basic(e, i);
    }
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::create_matrix_portrait(const std::vector<theory_t>& theories,
                                                                  const std::array<bool, 2> is_first_kind, const bool is_neumann) {
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, size_t(_base::matrix_inner().rows())})
            _base::matrix_inner().outerIndexPtr()[row + 1] = 1;
    _base::create_matrix_portrait(theories, is_first_kind);
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, size_t(_base::matrix_inner().rows())})
            _base::matrix_inner().innerIndexPtr()[_base::matrix_inner().outerIndexPtr()[row + 1] - 1] = _base::mesh().nodes_count();
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::calc_matrix(const parameters_1d<T>& parameters, const std::array<bool, 2> is_first_kind, const bool is_neumann,
                                                       const std::optional<std::vector<T>>& solution) {
    if (parameters.size() != _base::mesh().segments_count())
        throw std::runtime_error{"The number of segments and the number of material parameters do not match."};
    _base::clear();
    const size_t matrix_size = _base::mesh().nodes_count() + is_neumann;
    _base::matrix_inner().resize(matrix_size, matrix_size);
    const std::vector<theory_t> theories = theories_types(parameters);
    create_matrix_portrait(theories, is_first_kind, is_neumann);
    _base::template calc_matrix(theories, is_first_kind,
        [this, &parameters, &solution](const size_t segment, const size_t e, const size_t i, const size_t j) {
            using enum coefficients_t;
            const auto& [model, physic] = parameters[segment];
            if (const auto* const parameter = parameter_cast<CONSTANTS>(physic.get()); parameter)
                return model.local_weight * integrate_loc(parameter->conductivity, e, i, j);
            if (const auto* const parameter = parameter_cast<SPACE_DEPENDENT>(physic.get()); parameter)
                return model.local_weight * integrate_loc(parameter->conductivity, e, i, j);
            if (const auto* const parameter = parameter_cast<SOLUTION_DEPENDENT>(physic.get()); parameter)
                return model.local_weight * integrate_loc(parameter->conductivity, *solution, e, i, j);
            return std::numeric_limits<T>::quiet_NaN();
        },
        [this, &parameters, &solution](const size_t segment, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            using enum coefficients_t;
            const auto& [model, physic] = parameters[segment];
            const T nonlocal_weight = nonlocal::nonlocal_weight(model.local_weight);
            if (const auto* const parameter = parameter_cast<CONSTANTS>(physic.get()); parameter)
                return nonlocal_weight * integrate_nonloc(parameter->conductivity, model.influence, eL, eNL, iL, jNL);
            if (const auto* const parameter = parameter_cast<SPACE_DEPENDENT>(physic.get()); parameter)
                return nonlocal_weight * integrate_nonloc(parameter->conductivity, model.influence, eL, eNL, iL, jNL);
            if (const auto* const parameter = parameter_cast<SOLUTION_DEPENDENT>(physic.get()); parameter)
                return nonlocal_weight * integrate_nonloc(parameter->conductivity, model.influence, *solution, eL, eNL, iL, jNL);
            return std::numeric_limits<T>::quiet_NaN();
        }
    );
    if (is_neumann)
        neumann_problem_col_fill();
}

}

#endif