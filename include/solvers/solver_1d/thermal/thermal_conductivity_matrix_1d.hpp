#ifndef THERMAL_CONDUCTIVITY_MATRIX_1D_HPP
#define THERMAL_CONDUCTIVITY_MATRIX_1D_HPP

#include "finite_element_matrix_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class thermal_conductivity_matrix_1d : public finite_element_matrix_1d<T, I> {
    using _base = finite_element_matrix_1d<T, I>;

protected:
    T integrate_basic(const size_t e, const size_t i) const;
    template<template<class, auto...> class Physical, auto... Args>
    T integrate_loc(const equation_parameters<1, T, Physical, Args...>& parameters, 
                    const size_t e, const size_t i, const size_t j) const;
    template<template<class, auto...> class Physical, auto... Args>
    T integrate_nonloc(const equation_parameters<1, T, Physical, Args...>& parameters, 
                       const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL) const;

    void neumann_problem_col_fill();
    void create_matrix_portrait(const std::array<bool, 2> is_first_kind, const bool is_neumann);

public:
    explicit thermal_conductivity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~thermal_conductivity_matrix_1d() override = default;

    template<template<class, auto...> class Physical, auto... Args>
    void calc_matrix(const std::vector<equation_parameters<1, T, Physical, Args...>>& parameters,
                     const std::array<bool, 2> is_first_kind, const bool is_neumann = false);
};

template<class T, class I>
thermal_conductivity_matrix_1d<T, I>::thermal_conductivity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _base{mesh} {}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : std::ranges::iota_view{size_t{0}, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q);
    return integral * _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
template<template<class, auto...> class Physical, auto... Args>
T thermal_conductivity_matrix_1d<T, I>::integrate_loc(const equation_parameters<1, T, Physical, Args...>& parameters, 
                                                      const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t q : std::ranges::iota_view{size_t{0}, el.qnodes_count()})
        integral += el.weight(q) * el.qNxi(i, q) * el.qNxi(j, q);
    return parameters.model.local_weight * parameters.physical.conductivity * integral /
           _base::mesh().jacobian(_base::mesh().segment_number(e));
}

template<class T, class I>
template<template<class, auto...> class Physical, auto... Args>
T thermal_conductivity_matrix_1d<T, I>::integrate_nonloc(const equation_parameters<1, T, Physical, Args...>& parameters, 
                                                         const size_t eL, const size_t eNL,
                                                         const size_t iL, const size_t jNL) const {
    T integral = T{0};
    const auto& el = _base::mesh().element();
    for(const size_t qL : std::ranges::iota_view{size_t{0}, el.qnodes_count()}) {
        T inner_integral = T{0};
        const T qcoordL = _base::mesh().qnode_coord(eL, qL);
        for(const size_t qNL : std::ranges::iota_view{size_t{0}, el.qnodes_count()}) {
            const T qcoordNL = _base::mesh().qnode_coord(eNL, qNL);
            inner_integral += el.weight(qNL) * parameters.model.influence(qcoordL, qcoordNL) * el.qNxi(jNL, qNL);
        }
        integral += el.weight(qL) * el.qNxi(iL, qL) * inner_integral;
    }
    return nonlocal_weight(parameters.model.local_weight) * parameters.physical.conductivity * integral;
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::neumann_problem_col_fill() {
    #pragma omp parallel for default(none)
    for(size_t node = 0; node < _base::mesh().nodes_count(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node, _base::mesh().nodes_count());
        for(const auto& [e, i] : _base::mesh().node_elements(node).to_array())
            if (e)
                val += integrate_basic(e, i);
    }
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::create_matrix_portrait(const std::array<bool, 2> is_first_kind, const bool is_neumann) {
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{size_t{0}, size_t(_base::matrix_inner().rows())})
            _base::matrix_inner().outerIndexPtr()[row + 1] = 1;
    _base::create_matrix_portrait(is_first_kind);
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{size_t{0}, size_t(_base::matrix_inner().rows())})
            _base::matrix_inner().innerIndexPtr()[_base::matrix_inner().outerIndexPtr()[row + 1] - 1] = _base::mesh().nodes_count();
}

template<class T, class I>
template<template<class, auto...> class Physical, auto... Args>
void thermal_conductivity_matrix_1d<T, I>::calc_matrix(const std::vector<equation_parameters<1, T, Physical, Args...>>& parameters,
                                                       const std::array<bool, 2> is_first_kind, const bool is_neumann) {
    if (parameters.size() != _base::mesh().segments_count())
        throw std::runtime_error{"The number of segments and the number of material parameters do not match."};
    _base::clear();
    const size_t matrix_size = _base::mesh().nodes_count() + is_neumann;
    _base::matrix_inner().resize(matrix_size, matrix_size);
    create_matrix_portrait(is_first_kind, is_neumann);
    _base::template calc_matrix(is_first_kind, theories_types(parameters),
        [this, &parameters](const size_t segment, const size_t e, const size_t i, const size_t j) {
            return integrate_loc(parameters[segment], e, i, j);
        },
        [this, &parameters](const size_t segment, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            return integrate_nonloc(parameters[segment], eL, eNL, iL, jNL);
        }
    );
    if (is_neumann)
        neumann_problem_col_fill();
}

}

#endif