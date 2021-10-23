#ifndef THERMAL_CONDUCTIVITY_MATRIX_1D_HPP
#define THERMAL_CONDUCTIVITY_MATRIX_1D_HPP

#include "finite_element_matrix_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class thermal_conductivity_matrix_1d : public finite_element_matrix_1d<T, I> {
    using _base = finite_element_matrix_1d<T, I>;

protected:
    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;
    template<class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(const std::array<boundary_condition_t, 2> bound_cond,
                                const theory_t theory, const bool is_neumann);

    void neumann_problem_col_fill();

public:
    explicit thermal_conductivity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~thermal_conductivity_matrix_1d() override = default;

    template<class Influence_Function>
    void calc_matrix(const T lambda, const T p1, const Influence_Function& influence_function,
                     const std::array<boundary_condition_t, 2> bound_cond, const bool is_neumann = false);
};

template<class T, class I>
thermal_conductivity_matrix_1d<T, I>::thermal_conductivity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : finite_element_matrix_1d<T, I>{mesh} {}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = _base::mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q);
    return integral * _base::mesh()->jacobian();
}

template<class T, class I>
T thermal_conductivity_matrix_1d<T, I>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = _base::mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qNxi(i, q) * el->qNxi(j, q);
    return integral / _base::mesh()->jacobian();
}

template<class T, class I>
template<class Influence_Function>
T thermal_conductivity_matrix_1d<T, I>::integrate_nonloc(const size_t eL, const size_t eNL,
                                                         const size_t iL, const size_t jNL,
                                                         const Influence_Function& influence_function) const {
    T integral = T{0};
    const auto& el = _base::mesh()->element();
    for(size_t qL = 0; qL < el->qnodes_count(); ++qL) {
        T inner_integral = T{0};
        const T qcoordL = _base::mesh()->quad_coord(eL, qL);
        for(size_t qNL = 0; qNL < el->qnodes_count(); ++qNL) {
            const T qcoordNL = _base::mesh()->quad_coord(eNL, qNL);
            inner_integral += el->weight(qNL) * influence_function(qcoordL, qcoordNL) * el->qNxi(jNL, qNL);
        }
        integral += el->weight(qL) * el->qNxi(iL, qL) * inner_integral;
    }
    return integral;
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::create_matrix_portrait(const std::array<boundary_condition_t, 2> bound_cond,
                                                                  const theory_t theory, const bool is_neumann) {
    if (is_neumann)
        for(size_t row = 0; row < _base::matrix_inner().rows(); ++row)
            _base::matrix_inner().outerIndexPtr()[row+1] = 1;
    _base::create_matrix_portrait(bound_cond, theory);
    if (is_neumann)
        for(size_t row = 0; row < _base::matrix_inner().rows(); ++row)
            _base::matrix_inner().innerIndexPtr()[_base::matrix_inner().outerIndexPtr()[row+1]-1] = _base::mesh()->nodes_count();
}

template<class T, class I>
void thermal_conductivity_matrix_1d<T, I>::neumann_problem_col_fill() {
#pragma omp parallel for default(none)
    for(size_t node = 0; node < _base::mesh()->nodes_count(); ++node) {
        T& val = _base::matrix_inner().coeffRef(node, _base::mesh()->nodes_count());
        for(const auto& [e, i] : _base::mesh()->node_elements(node).arr)
            if(e != std::numeric_limits<size_t>::max())
                val += integrate_basic(e, i);
    }
}

template<class T, class I>
template<class Influence_Function>
void thermal_conductivity_matrix_1d<T, I>::calc_matrix(const T lambda, const T p1, const Influence_Function& influence_function,
                                                       const std::array<boundary_condition_t, 2> bound_cond, const bool is_neumann) {
    _base::clear_matrix();
    const theory_t theory = p1 < MAX_NONLOCAL_WEIGHT<T> ? theory_t::NONLOCAL : theory_t::LOCAL;
    const size_t matrix_size = _base::mesh()->nodes_count() + is_neumann;
    _base::matrix_inner().resize(matrix_size, matrix_size);
    create_matrix_portrait(bound_cond, theory, is_neumann);
    _base::template calc_matrix(bound_cond, theory, influence_function,
        [this, factor = lambda * p1](const size_t e, const size_t i, const size_t j) {
            return factor * integrate_loc(e, i, j);
        },
        [this, factor = lambda * (T{1} - p1)]
        (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            return factor * integrate_nonloc(eL, eNL, iL, jNL, influence_function);
        }
    );
    if (is_neumann)
        neumann_problem_col_fill();
}

}

#endif