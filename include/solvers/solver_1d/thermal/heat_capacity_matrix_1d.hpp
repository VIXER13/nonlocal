#ifndef HEAT_CAPACITY_MATRIX_1D_HPP
#define HEAT_CAPACITY_MATRIX_1D_HPP

#include "finite_element_matrix_1d.hpp"
#include "boundary_condition_1d.hpp"
#include "heat_equation_parameters_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class heat_capacity_matrix_1d : public finite_element_matrix_1d<T, I> {
    using _base = finite_element_matrix_1d<T, I>;

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

public:
    explicit heat_capacity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~heat_capacity_matrix_1d() override = default;

    void calc_matrix(const T c, const T rho, const std::array<boundary_condition_t, 2> bound_cond);
};

template<class T, class I>
heat_capacity_matrix_1d<T, I>::heat_capacity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _base{mesh} {}

template<class T, class I>
T heat_capacity_matrix_1d<T, I>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = _base::mesh()->element();
    for(const size_t q : std::views::iota(size_t{0}, el->qnodes_count()))
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q);
    return integral * _base::mesh()->jacobian();
}

template<class T, class I>
void heat_capacity_matrix_1d<T, I>::calc_matrix(const T c, const T rho, const std::array<boundary_condition_t, 2> bound_cond) {
    _base::clear_matrix();
    _base::matrix_inner().resize(_base::mesh()->nodes_count(), _base::mesh()->nodes_count());
    _base::create_matrix_portrait(utils::to_general_condition(bound_cond), theory_t::LOCAL);
    static constexpr auto NO_INFLUENCE = []() constexpr noexcept {};
    _base::template calc_matrix(utils::to_general_condition(bound_cond), theory_t::LOCAL, NO_INFLUENCE,
        [this, factor = c * rho](const size_t e, const size_t i, const size_t j) { return factor * integrate_basic_pair(e, i, j); },
        [](const size_t, const size_t, const size_t, const size_t, const decltype(NO_INFLUENCE)) constexpr noexcept { return 0; }
    );
}

}

#endif