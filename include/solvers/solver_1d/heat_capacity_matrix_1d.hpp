#ifndef HEAT_CAPACITY_MATRIX_1D_HPP
#define HEAT_CAPACITY_MATRIX_1D_HPP

#include "finite_element_matrix_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
class heat_capacity_matrix_1d : public finite_element_matrix_1d<T, I> {
    using _base = finite_element_matrix_base_1d<T, I>;

protected:
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

public:
    explicit heat_capacity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    ~heat_capacity_matrix_1d() override = default;

    void calc_matrix(const T c, const T rho, const std::array<boundary_condition_t, 2> bound_cond);
};

template<class T, class I>
heat_capacity_matrix_1d<T, I>::heat_capacity_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh) :
    : finite_element_matrix_1d<T, I>{mesh} {}

template<class T, class I>
T heat_capacity_matrix_1d<T, I>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->nodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q);
    return integral * mesh()->jacobian();
}

template<class T, class I>
void heat_capacity_matrix_1d<T, I>::calc_matrix(const T c, const T rho, const std::array<boundary_condition_t, 2> bound_cond) {
    _base::clear_matrix();
    _base::matrix_inner().resize(_base::mesh()->nodes_count(), _base::mesh()->nodes_count());
    create_matrix_portrait(bound_cond, theory_t::LOCAL);
    static constexpr auto no_influence = []() constexpr noexcept {};
    _base::template calc_matrix(bound_cond, theory_t::LOCAL, no_influence,
        [this, factor = c * rho](const size_t e, const size_t i, const size_t j) { return factor * integrate_basic_pair(e, i, j); },
        [](const size_t, const size_t, const size_t, const size_t, const decltype(no_influence)) constexpr noexcept { return 0; }
    );
}

}

#endif