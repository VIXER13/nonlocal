#ifndef NONLOCAL_STIFFNESS_MATRIX_2D_HPP
#define NONLOCAL_STIFFNESS_MATRIX_2D_HPP

#include "../finite_element_matrix_2d.hpp"

namespace nonlocal::mechanical {

template<class T, class I, class Matrix_Index>
class stiffness_matrix : public finite_element_matrix_2d<1, T, I, Matrix_Index> {
    using _base = finite_element_matrix_2d<1, T, I, Matrix_Index>;
    using _base::component::X;
    using _base::component::Y;

protected:
    static void add_to_pair(std::array<T, 4>& pairs, const std::array<T, 2>& wdNd, const std::array<T, 2>& dNd) noexcept;
    static std::array<T, 4> calc_block(const std::array<T, 3>& D, const std::array<T, 4>& pairs) noexcept;

    std::array<T, 4> integrate_loc(const std::array<T, 3>& D, const size_t e, const size_t i, const size_t j) const;

    template<class Influence_Function>
    std::array<T, 4> integrate_nonloc(const std::array<T, 3>& D,
                                      const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                      const Influence_Function& influence_function) const;

public:
    explicit stiffness_matrix(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);
    ~stiffness_matrix() noexcept override = default;

    template<class Influence_Function>
    void calc_matrix(const std::array<T, 3>& D,
                     const std::vector<bool>& is_inner,
                     const T p1, const Influence_Function& influence_fun);
};

template<class T, class I, class Matrix_Index>
stiffness_matrix<T, I, Matrix_Index>::stiffness_matrix(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy)
    : _base{mesh_proxy} {}

template<class T, class I, class Matrix_Index>
void stiffness_matrix<T, I, Matrix_Index>::add_to_pair(std::array<T, 4>& pairs, const std::array<T, 2>& wdNd, const std::array<T, 2>& dNd) noexcept {
    pairs[0] += wdNd[X] * dNd[X];
    pairs[1] += wdNd[X] * dNd[Y];
    pairs[2] += wdNd[Y] * dNd[X];
    pairs[3] += wdNd[Y] * dNd[Y];
}

template<class T, class I, class Matrix_Index>
std::array<T, 4> stiffness_matrix<T, I, Matrix_Index>::calc_block(const std::array<T, 3>& D, const std::array<T, 4>& pairs) noexcept {
    return {
        D[0] * pairs[0] + D[2] * pairs[3],
        D[1] * pairs[1] + D[2] * pairs[2],
        D[1] * pairs[2] + D[2] * pairs[1],
        D[0] * pairs[3] + D[2] * pairs[0],
    };
}

template<class T, class I, class Matrix_Index>
std::array<T, 4> stiffness_matrix<T, I, Matrix_Index>::integrate_loc(const std::array<T, 3>& D, const size_t e, const size_t i, const size_t j) const {
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    std::array<T, 4> pairs = {};
    for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj) {
        const T weight = el->weight(q) / _base::jacobian(*J);
        const std::array<T, 2> wdNdi = {weight * (*dNdi)[X], weight * (*dNdi)[Y]};
        add_to_pair(pairs, wdNdi, *dNdj);
    }
    return calc_block(D, pairs);
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
std::array<T, 4> stiffness_matrix<T, I, Matrix_Index>::integrate_nonloc(const std::array<T, 3>& D,
                                                                        const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                                                        const Influence_Function& influence_function) const {
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
    const auto  qcoordNL_begin = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_begin    = _base::mesh_proxy()->dNdX(eNL, jNL);
    std::array<T, 4> pairs = {};
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto dNdNL    = dNdNL_begin;
        auto qcoordNL = qcoordNL_begin;
        std::array<T, 2> inner_int = {};
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
            const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
            inner_int[X] += influence_weight * (*dNdNL)[X];
            inner_int[Y] += influence_weight * (*dNdNL)[Y];
        }
        const std::array<T, 2> wdNdL = {elL->weight(qL) * (*dNdL)[X], elL->weight(qL) * (*dNdL)[Y]};
        add_to_pair(pairs, wdNdL, inner_int);
    }
    return calc_block(D, pairs);
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
void stiffness_matrix<T, I, Matrix_Index>::calc_matrix(const std::array<T, 3>& D,
                                                       const std::vector<bool>& is_inner,
                                                       const T p1, const Influence_Function& influence_fun) {
    const theory_t theory = p1 < MAX_NONLOCAL_WEIGHT<T> ? theory_t::NONLOCAL : theory_t::LOCAL;
    const size_t rows = 2 * (_base::mesh_proxy()->last_node() - _base::mesh_proxy()->first_node()),
                 cols = 2 * _base::mesh().nodes_count();
    _base::matrix_inner().resize(rows, cols);
    _base::matrix_bound().resize(rows, cols);
    if (theory == theory_t::LOCAL)
        _base::template create_matrix_portrait<theory_t::LOCAL>(is_inner);
    else if (theory == theory_t::NONLOCAL)
        _base::template create_matrix_portrait<theory_t::NONLOCAL>(is_inner);
    _base::template calc_matrix<2>(is_inner, theory, influence_fun,
        [this, p1, &D](const size_t e, const size_t i, const size_t j) {
            using namespace metamath::function;
            return p1 * integrate_loc(D, e, i, j);
        },
        [this, p2 = 1 - p1, &D](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            using namespace metamath::function;
            return p2 * integrate_nonloc(D, eL, eNL, iL, jNL, influence_function);
        }
    );
}

}

#endif