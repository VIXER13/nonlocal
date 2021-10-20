#ifndef FINITE_ELEMENT_MATRIX_BASE_HPP
#define FINITE_ELEMENT_MATRIX_BASE_HPP

#include "../solvers_constants.hpp"
#include "../solvers_utils.hpp"
#include "mesh_1d.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace nonlocal {

template<class T, class I>
class finite_element_matrix_1d {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> _matrix_inner;
    std::array<std::unordered_map<size_t, T>, 2> _matrix_bound;

protected:
    explicit finite_element_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);

    void create_matrix_portrait(const std::array<boundary_condition_t, 2> bound_cond, const bool is_nonlocal);

    template<theory_t Theory, class Callback>
    void mesh_run(const Callback& callback) const;

    template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
    void calc_matrix(const std::array<boundary_condition_t, 2> bound_cond,
                     const bool is_nonlocal,
                     const Influence_Function& influence_fun,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc);

public:
    virtual ~finite_element_matrix_1d() = default;

    const std::shared_ptr<mesh::mesh_1d<T>>& mesh() const;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix_inner();
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix_inner() const;
    std::array<std::unordered_map<size_t, T>, 2>& matrix_bound();
    const std::array<std::unordered_map<size_t, T>, 2>& matrix_bound() const;

    void clear();
};

template<class T, class I>
finite_element_matrix_1d<T, I>::finite_element_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& finite_element_matrix_1d<T, I>::mesh() const {
    return _mesh;
}

template<class T, class I>
Eigen::SparseMatrix<T, Eigen::RowMajor, I>& finite_element_matrix_1d<T, I>::matrix_inner() {
    return _matrix_inner;
}

template<class T, class I>
const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& finite_element_matrix_1d<T, I>::matrix_inner() const {
    return _matrix_inner;
}

template<class T, class I>
std::array<std::unordered_map<size_t, T>, 2>& finite_element_matrix_1d<T, I>::matrix_bound() {
    return _matrix_bound;
}

template<class T, class I>
const std::array<std::unordered_map<size_t, T>, 2>& finite_element_matrix_1d<T, I>::matrix_bound() const {
    return _matrix_bound;
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::clear() {
    _matrix_inner = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{};
    _matrix_bound = std::array<std::unordered_map<size_t, T>, 2>{};
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::create_matrix_portrait(const std::array<boundary_condition_t, 2> bound_cond, const bool is_nonlocal) {
#pragma omp parallel for default(none) shared(bound_cond, is_nonlocal)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        if (bound_cond.front() == boundary_condition_t::FIRST_KIND && node == 0 ||
            bound_cond.back () == boundary_condition_t::FIRST_KIND && node == mesh()->nodes_count()-1)
            _matrix_inner.outerIndexPtr()[node+1] = 1;
        else {
            const auto [eL, iL, eR, iR] = mesh()->node_elements(node).named;
            const size_t e = eR == std::numeric_limits<size_t>::max() ? eL : eR,
                         i = iR == std::numeric_limits<size_t>::max() ? iL : iR,
                         right_neighbour = is_nonlocal ? mesh()->right_neighbour(e) : e + 1;
            const bool last_node_first_kind = bound_cond.back() == boundary_condition_t::FIRST_KIND &&
                                              right_neighbour * (mesh()->element()->nodes_count() - 1) == mesh()->nodes_count()-1;
            _matrix_inner.outerIndexPtr()[node+1] += (right_neighbour - e) * (mesh()->element()->nodes_count() - 1) - i + 1 - last_node_first_kind;
        }

    utils::prepare_memory(_matrix_inner);

    for(size_t i = 0; i < _matrix_inner.rows(); ++i)
        for(size_t j = _matrix_inner.outerIndexPtr()[i], k = i; j < _matrix_inner.outerIndexPtr()[i+1]; ++j, ++k)
            _matrix_inner.innerIndexPtr()[j] = k;
}

template<class T, class I>
template<theory_t Theory, class Callback>
void finite_element_matrix_1d<T, I>::mesh_run(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback) schedule(dynamic)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        for(const auto& [eL, iL] : mesh()->node_elements(node).arr)
            if(eL != std::numeric_limits<size_t>::max()) {
                if constexpr (Theory == theory_t::LOCAL)
                    for(size_t jL = 0; jL < mesh()->element()->nodes_count(); ++jL)
                        callback(eL, iL, jL);
                if constexpr (Theory == theory_t::NONLOCAL) {
                    const size_t finish = mesh()->right_neighbour(eL);
                    for(size_t eNL = mesh()->left_neighbour(eL); eNL < finish; ++eNL)
                        for(size_t jNL = 0; jNL < mesh()->element()->nodes_count(); ++jNL)
                            callback(eL, eNL, iL, jNL);
                }
            }
}

template<class T, class I>
template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
void finite_element_matrix_1d<T, I>::calc_matrix(const std::array<boundary_condition_t, 2> bound_cond,
                                                 const bool is_nonlocal,
                                                 const Influence_Function& influence_fun,
                                                 const Integrate_Loc& integrate_rule_loc,
                                                 const Integrate_Nonloc& integrate_rule_nonloc) {
    enum class coeff_destination : uint8_t {
        LEFT_BOUND,
        RIGHT_BOUND,
        REGULAR,
        NO
    };

    const auto calc_predicate = [this, bound_cond](const size_t row, const size_t col) {
        if (bound_cond.front() == boundary_condition_t::FIRST_KIND && (row == 0 || col == 0)) {
            if (row == 0)
                return coeff_destination::LEFT_BOUND;
        } else if (bound_cond.back() == boundary_condition_t::FIRST_KIND &&
                   (row == mesh()->nodes_count()-1 || col == mesh()->nodes_count()-1)) {
            if (row == mesh()->nodes_count()-1)
                return coeff_destination::RIGHT_BOUND;
        } else if (row <= col)
            return coeff_destination::REGULAR;
        return coeff_destination::NO;
    };

    const auto boundary_calc = [this](const size_t b, const size_t row, const size_t col, const T integral) {
        if (col == row)
            _matrix_inner.coeffRef(row, col) = T{1};
        else {
            const auto [it, flag] = _matrix_bound[b].template try_emplace(col, integral);
            if (!flag) it->second += integral;
        }
    };

    mesh_run<theory_t::LOCAL>(
        [this, &calc_predicate, &boundary_calc, &integrate_rule_loc](const size_t e, const size_t i, const size_t j) {
            const size_t row = mesh()->node_number(e, i),
                         col = mesh()->node_number(e, j);
            if (const coeff_destination dst = calc_predicate(row, col); dst != coeff_destination::NO) {
                const T integral = integrate_rule_loc(e, i, j);
                if (dst == coeff_destination::REGULAR) _matrix_inner.coeffRef(row, col) += integral;
                else                                   boundary_calc(size_t(dst), row, col, integral);
            }
        }
    );

    if (is_nonlocal) {
        mesh_run<theory_t::NONLOCAL>(
            [this, &calc_predicate, &boundary_calc, &integrate_rule_nonloc, &influence_fun]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const size_t row = mesh()->node_number(eL,  iL ),
                             col = mesh()->node_number(eNL, jNL);
                if (const coeff_destination dst = calc_predicate(row, col); dst != coeff_destination::NO) {
                    const T integral = integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
                    if (dst == coeff_destination::REGULAR) _matrix_inner.coeffRef(row, col) += integral;
                    else                                   boundary_calc(size_t(dst), row, col, integral);
                }
            }
        );
    }
}

}

#endif