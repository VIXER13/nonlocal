#ifndef NONLOCAL_FINITE_ELEMENT_MATRIX_1D_HPP
#define NONLOCAL_FINITE_ELEMENT_MATRIX_1D_HPP

#include "../equation_parameters.hpp"
#include "../solvers_utils.hpp"

#include "mesh_1d.hpp"
#include <eigen3/Eigen/Sparse>
#include <iostream>

namespace nonlocal {

template<class T, class I>
class finite_element_matrix_1d {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> _matrix_inner;
    std::array<std::unordered_map<size_t, T>, 2> _matrix_bound;

protected:
    explicit finite_element_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);

    void calc_shifts(const std::array<bool, 2> is_first_kind);
    void init_indices();
    void create_matrix_portrait(const std::array<bool, 2> is_first_kind);

    template<theory_t Theory, class Callback>
    void mesh_run(const size_t segment, const Callback& callback) const;
    template<template<class, auto...> class Physical, auto... Args, class Integrate_Loc, class Integrate_Nonloc>
    void calc_matrix(const std::array<bool, 2> is_first_kind,
                     const std::vector<equation_parameters<1, T, Physical, Args...>>& parameters,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc);

public:
    virtual ~finite_element_matrix_1d() = default;

    const mesh::mesh_1d<T>& mesh() const noexcept;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix_inner() noexcept;
    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix_inner() const noexcept;
    std::array<std::unordered_map<size_t, T>, 2>& matrix_bound() noexcept;
    const std::array<std::unordered_map<size_t, T>, 2>& matrix_bound() const noexcept;

    void clear();
};

template<class T, class I>
finite_element_matrix_1d<T, I>::finite_element_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
const mesh::mesh_1d<T>& finite_element_matrix_1d<T, I>::mesh() const noexcept {
    return *_mesh;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& finite_element_matrix_1d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
Eigen::SparseMatrix<T, Eigen::RowMajor, I>& finite_element_matrix_1d<T, I>::matrix_inner() noexcept {
    return _matrix_inner;
}

template<class T, class I>
const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& finite_element_matrix_1d<T, I>::matrix_inner() const noexcept {
    return _matrix_inner;
}

template<class T, class I>
std::array<std::unordered_map<size_t, T>, 2>& finite_element_matrix_1d<T, I>::matrix_bound() noexcept {
    return _matrix_bound;
}

template<class T, class I>
const std::array<std::unordered_map<size_t, T>, 2>& finite_element_matrix_1d<T, I>::matrix_bound() const noexcept {
    return _matrix_bound;
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::clear() {
    _matrix_inner = {};
    _matrix_bound = {};
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::calc_shifts(const std::array<bool, 2> is_first_kind) {
    const size_t last_node = mesh().nodes_count() - 1;
    const size_t nodes_in_element = mesh().element().nodes_count() - 1;
    for(size_t node = 0; node < mesh().nodes_count(); ++node) {
        if (is_first_kind.front() && node == 0 || is_first_kind.back() && node == last_node)
            matrix_inner().outerIndexPtr()[node + 1] = 1;
        else {
            const auto [left, right] = mesh().node_elements(node);
            const auto [e, i] = right ? right : left;
            const size_t neighbour = *mesh().neighbours(e).end();
            const bool is_last_first_kind = is_first_kind.back() && neighbour * nodes_in_element == last_node;
            matrix_inner().outerIndexPtr()[node + 1] = (neighbour - e) * nodes_in_element - i + 1 - is_last_first_kind;
        }
    }
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::init_indices() {
    for(const size_t i : std::ranges::iota_view(size_t{0}, size_t(matrix_inner().rows())))
        for(size_t j = matrix_inner().outerIndexPtr()[i], k = i; j < matrix_inner().outerIndexPtr()[i+1]; ++j, ++k)
            matrix_inner().innerIndexPtr()[j] = k;
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::create_matrix_portrait(const std::array<bool, 2> is_first_kind) {
    calc_shifts(is_first_kind);
    utils::accumulate_shifts(matrix_inner());
    std::cout << "Non-zero elements count: " << matrix_inner().outerIndexPtr()[matrix_inner().rows()] << std::endl;
    utils::allocate_matrix(matrix_inner());
    init_indices();
}

template<class T, class I>
template<theory_t Theory, class Callback>
void finite_element_matrix_1d<T, I>::mesh_run(const size_t segment, const Callback& callback) const {
    const auto correct_node_data =
        [this](const size_t node, const std::ranges::iota_view<size_t, size_t> segment_nodes) constexpr noexcept {
            auto node_elements = mesh().node_elements(node).to_array();
            if (node == segment_nodes.front())
                node_elements.front() = std::nullopt;
            else if (node == segment_nodes.back())
                node_elements.back() = std::nullopt;
            return node_elements;
        };

    const auto segment_nodes = mesh().segment_nodes(segment);
    for(const size_t node : segment_nodes) {
        for(const auto node_data : correct_node_data(node, segment_nodes))
            if (node_data) {
                const auto& [eL, iL] = node_data;
                if constexpr (Theory == theory_t::LOCAL)
                    for(const size_t jL : std::ranges::iota_view{size_t{0}, mesh().element().nodes_count()})
                        callback(eL, iL, jL);
                if constexpr (Theory == theory_t::NONLOCAL)
                    for(const size_t eNL : mesh().neighbours(eL))
                        for(const size_t jNL : std::ranges::iota_view{size_t{0}, mesh().element().nodes_count()})
                            callback(eL, eNL, iL, jNL);
            }            
    }
}

template<class T, class I>
template<template<class, auto...> class Physical, auto... Args, class Integrate_Loc, class Integrate_Nonloc>
void finite_element_matrix_1d<T, I>::calc_matrix(const std::array<bool, 2> is_first_kind,
                                                 const std::vector<equation_parameters<1, T, Physical, Args...>>& parameters,
                                                 const Integrate_Loc& integrate_rule_loc,
                                                 const Integrate_Nonloc& integrate_rule_nonloc) {
    const auto assemble_bound = [this](std::unordered_map<size_t, T>& matrix_bound, const size_t row, const size_t col, const T integral) {
        if (col == row)
            matrix_inner().coeffRef(row, col) = T{1};
        else if (const auto [it, flag] = matrix_bound.template try_emplace(col, integral); !flag)
            it->second += integral;
    };

    const auto assemble = [this, is_first_kind, &assemble_bound, last_node = mesh().nodes_count() - 1]
                          <class Integrate, class... Model_Args>
                          (const size_t row, const size_t col, const Integrate& integrate_rule, const Model_Args&... args) {
        if (is_first_kind.front() && (row == 0 || col == 0)) {
            if (row == 0)
                assemble_bound(matrix_bound().front(), row, col, integrate_rule(args...));
        } else if (is_first_kind.back() && (row == last_node || col == last_node)) {
            if (row == last_node)
                assemble_bound(matrix_bound().back(), row, col, integrate_rule(args...));
        } else if (row <= col)
            matrix_inner().coeffRef(row, col) += integrate_rule(args...);
    };

    for(const size_t segment : std::ranges::iota_view{size_t{0}, mesh().segments_count()}) {
        const auto& parameter = parameters[segment];
        switch (is_nonlocal(parameters[segment].model.local_weight)) {
            case theory_t::LOCAL:
                mesh_run<theory_t::LOCAL>(segment,
                    [this, integrate_rule_loc, &parameter, &assemble](const size_t e, const size_t i, const size_t j) {
                        assemble(mesh().node_number(e, i), mesh().node_number(e, j), integrate_rule_loc, parameter, e, i, j);
                    }
                );
            break;

            case theory_t::NONLOCAL:
                mesh_run<theory_t::NONLOCAL>(segment,
                    [this, &integrate_rule_nonloc, &parameter, &assemble](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                        assemble(mesh().node_number(eL, iL), mesh().node_number(eNL, jNL), integrate_rule_nonloc, parameter, eL, eNL, iL, jNL);
                    }
                );
            break;
            
            default:
                throw std::runtime_error{"Unknown theory type."};
        }
    }
}

/*
template<class T, class I>
class finite_element_matrix_1d {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> _matrix_inner;
    std::array<std::unordered_map<size_t, T>, 2> _matrix_bound;Ы

protected:
    explicit finite_element_matrix_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh);

    void create_matrix_portrait(const std::array<boundary_condition_t, 2> bound_cond, const theory_t theory);

    template<theory_t Theory, class Callback>
    void mesh_run(const Callback& callback) const;

    template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
    void calc_matrix(const std::array<boundary_condition_t, 2> bound_cond,
                     const theory_t theory,
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

    void clear_matrix();
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
void finite_element_matrix_1d<T, I>::clear_matrix() {
    _matrix_inner = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{};
    _matrix_bound = std::array<std::unordered_map<size_t, T>, 2>{};
}

template<class T, class I>
void finite_element_matrix_1d<T, I>::create_matrix_portrait(const std::array<boundary_condition_t, 2> bound_cond, const theory_t theory) {
#pragma omp parallel for default(none) shared(bound_cond, theory)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        if (bound_cond.front() == boundary_condition_t::FIRST_KIND && node == 0 ||
            bound_cond.back () == boundary_condition_t::FIRST_KIND && node == mesh()->nodes_count()-1)
            matrix_inner().outerIndexPtr()[node+1] = 1;
        else {
            const auto [eL, iL, eR, iR] = mesh()->node_elements(node).named;
            const size_t e = eR == std::numeric_limits<size_t>::max() ? eL : eR,
                         i = iR == std::numeric_limits<size_t>::max() ? iL : iR,
                         right_neighbour = theory == theory_t::NONLOCAL ? mesh()->right_neighbour(e) : e + 1;
            const bool last_node_first_kind = bound_cond.back() == boundary_condition_t::FIRST_KIND &&
                                              right_neighbour * (mesh()->element().nodes_count() - 1) == mesh()->nodes_count()-1;
            matrix_inner().outerIndexPtr()[node+1] += (right_neighbour - e) * (mesh()->element().nodes_count() - 1) - i + 1 - last_node_first_kind;
        }

    utils::prepare_memory(matrix_inner());

    for(const size_t i : std::views::iota(size_t{0}, size_t(matrix_inner().rows())))
        for(size_t j = matrix_inner().outerIndexPtr()[i], k = i; j < matrix_inner().outerIndexPtr()[i+1]; ++j, ++k)
            matrix_inner().innerIndexPtr()[j] = k;
}

template<class T, class I>
template<theory_t Theory, class Callback>
void finite_element_matrix_1d<T, I>::mesh_run(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback) shared(std::views::iota) schedule(dynamic)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        for(const auto& [eL, iL] : mesh()->node_elements(node).arr)
            if(eL != std::numeric_limits<size_t>::max()) {
                if constexpr (Theory == theory_t::LOCAL)
                    for(const size_t jL : std::views::iota(size_t{0}, mesh()->element().nodes_count()))
                        callback(eL, iL, jL);
                if constexpr (Theory == theory_t::NONLOCAL)
                    for(const size_t eNL : std::views::iota(mesh()->left_neighbour(eL), mesh()->right_neighbour(eL)))
                        for(const size_t jNL : std::views::iota(size_t{0}, mesh()->element().nodes_count()))
                            callback(eL, eNL, iL, jNL);
            }
}

template<class T, class I>
template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
void finite_element_matrix_1d<T, I>::calc_matrix(const std::array<boundary_condition_t, 2> bound_cond,
                                                 const theory_t theory,
                                                 const Influence_Function& influence_fun,
                                                 const Integrate_Loc& integrate_rule_loc,
                                                 const Integrate_Nonloc& integrate_rule_nonloc) {
    const auto assemble_bound = [this](std::unordered_map<size_t, T>& matrix_bound, const size_t row, const size_t col, const T integral) {
        if (col == row)
            matrix_inner().coeffRef(row, col) = 1;
        else if (const auto [it, flag] = matrix_bound.template try_emplace(col, integral); !flag)
            it->second += integral;
    };

    const auto assemble = [this, bound_cond, &assemble_bound, last_node = mesh()->nodes_count()-1]<class Integrate, class... Args>
                          (const size_t row, const size_t col, const Integrate& integrate_rule, const Args&... args) {
        if (bound_cond.front() == boundary_condition_t::FIRST_KIND && (row == 0 || col == 0)) {
            if (row == 0)
                assemble_bound(matrix_bound().front(), row, col, integrate_rule(args...));
        } else if (bound_cond.back() == boundary_condition_t::FIRST_KIND && (row == last_node || col == last_node)) {
            if (row == last_node)
                assemble_bound(matrix_bound().back(), row, col, integrate_rule(args...));
        } else if (row <= col)
            matrix_inner().coeffRef(row, col) += integrate_rule(args...);
    };

    mesh_run<theory_t::LOCAL>(
        [this, &integrate_rule_loc, &assemble](const size_t e, const size_t i, const size_t j) {
            assemble(mesh()->node_number(e, i), mesh()->node_number(e, j), integrate_rule_loc, e, i, j);
        }
    );

    if (theory == theory_t::NONLOCAL)
        mesh_run<theory_t::NONLOCAL>(
            [this, &integrate_rule_nonloc, &influence_fun, &assemble]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                assemble(mesh()->node_number(eL, iL), mesh()->node_number(eNL, jNL), integrate_rule_nonloc, eL, eNL, iL, jNL, influence_fun);
            }
        );
}
*/

}

#endif