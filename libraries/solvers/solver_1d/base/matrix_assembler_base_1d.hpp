#pragma once

#include "finite_element_matrix_1d.hpp"

#include "../solvers_utils.hpp"
#include "mesh_1d.hpp"

namespace nonlocal {

template<class T, class I>
class matrix_assembler_base_1d {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    finite_element_matrix_1d<T, I>& _matrix;

    auto correct_node_data(const size_t node, const std::ranges::iota_view<size_t, size_t> segment_nodes) const;
    template<theory_t Theory, class Callback>
    void mesh_run(const size_t segment, const Callback& callback) const;

    T integrate_basic(const size_t e, const size_t i) const;
    void integral_condition(const bool is_symmetric); // for Neumann problem

protected:
    explicit matrix_assembler_base_1d(finite_element_matrix_1d<T, I>& matrix,
                                      const std::shared_ptr<mesh::mesh_1d<T>>& mesh);
    virtual ~matrix_assembler_base_1d() = default;

    template<class Integrate_Loc, class Integrate_Nonloc>
    void compute(const std::vector<theory_t>& theories, const std::array<bool, 2> is_first_kind,
                 const bool is_symmetric, const bool is_neumann,
                 const Integrate_Loc& integrate_rule_loc, const Integrate_Nonloc& integrate_rule_nonloc);

public:
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    finite_element_matrix_1d<T, I>& matrix() noexcept;
    const finite_element_matrix_1d<T, I>& matrix() const noexcept;
};

template<class T, class I>
matrix_assembler_base_1d<T, I>::matrix_assembler_base_1d(finite_element_matrix_1d<T, I>& matrix,
                                                         const std::shared_ptr<mesh::mesh_1d<T>>& mesh)
    : _mesh{mesh}
    , _matrix{matrix} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& matrix_assembler_base_1d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
const mesh::mesh_1d<T>& matrix_assembler_base_1d<T, I>::mesh() const {
    return *mesh_ptr();
}

template<class T, class I>
finite_element_matrix_1d<T, I>& matrix_assembler_base_1d<T, I>::matrix() noexcept {
    return _matrix;
}

template<class T, class I>
const finite_element_matrix_1d<T, I>& matrix_assembler_base_1d<T, I>::matrix() const noexcept {
    return _matrix;
}

template<class T, class I>
auto matrix_assembler_base_1d<T, I>::correct_node_data(const size_t node, const std::ranges::iota_view<size_t, size_t> segment_nodes) const {
    auto node_elements = mesh().node_elements(node).to_array();
    if (node == segment_nodes.front() && node != 0)
        node_elements.front() = std::nullopt;
    else if (node == segment_nodes.back() && node != mesh().nodes_count() - 1)
        node_elements.back() = std::nullopt;
    return node_elements;
}

template<class T, class I>
template<theory_t Theory, class Callback>
void matrix_assembler_base_1d<T, I>::mesh_run(const size_t segment, const Callback& callback) const {
    const auto segment_nodes = mesh().nodes(segment);
#pragma omp parallel for default(none) shared(segment_nodes) firstprivate(callback)
    for(size_t node = segment_nodes.front(); node < *segment_nodes.end(); ++node) {
        for(const auto node_data : correct_node_data(node, segment_nodes))
            if (node_data) {
                const auto& [eL, iL] = node_data;
                if constexpr (Theory == theory_t::LOCAL)
                    for(const size_t jL : mesh().element().nodes())
                        callback(eL, iL, jL);
                if constexpr (Theory == theory_t::NONLOCAL)
                    for(const size_t eNL : mesh().neighbours(eL))
                        for(const size_t jNL : mesh().element().nodes())
                            callback(eL, eNL, iL, jNL);
            }            
    }
}

template<class T, class I>
T matrix_assembler_base_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = mesh().element();
    for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
        integral += el.weight(q) * el.qN(i, q);
    return integral * mesh().jacobian(mesh().segment_number(e));
}

template<class T, class I>
void matrix_assembler_base_1d<T, I>::integral_condition(const bool is_symmetric) {
#pragma omp parallel for default(none) shared(is_symmetric)
    for(size_t node = 0; node < mesh().nodes_count(); ++node) {
        T& val = matrix().inner().coeffRef(node, mesh().nodes_count());
        for(const auto& [e, i] : mesh().node_elements(node).to_array())
            if (e) val += integrate_basic(e, i);
        if (!is_symmetric)
            matrix().inner().coeffRef(mesh().nodes_count(), node) = val;
    }
}

template<class T, class I>
template<class Integrate_Loc, class Integrate_Nonloc>
void matrix_assembler_base_1d<T, I>::compute(const std::vector<theory_t>& theories, const std::array<bool, 2> is_first_kind,
                                             const bool is_symmetric, const bool is_neumann,
                                             const Integrate_Loc& integrate_rule_loc, const Integrate_Nonloc& integrate_rule_nonloc) {
    const auto assemble_bound = [this](std::unordered_map<size_t, T>& matrix_bound, const size_t row, const size_t col, const T integral) {
        if (col == row)
            matrix().inner().coeffRef(row, col) = T{1};
        else if (const auto [it, flag] = matrix_bound.template try_emplace(col, integral); !flag)
            it->second += integral;
    };

    const auto assemble = [this, is_first_kind, is_symmetric, &assemble_bound, last_node = mesh().nodes_count() - 1]
                          <class Integrate, class... Model_Args>
                          (const size_t row, const size_t col, const Integrate& integrate_rule, const Model_Args&... args) {
        if (is_first_kind.front() && (row == 0 || col == 0)) {
            if (row == 0)
                assemble_bound(matrix().bound().front(), row, col, integrate_rule(args...));
        } else if (is_first_kind.back() && (row == last_node || col == last_node)) {
            if (row == last_node)
                assemble_bound(matrix().bound().back(), row, col, integrate_rule(args...));
        } else if (!is_symmetric || row <= col)
            matrix().inner().coeffRef(row, col) += integrate_rule(args...);
    };

    matrix().bound() = {};
    utils::clear_matrix_coefficients(matrix().inner());
    if (is_neumann)
        integral_condition(is_symmetric);

    for(const size_t segment : mesh().segments())
        switch (theories[segment]) {
            case theory_t::NONLOCAL:
                mesh_run<theory_t::NONLOCAL>(segment,
                    [this, &assemble, &integrate_rule_nonloc, segment](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                        assemble(mesh().node_number(eL, iL), mesh().node_number(eNL, jNL), integrate_rule_nonloc, segment, eL, eNL, iL, jNL);
                    }
                );
            case theory_t::LOCAL:
                mesh_run<theory_t::LOCAL>(segment,
                    [this, &assemble, integrate_rule_loc, segment](const size_t e, const size_t i, const size_t j) {
                        assemble(mesh().node_number(e, i), mesh().node_number(e, j), integrate_rule_loc, segment, e, i, j);
                    }
                );
            break;
            
            default:
                throw std::runtime_error{"Unknown theory type."};
        }
}

}