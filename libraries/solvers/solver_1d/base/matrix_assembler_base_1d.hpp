#pragma once

#include "finite_element_matrix_1d.hpp"

#include "../solvers_utils.hpp"
#include "mesh_1d.hpp"

namespace nonlocal {

template<class T, class I>
class matrix_assembler_base_1d {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    finite_element_matrix_1d<T, I>& _matrix;
    utils::nodes_sequence _nodes_for_processing;

    auto correct_node_data(const size_t node, const std::ranges::iota_view<size_t, size_t> segment_nodes) const;
    template<class Nodes, class Local_Runner, class Nonlocal_Runner>
    void mesh_run(const Nodes& nodes, const std::vector<theory_t>& theories, 
                  Local_Runner&& local_runner, Nonlocal_Runner&& nonlocal_runner) const;
    template<class Local_Runner, class Nonlocal_Runner>
    void mesh_run(const std::vector<theory_t>& theories, Local_Runner&& local_runner, Nonlocal_Runner&& nonlocal_runner) const;

    T integrate_basic(const size_t e, const size_t i) const;
    template<class Nodes>
    void integral_condition(const Nodes& nodes, const bool is_symmetric);
    void integral_condition(const bool is_symmetric); // for Neumann problem

protected:
    explicit matrix_assembler_base_1d(finite_element_matrix_1d<T, I>& matrix,
                                      const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                      const std::optional<utils::nodes_sequence>& nodes_for_processing = std::nullopt);
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
    const utils::nodes_sequence& nodes_for_processing() const noexcept;
    void nodes_for_processing(const utils::nodes_sequence& nodes);
};

template<class T, class I>
matrix_assembler_base_1d<T, I>::matrix_assembler_base_1d(finite_element_matrix_1d<T, I>& matrix,
                                                         const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                         const std::optional<utils::nodes_sequence>& nodes_for_processing)
    : _mesh{mesh}
    , _matrix{matrix}
    , _nodes_for_processing{nodes_for_processing ? *nodes_for_processing : 
                            std::ranges::iota_view<size_t, size_t>{0u, _mesh->nodes_count()}}
    {}

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
const utils::nodes_sequence& matrix_assembler_base_1d<T, I>::nodes_for_processing() const noexcept {
    return _nodes_for_processing;
}

template<class T, class I>
void matrix_assembler_base_1d<T, I>::nodes_for_processing(const utils::nodes_sequence& nodes) {
    _nodes_for_processing = nodes;
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
template<class Nodes, class Local_Runner, class Nonlocal_Runner>
void matrix_assembler_base_1d<T, I>::mesh_run(const Nodes& nodes, const std::vector<theory_t>& theories,
                                              Local_Runner&& local_runner, Nonlocal_Runner&& nonlocal_runner) const {
#pragma omp parallel for default(none) shared(nodes, theories) firstprivate(local_runner, nonlocal_runner)
    for(size_t i = 0; i < nodes.size(); ++i) {
        const size_t node = nodes[i];
        for(const auto node_data : mesh().node_elements(node).to_array())
            if (node_data) {
                const auto& [eL, iL] = node_data;
                const size_t segment = mesh().segment_number(eL);
                switch (const theory_t theory = theories[segment]) {
                case theory_t::NONLOCAL:
                    for(const size_t eNL : mesh().neighbours(eL))
                        for(const size_t jNL : mesh().element().nodes())
                            nonlocal_runner(segment, eL, eNL, iL, jNL);
                case theory_t::LOCAL:
                    for(const size_t jL : mesh().element().nodes())
                        local_runner(segment, eL, iL, jL);
                break;
                default:
                    throw std::domain_error{"Unknown theory."};
                }
            }
    }
}

template<class T, class I>
template<class Local_Runner, class Nonlocal_Runner>
void matrix_assembler_base_1d<T, I>::mesh_run(const std::vector<theory_t>& theories, 
                                              Local_Runner&& local_runner,
                                              Nonlocal_Runner&& nonlocal_runner) const {
    if (std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(nodes_for_processing()))
        mesh_run(std::get<std::ranges::iota_view<size_t, size_t>>(nodes_for_processing()), theories, 
                 std::forward<Local_Runner>(local_runner), std::forward<Nonlocal_Runner>(nonlocal_runner));
    else
        mesh_run(std::get<std::vector<size_t>>(nodes_for_processing()), theories, 
                 std::forward<Local_Runner>(local_runner), std::forward<Nonlocal_Runner>(nonlocal_runner));
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
template<class Nodes>
void matrix_assembler_base_1d<T, I>::integral_condition(const Nodes& nodes, const bool is_symmetric) {
#pragma omp parallel for default(none) shared(nodes, is_symmetric)
    for(size_t k = 0; k < nodes.size(); ++k) {
        const size_t node = nodes[k];
        T& val = matrix().inner().coeffRef(node, mesh().nodes_count());
        for(const auto& [e, i] : mesh().node_elements(node).to_array())
            if (e) val += integrate_basic(e, i);
        if (!is_symmetric)
            matrix().inner().coeffRef(mesh().nodes_count(), node) = val;
    }
}

template<class T, class I>
void matrix_assembler_base_1d<T, I>::integral_condition(const bool is_symmetric) {
    if (std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(nodes_for_processing()))
        integral_condition(std::get<std::ranges::iota_view<size_t, size_t>>(nodes_for_processing()), is_symmetric);
    else
        integral_condition(std::get<std::vector<size_t>>(nodes_for_processing()), is_symmetric);
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
    utils::clear_matrix_rows(matrix().inner(), nodes_for_processing());
    if (is_neumann)
        integral_condition(is_symmetric);

    mesh_run(theories,
        [this, &assemble, &integrate_rule_loc](const size_t segment, const size_t e, const size_t i, const size_t j) {
            assemble(mesh().node_number(e, i), mesh().node_number(e, j), integrate_rule_loc, segment, e, i, j);
        },
        [this, &assemble, &integrate_rule_nonloc](const size_t segment, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            assemble(mesh().node_number(eL, iL), mesh().node_number(eNL, jNL), integrate_rule_nonloc, segment, eL, eNL, iL, jNL);
        }
    );
}

}