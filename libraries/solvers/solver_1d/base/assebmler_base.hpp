#pragma once

#include "finite_element_matrix.hpp"
#include "problem_settings.hpp"

#include <logger/logger.hpp>
#include <mesh/mesh_1d/mesh_1d.hpp>
#include <parallel/OMP_utils.hpp>
#include <solvers/base/utils.hpp>

namespace nonlocal {

template<class T, class I>
class assembler_base_1d {
    static constexpr auto Empty_Function = [](const size_t, const size_t, const size_t, const size_t, const size_t) constexpr noexcept { return T{0}; };

    finite_element_matrix_1d<T, I>& _matrix;
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    utils::nodes_sequence _nodes_to_assemble;

    std::vector<size_t> get_subsequence(const size_t segment) const;

protected:
    explicit assembler_base_1d(finite_element_matrix_1d<T, I>& matrix, const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                               const std::optional<utils::nodes_sequence>& nodes_to_assemble = std::nullopt);

    template<theory_t Theory, class Callback>
    void mesh_run(const size_t segment, const Callback& callback) const;
    template<class Integrate_Loc, class Integrate_Nonloc = decltype(Empty_Function)>
    void calc_matrix(const problem_settings& settings, const Integrate_Loc& integrate_loc,
                     const Integrate_Nonloc& integrate_nonloc = Empty_Function);

public:
    virtual ~assembler_base_1d() = default;

    const mesh::mesh_1d<T>& mesh() const;
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    finite_element_matrix_1d<T, I>& matrix() noexcept;
    const finite_element_matrix_1d<T, I>& matrix() const noexcept;
};

template<class T, class I>
assembler_base_1d<T, I>::assembler_base_1d(finite_element_matrix_1d<T, I>& matrix, 
                                           const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                           const std::optional<utils::nodes_sequence>& nodes_to_assemble)
    : _matrix{matrix}
    , _mesh{mesh}
    , _nodes_to_assemble{nodes_to_assemble ? *nodes_to_assemble : std::ranges::iota_view<size_t, size_t>{0u, mesh->nodes_count()}}
    {}

template<class T, class I>
const mesh::mesh_1d<T>& assembler_base_1d<T, I>::mesh() const {
    return *_mesh;
}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& assembler_base_1d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
finite_element_matrix_1d<T, I>& assembler_base_1d<T, I>::matrix() noexcept {
    return _matrix;
}

template<class T, class I>
const finite_element_matrix_1d<T, I>& assembler_base_1d<T, I>::matrix() const noexcept {
    return _matrix;
}

template<class T, class I>
std::vector<size_t> assembler_base_1d<T, I>::get_subsequence(const size_t segment) const {
    std::vector<size_t> subsequence;
    const auto segment_nodes = mesh().nodes(segment);
    subsequence.reserve(segment_nodes.size());
    utils::iterate(_nodes_to_assemble, [&subsequence, &segment_nodes](const size_t row) {
        if (row >= segment_nodes.front() && row <= segment_nodes.back())
            subsequence.push_back(row);
    });
    subsequence.shrink_to_fit();
    return subsequence;
}

template<class T, class I>
template<theory_t Theory, class Callback>
void assembler_base_1d<T, I>::mesh_run(const size_t segment, const Callback& callback) const {
    const auto correct_node_data = [this, segment_nodes = mesh().nodes(segment)](const size_t node) {
        auto node_elements = mesh().node_elements(node).to_array();
        if (node == segment_nodes.front() && node != 0)
            node_elements.front() = std::nullopt;
        else if (node == segment_nodes.back() && node != mesh().nodes_count() - 1)
            node_elements.back() = std::nullopt;
        return node_elements;
    };

    utils::iterate(get_subsequence(segment), [this, &correct_node_data, &callback](const size_t node) {
        for(const auto node_data : correct_node_data(node))
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
    }, parallel::threads_count());
}

template<class T, class I>
template<class Integrate_Loc, class Integrate_Nonloc>
void assembler_base_1d<T, I>::calc_matrix(const problem_settings& settings,
                                          const Integrate_Loc& integrate_loc, 
                                          const Integrate_Nonloc& integrate_nonloc) {
    const auto assemble_bound = [this](std::unordered_map<size_t, T>& matrix_bound, const size_t row, const size_t col, const T integral) {
        if (col == row)
            matrix().inner.coeffRef(row, col) = T{1};
        else if (const auto [it, flag] = matrix_bound.template try_emplace(col, integral); !flag)
            it->second += integral;
    };

    const auto assemble = [this, &settings, &assemble_bound, is_symmetric = settings.is_symmetric(), last_node = mesh().nodes_count() - 1]
                          <class Integrate, class... Model_Args>
                          (const size_t row, const size_t col, const Integrate& integrate, const Model_Args&... args) {
        if (settings.is_first_kind.front() && (row == 0 || col == 0)) {
            if (row == 0)
                assemble_bound(matrix().bound.front(), row, col, integrate(args...));
        } else if (settings.is_first_kind.back() && (row == last_node || col == last_node)) {
            if (row == last_node)
                assemble_bound(matrix().bound.back(), row, col, integrate(args...));
        } else if (!is_symmetric || row <= col)
            matrix().inner.coeffRef(row, col) += integrate(args...);
    };

    for(const size_t segment : mesh().segments())
        switch (settings.theories[segment]) {
            case theory_t::NONLOCAL:
                mesh_run<theory_t::NONLOCAL>(segment,
                    [this, &assemble, &integrate_nonloc, segment](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                        assemble(mesh().node_number(eL, iL), mesh().node_number(eNL, jNL), integrate_nonloc, segment, eL, eNL, iL, jNL);
                    }
                );
            case theory_t::LOCAL:
                mesh_run<theory_t::LOCAL>(segment,
                    [this, &assemble, integrate_loc, segment](const size_t e, const size_t i, const size_t j) {
                        assemble(mesh().node_number(e, i), mesh().node_number(e, j), integrate_loc, segment, e, i, j);
                    }
                );
            break;
            
            default:
                throw std::runtime_error{"Unknown theory type."};
        }
}

}