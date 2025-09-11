#pragma once

#include "problem_settings.hpp"

#include <mesh/mesh_1d/mesh_1d.hpp>
#include <solvers/base/utils.hpp>

#include <array>
#include <vector>

namespace nonlocal::solver_1d {

// Supported only upper part
template<class T, class I>
void init_shifts(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                 const mesh::mesh_1d<T>& mesh,
                 const problem_settings& settings) {
    if (settings.theories.size() != mesh.segments_count())
        throw std::domain_error{"The number of segments and the theories number do not match."};
    const size_t matrix_size = mesh.nodes_count() + settings.is_neumann;
    matrix.resize(matrix_size, matrix_size);
    const size_t last_node = mesh.nodes_count() - 1;
    const size_t nodes_in_element = mesh.element().nodes_count() - 1;
    if (settings.is_first_kind.front())
        matrix.outerIndexPtr()[1] = 1;
#pragma omp parallel for default(none) shared(matrix, mesh, settings, last_node, nodes_in_element)
    for(size_t node = settings.is_first_kind.front(); node < mesh.nodes_count() - settings.is_first_kind.back(); ++node) {
        const auto [left, right] = mesh.node_elements(node);
        const auto [e, i] = right ? right : left;
        const size_t neighbour = settings.theories[mesh.segment_number(e)] == theory_t::LOCAL ? e + 1 : *mesh.neighbours(e).end();
        const bool is_right_first_kind = settings.is_first_kind.back() && neighbour * nodes_in_element == last_node;
        const size_t count = (neighbour - e) * nodes_in_element + 1 - i - is_right_first_kind;
        matrix.outerIndexPtr()[node + 1] = count + settings.is_neumann;
    }
    if (settings.is_first_kind.back())
        matrix.outerIndexPtr()[mesh.nodes_count()] = 1;
    utils::accumulate_shifts(matrix);
    logger::info() << "Non-zero elements in matrix: " << matrix.nonZeros() << std::endl;
}

// Supported only upper part
template<class T, class I>
void init_indices(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, const bool is_neumann = false) {
    utils::allocate_matrix(matrix);
    for(const size_t i : std::ranges::iota_view{0u, size_t(matrix.rows())})
        for(size_t j = matrix.outerIndexPtr()[i], k = i; j < matrix.outerIndexPtr()[i+1]; ++j)
            matrix.innerIndexPtr()[j] = k++;
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, size_t(matrix.rows())})
            matrix.innerIndexPtr()[matrix.outerIndexPtr()[row + 1] - 1] = matrix.cols() - 1;
}

template<class T, class I>
void init_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                          const mesh::mesh_1d<T>& mesh,
                          const problem_settings& settings) {
    init_shifts(matrix, mesh, settings);
    init_indices(matrix, settings.is_neumann);
    if (!settings.is_symmetric())
        matrix = Eigen::SparseMatrix<T, Eigen::RowMajor, I>{matrix.template selfadjointView<Eigen::Upper>()};
}

}