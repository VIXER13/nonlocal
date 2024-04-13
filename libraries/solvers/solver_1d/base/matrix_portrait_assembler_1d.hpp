#pragma once

#include "finite_element_matrix_1d.hpp"

#include "../solvers_utils.hpp"
#include "mesh_1d.hpp"

namespace nonlocal {

template<class T, class I>
class matrix_portrait_assembler_1d final {
    std::shared_ptr<mesh::mesh_1d<T>> _mesh;
    finite_element_matrix_1d<T, I>& _matrix;
    utils::nodes_sequence _nodes_for_processing;

    void calc_shifts(const std::vector<theory_t>& theories, const std::array<bool, 2> is_first_kind, const bool is_neumann);
    void init_indices(const bool is_neumann);

public:
    explicit matrix_portrait_assembler_1d(finite_element_matrix_1d<T, I>& matrix,
                                          const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                          const std::optional<utils::nodes_sequence>& nodes_for_processing = std::nullopt);

    const std::shared_ptr<mesh::mesh_1d<T>>& mesh_ptr() const noexcept;
    const mesh::mesh_1d<T>& mesh() const;
    finite_element_matrix_1d<T, I>& matrix() noexcept;
    const finite_element_matrix_1d<T, I>& matrix() const noexcept;

    void compute(const std::vector<theory_t>& theories, const std::array<bool, 2> is_first_kind, 
                 const bool is_symmetric = true, const bool is_neumann = false);
};

template<class T, class I>
matrix_portrait_assembler_1d<T, I>::matrix_portrait_assembler_1d(finite_element_matrix_1d<T, I>& matrix,
                                                                 const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                 const std::optional<utils::nodes_sequence>& nodes_for_processing)
    : _mesh{mesh}
    , _matrix{matrix}
    , _nodes_for_processing{nodes_for_processing ? *nodes_for_processing : 
                            std::ranges::iota_view<size_t, size_t>{0u, _mesh->nodes_count()}}
    {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T>>& matrix_portrait_assembler_1d<T, I>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I>
const mesh::mesh_1d<T>& matrix_portrait_assembler_1d<T, I>::mesh() const {
    return *mesh_ptr();
}

template<class T, class I>
finite_element_matrix_1d<T, I>& matrix_portrait_assembler_1d<T, I>::matrix() noexcept {
    return _matrix;
}

template<class T, class I>
const finite_element_matrix_1d<T, I>& matrix_portrait_assembler_1d<T, I>::matrix() const noexcept {
    return _matrix;
}

template<class T, class I>
void matrix_portrait_assembler_1d<T, I>::calc_shifts(const std::vector<theory_t>& theories, 
                                                     const std::array<bool, 2> is_first_kind,
                                                     const bool is_neumann) {
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, size_t(matrix().inner().rows())})
            matrix().inner().outerIndexPtr()[row + 1] = 1;
    const size_t last_node = mesh().nodes_count() - 1;
    const size_t nodes_in_element = mesh().element().nodes_count() - 1;
    if (is_first_kind.front())
        matrix().inner().outerIndexPtr()[1] = 1;
#pragma omp parallel for default(none) shared(theories, is_first_kind, last_node, nodes_in_element)
    for(size_t node = is_first_kind.front(); node < mesh().nodes_count() - is_first_kind.back(); ++node) {
        const auto [left, right] = mesh().node_elements(node);
        const auto [e, i] = right ? right : left;
        const size_t neighbour = theories[mesh().segment_number(e)] == theory_t::LOCAL ? e + 1 : *mesh().neighbours(e).end();
        const bool is_last_first_kind = is_first_kind.back() && neighbour * nodes_in_element == last_node;
        const size_t count = (neighbour - e) * nodes_in_element + 1 - i - is_last_first_kind;
        matrix().inner().outerIndexPtr()[node + 1] += count;
    }
    if (is_first_kind.back())
        matrix().inner().outerIndexPtr()[mesh().nodes_count()] = 1;
}

template<class T, class I>
void matrix_portrait_assembler_1d<T, I>::init_indices(const bool is_neumann) {
    for(const size_t i : std::ranges::iota_view{0u, size_t(matrix().inner().rows())})
        for(size_t j = matrix().inner().outerIndexPtr()[i], k = i; j < matrix().inner().outerIndexPtr()[i+1]; ++j)
            matrix().inner().innerIndexPtr()[j] = k++;
    if (is_neumann)
        for(const size_t row : std::ranges::iota_view{0u, size_t(matrix().inner().rows())})
            matrix().inner().innerIndexPtr()[matrix().inner().outerIndexPtr()[row + 1] - 1] = mesh().nodes_count();
}

template<class T, class I>
void matrix_portrait_assembler_1d<T, I>::compute(const std::vector<theory_t>& theories, const std::array<bool, 2> is_first_kind,
                                                 const bool is_symmetric, const bool is_neumann) {
    matrix().clear();
    const size_t matrix_size = mesh().nodes_count() + is_neumann;
    matrix().inner().resize(matrix_size, matrix_size);
    calc_shifts(theories, is_first_kind, is_neumann);
    utils::accumulate_shifts(matrix().inner());
    logger::get().log() << "Non-zero elements count: " << matrix().inner().nonZeros() << std::endl;
    utils::allocate_matrix(matrix().inner());
    init_indices(is_neumann);
    if (!is_symmetric) // TODO: optimize code for nonsymmetric matrices
        matrix().inner() = Eigen::SparseMatrix<T, Eigen::RowMajor, I>(matrix().inner().template selfadjointView<Eigen::Upper>());
}

}