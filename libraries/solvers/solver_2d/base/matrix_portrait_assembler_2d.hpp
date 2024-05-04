#pragma once

#include "../solvers_utils.hpp"

#include "shift_initializer.hpp"
#include "matrix_index_initializer.hpp"

#include "mesh_2d_utils.hpp"

namespace nonlocal {

template<class T, class I, class J, size_t DoF>
class matrix_portrait_assembler_2d final {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    finite_element_matrix<T, J>& _matrix;
    std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    std::ranges::iota_view<size_t, size_t> _nodes_for_processing;

    size_t cols() const noexcept;
    size_t rows() const noexcept;

    void init_shifts(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, 
                     const std::array<bool, DoF> is_neumann, const bool is_symmetric);
    void init_indices(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, 
                      const std::array<bool, DoF> is_neumann, const bool is_symmetric);

public:
    explicit matrix_portrait_assembler_2d(finite_element_matrix<T, J>& matrix,
                                          const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                          const std::optional<std::ranges::iota_view<size_t, size_t>>& nodes_for_processing = std::nullopt);

    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const noexcept;
    const mesh::mesh_2d<T, I>& mesh() const;
    finite_element_matrix<T, J>& matrix() noexcept;
    const finite_element_matrix<T, J>& matrix() const noexcept;

    void compute(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, 
                 const std::array<bool, DoF> is_neumann = {}, const bool is_symmetric = true);
};

template<class T, class I, class J, size_t DoF>
matrix_portrait_assembler_2d<T, I, J, DoF>::matrix_portrait_assembler_2d(
    finite_element_matrix<T, J>& matrix,
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
    const std::optional<std::ranges::iota_view<size_t, size_t>>& nodes_for_processing)
    : _matrix{matrix}
    , _mesh{mesh}
    , _nodes_for_processing{nodes_for_processing ? *nodes_for_processing : _mesh->process_nodes()} {}

template<class T, class I, class J, size_t DoF>
const std::shared_ptr<mesh::mesh_2d<T, I>>& matrix_portrait_assembler_2d<T, I, J, DoF>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I, class J, size_t DoF>
const mesh::mesh_2d<T, I>& matrix_portrait_assembler_2d<T, I, J, DoF>::mesh() const {
    return *mesh_ptr();
}

template<class T, class I, class J, size_t DoF>
finite_element_matrix<T, J>& matrix_portrait_assembler_2d<T, I, J, DoF>::matrix() noexcept {
    return _matrix;
}

template<class T, class I, class J, size_t DoF>
const finite_element_matrix<T, J>& matrix_portrait_assembler_2d<T, I, J, DoF>::matrix() const noexcept {
    return _matrix;
}

template<class T, class I, class J, size_t DoF>
size_t cols() const noexcept {
    return DoF * mesh().container().nodes_count();
}

template<class T, class I, class J, size_t DoF>
size_t rows() const noexcept {
    return DoF * _nodes_for_processing.size();
}

template<class T, class I, class J, size_t DoF>
void matrix_portrait_assembler_2d<T, I, J, DoF>::init_shifts(
    const std::unordered_map<std::string, theory_t>& theories,
    const std::vector<bool>& is_inner,
    const std::array<bool, DoF> is_neumann,
    const bool is_symmetric) {
    const auto process_rows = std::ranges::iota_view{DoF * _nodes_for_processing.front(), DoF * *_nodes_for_processing.end()};
    mesh_run(theories, shift_initializer<T, J, DoF>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { ++_matrix.inner().outerIndexPtr()[row + 1]; });
    utils::accumulate_shifts(matrix().inner());
    utils::accumulate_shifts(matrix().bound());
    logger::get().log() << "Non-zero elements count: " << 
        matrix().inner().nonZeros() + matrix().bound().nonZeros() << std::endl;
}

template<class T, class I, class J, size_t DoF>
void matrix_portrait_assembler_2d<T, I, J, DoF>::init_indices(
    const std::unordered_map<std::string, theory_t>& theories,
    const std::vector<bool>& is_inner, 
    const std::array<bool, DoF> is_neumann,
    const bool is_symmetric) {
    
}

template<class T, class I, class J, size_t DoF>
void matrix_portrait_assembler_2d<T, I, J, DoF>::compute(const std::unordered_map<std::string, theory_t>& theories, 
                                                         const std::vector<bool>& is_inner, 
                                                         const std::array<bool, DoF> is_neumann,
                                                         const bool is_symmetric) {
    if (std::accumulate(is_neumann.begin(), is_neumann.end(), 0u) > 1u)
        throw std::logical_error{"Methods for solving problems in which Neumann boundary conditions are specified "
                                 "at more than one degree of freedom are currently not available."};
    matrix().clear();
    const bool neumann = std::any_of(is_neumann.begin(), is_neumann.end(), 
        [](const bool is_neumann) constexpr noexcept { return is_neumann; });
    const size_t cols = cols() + neumann;
    const size_t rows = rows() + (neumann && parallel::is_last_process());
    matrix().inner().resize(rows, cols);
    matrix().bound().resize(rows, cols);
    init_shifts(theories, is_inner, is_neumann, is_symmetric);
    init_indices(theories, is_inner, is_neumann, is_symmetric);
}

}