#ifndef NONLOCAL_MATRIX_ASSEMBLER_2D_HPP
#define NONLOCAL_MATRIX_ASSEMBLER_2D_HPP

#include "../solvers_utils.hpp"

#include "shift_initializer.hpp"
#include "matrix_index_initializer.hpp"
#include "integrator.hpp"

#include "mesh_2d_utils.hpp"

#include <variant>
#include <iostream>

namespace nonlocal {

enum class assemble_part : uint8_t {
    LOCAL,
    NONLOCAL,
    FULL
};

template<class T, class I>
std::unordered_map<std::string, theory_t> local_theories(const mesh::mesh_container_2d<T, I>& mesh) {
    const auto theroires_setter = std::views::all(mesh.groups_2d()) |
                                  std::views::transform([](const std::string& group) { return std::pair{group, theory_t::LOCAL}; });
    const std::unordered_map<std::string, theory_t> theories(theroires_setter.begin(), theroires_setter.end());
    return theories;
}

template<class Callback>
void first_kind_filler(const std::ranges::iota_view<size_t, size_t> rows, 
                       const std::vector<bool>& is_inner, const Callback& callback) {
    for(const size_t row : rows)
        if (!is_inner[row])
            callback(row);
}

template<class T, class I, class J, size_t DoF>
class matrix_assembler_2d {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    using nodes_sequence = std::variant<
        std::ranges::iota_view<size_t, size_t>,
        std::vector<size_t>
    >;

    finite_element_matrix<T, J> _matrix;
    std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    nodes_sequence _nodes_for_processing;

protected:
    explicit matrix_assembler_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);

    template<class Runner>
    void mesh_run(const std::unordered_map<std::string, theory_t>& theories, Runner&& runner);

    void init_shifts(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric);
    void init_indices(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner,
                      const bool is_symmetric, const bool sort_indices = true);
    template<class Local_Integrator, class Nonlocal_Integrator>
    void calc_coeffs(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric,
                     Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator);

public:
    virtual ~matrix_assembler_2d() noexcept = default;

    size_t cols() const noexcept;
    size_t rows() const noexcept;

    const mesh::mesh_2d<T, I>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const noexcept;
    finite_element_matrix<T, J>& matrix() noexcept;
    const finite_element_matrix<T, J>& matrix() const noexcept;
    const nodes_sequence& nodes_for_processing() const noexcept;

    void nodes_for_processing(const nodes_sequence& nodes);

    void clear();
};

template<class T, class I, class J, size_t DoF>
matrix_assembler_2d<T, I, J, DoF>::matrix_assembler_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _mesh{mesh}
    , _nodes_for_processing{mesh->process_nodes()} {}

template<class T, class I, class J, size_t DoF>
size_t matrix_assembler_2d<T, I, J, DoF>::cols() const noexcept {
    if (std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing))
        return DoF * mesh().container().nodes_count();
    return matrix().inner().cols();
}

template<class T, class I, class J, size_t DoF>
size_t matrix_assembler_2d<T, I, J, DoF>::rows() const noexcept {
    if (std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing))
        return DoF * std::get<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing).size();
    return matrix().inner().rows();
}

template<class T, class I, class J, size_t DoF>
const mesh::mesh_2d<T, I>& matrix_assembler_2d<T, I, J, DoF>::mesh() const {
    return *mesh_ptr();
}

template<class T, class I, class J, size_t DoF>
const std::shared_ptr<mesh::mesh_2d<T, I>>& matrix_assembler_2d<T, I, J, DoF>::mesh_ptr() const noexcept {
    return _mesh;
}

template<class T, class I, class J, size_t DoF>
finite_element_matrix<T, J>& matrix_assembler_2d<T, I, J, DoF>::matrix() noexcept {
    return _matrix;
}

template<class T, class I, class J, size_t DoF>
const finite_element_matrix<T, J>& matrix_assembler_2d<T, I, J, DoF>::matrix() const noexcept {
    return _matrix;
}

template<class T, class I, class J, size_t DoF>
const matrix_assembler_2d<T, I, J, DoF>::nodes_sequence& matrix_assembler_2d<T, I, J, DoF>::nodes_for_processing() const noexcept {
    return _nodes_for_processing;
}

template<class T, class I, class J, size_t DoF>
void matrix_assembler_2d<T, I, J, DoF>::nodes_for_processing(const nodes_sequence& nodes) {
    _nodes_for_processing = nodes;
}

template<class T, class I, class J, size_t DoF>
void matrix_assembler_2d<T, I, J, DoF>::clear() {
    _matrix.inner() = Eigen::SparseMatrix<T, Eigen::RowMajor, J>{};
    _matrix.bound() = Eigen::SparseMatrix<T, Eigen::RowMajor, J>{};
}

template<class T, class I, class J, size_t DoF>
template<class Runner>
void matrix_assembler_2d<T, I, J, DoF>::mesh_run(const std::unordered_map<std::string, theory_t>& theories, Runner&& runner) {
    if (std::holds_alternative<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing))
        mesh::utils::mesh_run(mesh(), std::get<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing), theories, std::forward<Runner>(runner));
    else
        mesh::utils::mesh_run(mesh(), std::get<std::vector<size_t>>(_nodes_for_processing), theories, std::forward<Runner>(runner));
}

template<class T, class I, class J, size_t DoF>
void matrix_assembler_2d<T, I, J, DoF>::init_shifts(
    const std::unordered_map<std::string, theory_t>& theories, 
    const std::vector<bool>& is_inner, const bool is_symmetric) {
    const auto process_nodes = std::get<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing);
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, shift_initializer<T, J, DoF>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { ++_matrix.inner().outerIndexPtr()[row + 1]; });
    utils::accumulate_shifts(matrix().inner());
    utils::accumulate_shifts(matrix().bound());
    logger::get().log() << "Non-zero elements count: " << 
        matrix().inner().nonZeros() + matrix().bound().nonZeros() << std::endl;
}

template<class T, class I, class J, size_t DoF>
void matrix_assembler_2d<T, I, J, DoF>::init_indices(
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, 
    const bool is_symmetric, const bool sort_indices) {
    const auto process_nodes = std::get<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing);
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    utils::allocate_matrix(matrix().inner());
    utils::allocate_matrix(matrix().bound());
    logger::get().log() << "Matrix allocated successfully" << std::endl;
    mesh_run(theories, matrix_index_initializer<T, J, DoF>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this, shift = process_rows.front()](const size_t row) { 
        matrix().inner().innerIndexPtr()[matrix().inner().outerIndexPtr()[row]] = row + shift; 
    });
    if (sort_indices) {
        utils::sort_indices(matrix().inner());
        utils::sort_indices(matrix().bound());
    }
}

template<class T, class I, class J, size_t DoF>
template<class Local_Integrator, class Nonlocal_Integrator>
void matrix_assembler_2d<T, I, J, DoF>::calc_coeffs(
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric,
    Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator) {
    const auto process_nodes = std::get<std::ranges::iota_view<size_t, size_t>>(_nodes_for_processing);
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, integrator<T, J, DoF, Local_Integrator, Nonlocal_Integrator>{
        _matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric, local_integrator, nonlocal_integrator});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { 
        matrix().inner().valuePtr()[matrix().inner().outerIndexPtr()[row]] = T{1};
    });
}

}

#endif