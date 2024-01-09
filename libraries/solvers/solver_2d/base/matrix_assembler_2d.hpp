#ifndef NONLOCAL_MATRIX_ASSEMBLER_2D_HPP
#define NONLOCAL_MATRIX_ASSEMBLER_2D_HPP

#include "../solvers_utils.hpp"

#include "shift_initializer.hpp"
#include "index_initializer.hpp"
#include "integrator.hpp"

#include "mesh_2d.hpp"

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

    std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    finite_element_matrix<T, J> _matrix;

protected:
    explicit matrix_assembler_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);

    template<class Initializer>
    void mesh_run(const std::unordered_map<std::string, theory_t>& theories, Initializer&& initializer);

    void init_shifts(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric);
    void init_indices(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner,
                      const bool is_symmetric, const bool sort_indices = true);
    template<class Local_Integrator, class Nonlocal_Integrator>
    void calc_coeffs(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric,
                     Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator);

public:
    virtual ~matrix_assembler_2d() noexcept = default;

    const mesh::mesh_2d<T, I>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const noexcept;
    finite_element_matrix<T, J>& matrix() noexcept;
    const finite_element_matrix<T, J>& matrix() const noexcept;

    void clear();
};

template<class T, class I, class J, size_t DoF>
matrix_assembler_2d<T, I, J, DoF>::matrix_assembler_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _mesh{mesh} {}

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
void matrix_assembler_2d<T, I, J, DoF>::clear() {
    _matrix[matrix_part::INNER] = Eigen::SparseMatrix<T, Eigen::RowMajor, J>{};
    _matrix[matrix_part::BOUND] = Eigen::SparseMatrix<T, Eigen::RowMajor, J>{};
}

template<class T, class I, class J, size_t DoF>
template<class Initializer>
void matrix_assembler_2d<T, I, J, DoF>::mesh_run(const std::unordered_map<std::string, theory_t>& theories,
                                                 Initializer&& initializer) {
    const auto process_nodes = mesh().process_nodes();
#pragma omp parallel for default(none) shared(theories, process_nodes) firstprivate(initializer) schedule(dynamic)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        if constexpr (std::is_base_of_v<indexator_base<DoF>, Initializer>)
            initializer.reset(node);
        for(const I eL : mesh().elements(node)) {
            const size_t iL = mesh().global_to_local(eL, node);
            const std::string& group = mesh().container().group(eL);
            if (const theory_t theory = theories.at(group); theory == theory_t::LOCAL)
                for(const size_t jL : std::ranges::iota_view{0u, mesh().container().nodes_count(eL)})
                    initializer(group, eL, iL, jL);
            else if (theory == theory_t::NONLOCAL)
                for(const I eNL : mesh().neighbours(eL))
                    for(const size_t jNL : std::ranges::iota_view{0u, mesh().container().nodes_count(eNL)})
                        initializer(group, eL, eNL, iL, jNL);
            else
                throw std::domain_error{"Unknown theory."};
        }
    }
}

template<class T, class I, class J, size_t DoF>
void matrix_assembler_2d<T, I, J, DoF>::init_shifts(
    const std::unordered_map<std::string, theory_t>& theories, 
    const std::vector<bool>& is_inner, const bool is_symmetric) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, shift_initializer<T, J, DoF>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { ++_matrix[matrix_part::INNER].outerIndexPtr()[row + 1]; });
    utils::accumulate_shifts(matrix()[matrix_part::INNER]);
    utils::accumulate_shifts(matrix()[matrix_part::BOUND]);
    std::cout << "Non-zero elements count: " << matrix()[matrix_part::INNER].nonZeros() + matrix()[matrix_part::BOUND].nonZeros() << std::endl;
}

template<class T, class I, class J, size_t DoF>
void matrix_assembler_2d<T, I, J, DoF>::init_indices(
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, 
    const bool is_symmetric, const bool sort_indices) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    utils::allocate_matrix(matrix()[matrix_part::INNER]);
    utils::allocate_matrix(matrix()[matrix_part::BOUND]);
    mesh_run(theories, index_initializer<T, J, DoF>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this, shift = process_rows.front()](const size_t row) { 
        matrix()[matrix_part::INNER].innerIndexPtr()[matrix()[matrix_part::INNER].outerIndexPtr()[row]] = row + shift; 
    });
    if (sort_indices) {
        utils::sort_indices(matrix()[matrix_part::INNER]);
        utils::sort_indices(matrix()[matrix_part::BOUND]);
    }
}

template<class T, class I, class J, size_t DoF>
template<class Local_Integrator, class Nonlocal_Integrator>
void matrix_assembler_2d<T, I, J, DoF>::calc_coeffs(
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric,
    Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, integrator<T, J, DoF, Local_Integrator, Nonlocal_Integrator>{
        _matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric, local_integrator, nonlocal_integrator});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { 
        matrix()[matrix_part::INNER].valuePtr()[matrix()[matrix_part::INNER].outerIndexPtr()[row]] = T{1};
    });
}

}

#endif