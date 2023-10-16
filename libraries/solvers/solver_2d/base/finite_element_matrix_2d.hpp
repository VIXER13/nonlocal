#ifndef NONLOCAL_FINITE_ELEMENT_MATRIX_2D_HPP
#define NONLOCAL_FINITE_ELEMENT_MATRIX_2D_HPP

#include "../solvers_utils.hpp"

#include "shift_initializer.hpp"
#include "index_initializer.hpp"
#include "integrator.hpp"

#include "mesh_2d.hpp"

#include <iostream>

namespace nonlocal {

template<class Callback>
void first_kind_filler(const std::ranges::iota_view<size_t, size_t> rows, 
                       const std::vector<bool>& is_inner, const Callback& callback) {
    for(const size_t row : rows)
        if (!is_inner[row])
            callback(row);
}

template<size_t DoF, class T, class I, class Matrix_Index>
class finite_element_matrix_2d {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    matrix_parts_t<T, Matrix_Index> _matrix;

protected:
    explicit finite_element_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);

    template<class Initializer>
    void mesh_run(const std::unordered_map<std::string, theory_t>& theories, Initializer&& initializer);

    void init_shifts(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric);
    void init_indices(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner,
                      const bool is_symmetric, const bool sort_indices = true);
    template<class Integrate_Loc, class Integrate_Nonloc>
    void calc_coeffs(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric,
                     Integrate_Loc&& integrate_loc, Integrate_Nonloc&& integrate_nonloc);

public:
    virtual ~finite_element_matrix_2d() noexcept = default;

    const mesh::mesh_2d<T, I>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const noexcept;
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_inner() noexcept;
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_bound() noexcept;
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix(const matrix_part part) noexcept;
    const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_inner() const noexcept;
    const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_bound() const noexcept;
    const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix(const matrix_part part) const noexcept;

    void clear();
};

template<size_t DoF, class T, class I, class Matrix_Index>
finite_element_matrix_2d<DoF, T, I, Matrix_Index>::finite_element_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _mesh{mesh} {}

template<size_t DoF, class T, class I, class Matrix_Index>
const mesh::mesh_2d<T, I>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::mesh() const {
    return *mesh_ptr();
}

template<size_t DoF, class T, class I, class Matrix_Index>
const std::shared_ptr<mesh::mesh_2d<T, I>>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::mesh_ptr() const noexcept {
    return _mesh;
}

template<size_t DoF, class T, class I, class Matrix_Index>
Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_inner() noexcept {
    return matrix(matrix_part::INNER);
}

template<size_t DoF, class T, class I, class Matrix_Index>
Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_bound() noexcept {
    return matrix(matrix_part::BOUND);
}

template<size_t DoF, class T, class I, class Matrix_Index>
Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix(const matrix_part part) noexcept {
    return _matrix[size_t(part)];
}

template<size_t DoF, class T, class I, class Matrix_Index>
const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_inner() const noexcept {
    return matrix(matrix_part::INNER);
}

template<size_t DoF, class T, class I, class Matrix_Index>
const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_bound() const noexcept {
    return matrix(matrix_part::BOUND);
}

template<size_t DoF, class T, class I, class Matrix_Index>
const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix(const matrix_part part) const noexcept {
    return _matrix[size_t(part)];
}

template<size_t DoF, class T, class I, class Matrix_Index>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::clear() {
    matrix_inner() = Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>{};
    matrix_bound() = Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>{};
}

template<size_t DoF, class T, class I, class Matrix_Index>
template<class Initializer>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::mesh_run(const std::unordered_map<std::string, theory_t>& theories,
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

template<size_t DoF, class T, class I, class Matrix_Index>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::init_shifts(
    const std::unordered_map<std::string, theory_t>& theories, 
    const std::vector<bool>& is_inner, const bool is_symmetric) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, shift_initializer<DoF, T, Matrix_Index>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { ++matrix_inner().outerIndexPtr()[row + 1]; });
    utils::accumulate_shifts(matrix_inner());
    utils::accumulate_shifts(matrix_bound());
    std::cout << "Non-zero elements count: " << matrix_inner().nonZeros() + matrix_bound().nonZeros() << std::endl;
}

template<size_t DoF, class T, class I, class Matrix_Index>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::init_indices(
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, 
    const bool is_symmetric, const bool sort_indices) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    utils::allocate_matrix(matrix_inner());
    utils::allocate_matrix(matrix_bound());
    mesh_run(theories, index_initializer<DoF, T, Matrix_Index>{_matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric});
    first_kind_filler(process_rows, is_inner, [this, shift = process_rows.front()](const size_t row) { 
        matrix_inner().innerIndexPtr()[matrix_inner().outerIndexPtr()[row]] = row + shift; 
    });
    if (sort_indices) {
        utils::sort_indices(matrix_inner());
        utils::sort_indices(matrix_bound());
    }
}

template<size_t DoF, class T, class I, class Matrix_Index>
template<class Integrate_Loc, class Integrate_Nonloc>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::calc_coeffs(
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool is_symmetric,
    Integrate_Loc&& integrate_loc, Integrate_Nonloc&& integrate_nonloc) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, integrator<DoF, T, Matrix_Index, Integrate_Loc, Integrate_Nonloc>{
        _matrix, mesh().container(), is_inner, process_nodes.front(), is_symmetric, integrate_loc, integrate_nonloc});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { 
        matrix_inner().valuePtr()[matrix_inner().outerIndexPtr()[row]] = T{1};
    });
}

}

#endif