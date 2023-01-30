#ifndef NONLOCAL_FINITE_ELEMENT_MATRIX_2D_HPP
#define NONLOCAL_FINITE_ELEMENT_MATRIX_2D_HPP

#include "../solvers_utils.hpp"

#include "shift_initializer.hpp"
#include "index_initializer.hpp"
#include "integrator.hpp"

#include "mesh_2d.hpp"

#include <iostream>

namespace nonlocal {

template<size_t DoF, class T, class I, class Matrix_Index>
class finite_element_matrix_2d {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    matrix_parts_t<T, Matrix_Index> _matrix;

protected:
    template<class Callback>
    static void first_kind_filler(const std::ranges::iota_view<size_t, size_t> rows, 
                                  const std::vector<bool>& is_inner, const Callback& callback);

    explicit finite_element_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);

    template<class Initializer>
    void mesh_run(const std::unordered_map<std::string, theory_t>& theories, Initializer&& initializer);

    void init_shifts(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner);
    void init_indices(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner, const bool sort_indices = true);
    template<class Integrate_Loc, class Integrate_Nonloc>
    void calc_coeffs(const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner,
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
template<class Callback>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::first_kind_filler(const std::ranges::iota_view<size_t, size_t> rows,
                                                                          const std::vector<bool>& is_inner,
                                                                          const Callback& callback) {
    for(const size_t row : rows)
        if (!is_inner[row])
            callback(row);
}

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
#pragma omp parallel for default(none) shared(theories, process_nodes) firstprivate(initializer)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        if constexpr (std::is_base_of_v<indexator_base<DoF>, Initializer>)
            initializer.reset(node);
        for(const I eL : mesh().elements(node)) {
            const std::string& group = mesh().container().group(eL);
            if (const theory_t theory = theories.at(group); theory == theory_t::LOCAL)
                for(const size_t jL : std::ranges::iota_view{0u, mesh().container().nodes_count(eL)})
                    initializer(group, eL, node, mesh().container().node_number(eL, jL));
            else if (theory == theory_t::NONLOCAL)
                for(const I eNL : mesh().neighbours(eL))
                    for(const size_t jNL : std::ranges::iota_view{0u, mesh().container().nodes_count(eNL)})
                        initializer(group, eL, eNL, node, mesh().container().node_number(eNL, jNL));
            else
                throw std::domain_error{"Unknown theory."};
        }
    }
}

template<size_t DoF, class T, class I, class Matrix_Index>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::init_shifts(
    const std::unordered_map<std::string, theory_t>& theories, 
    const std::vector<bool>& is_inner) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, shift_initializer<DoF, T, Matrix_Index>{_matrix, is_inner, process_nodes.front()});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { ++matrix_inner().outerIndexPtr()[row + 1]; });
    utils::accumulate_shifts(matrix_inner());
    utils::accumulate_shifts(matrix_bound());
    std::cout << "Non-zero elements count: " << matrix_inner().nonZeros() + matrix_bound().nonZeros() << std::endl;
}

template<size_t DoF, class T, class I, class Matrix_Index>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::init_indices(
    const std::unordered_map<std::string, theory_t>& theories, 
    const std::vector<bool>& is_inner, const bool sort_indices) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    utils::allocate_matrix(matrix_inner());
    utils::allocate_matrix(matrix_bound());
    mesh_run(theories, index_initializer<DoF, T, Matrix_Index>{_matrix, is_inner, process_nodes.front()});
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
    const std::unordered_map<std::string, theory_t>& theories, const std::vector<bool>& is_inner,
    Integrate_Loc&& integrate_loc, Integrate_Nonloc&& integrate_nonloc) {
    const auto process_nodes = mesh().process_nodes();
    const auto process_rows = std::ranges::iota_view{DoF * process_nodes.front(), DoF * *process_nodes.end()};
    mesh_run(theories, integrator<DoF, T, Matrix_Index, Integrate_Loc, Integrate_Nonloc>{
        _matrix, is_inner, process_nodes.front(), integrate_loc, integrate_nonloc});
    first_kind_filler(process_rows, is_inner, [this](const size_t row) { 
        matrix_inner().valuePtr()[matrix_inner().outerIndexPtr()[row]] = T{1}; 
    });
}

/*
        const auto assemble_predicate =
        [&is_inner](const size_t glob_row, const size_t glob_col, const theory_t theory) {
            for(const size_t row : std::views::iota(glob_row, glob_row + DoF))
                for(const size_t col : std::views::iota(glob_col, glob_col + DoF))
                    if (is_inner[row] && is_inner[col] ? row <= col     :
                        row != col                     ? !is_inner[col] : theory == theory_t::LOCAL)
                        return true;
            return false;
        };

    using block_t = std::array<T, DoF * DoF>;
    const auto assemble_block =
        [this, &is_inner, shift = DoF * mesh().process_nodes().front()]
        (const block_t& block, const size_t glob_row, const size_t glob_col, const theory_t theory) {
            for(auto [row, component] = std::pair{glob_row, block.cbegin()}; row < glob_row + DoF; ++row) {
                T* data_ptr = nullptr;
                for(size_t col = glob_col; col < glob_col + DoF; ++col, ++component)
                    if (is_inner[row] && is_inner[col]) {
                        if (row <= col) {
                            if (!data_ptr) data_ptr = &matrix_inner().coeffRef(row - shift, col);
                            else         ++data_ptr;
                            *data_ptr += *component;
                        }
                    } else if (row != col) {
                        if (!is_inner[col])
                            matrix_bound().coeffRef(row - shift, col) += *component;
                    } else if (theory == theory_t::LOCAL)
                        matrix_inner().coeffRef(row - shift, col) = T{1};
            }
        };

    mesh_run<theory_t::LOCAL>(
        [this, &integrate_rule_loc, &assemble_predicate, &assemble_block](const size_t e, const size_t i, const size_t j) {
            const size_t glob_row = DoF * mesh().container().node_number(e, i),
                         glob_col = DoF * mesh().container().node_number(e, j);
            if (assemble_predicate(glob_row, glob_col, theory_t::LOCAL))
                assemble_block(block_t{integrate_rule_loc(e, i, j)}, glob_row, glob_col, theory_t::LOCAL);
        });

    if (theory == theory_t::NONLOCAL)
        mesh_run<theory_t::NONLOCAL>(
            [this, &influence_fun, &integrate_rule_nonloc, &assemble_predicate, &assemble_block]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const size_t glob_row = DoF * mesh().container().node_number(eL,   iL),
                             glob_col = DoF * mesh().container().node_number(eNL, jNL);
                if (assemble_predicate(glob_row, glob_col, theory_t::NONLOCAL))
                    assemble_block(block_t{integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun)}, glob_row, glob_col, theory_t::NONLOCAL);
            });
*/

/*
template<size_t DoF, class T, class I, class Matrix_Index>
class finite_element_matrix_2d {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    std::shared_ptr<mesh::mesh_2d<T, I>> _mesh;
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index> _matrix_inner;
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index> _matrix_bound;

protected:
    enum class index_stage : bool { SHIFTS, NONZERO };

    // TODO: indexator refactoring
    class indexator {
        const std::vector<bool>& _inner_nodes;
        const size_t _node_shift;
        Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& _K_inner,
                                                             & _K_bound;
        std::array<std::vector<bool>, DoF> _inner, _bound;
        std::array<size_t, DoF> _inner_index = {}, _bound_index = {};

    public:
        explicit indexator(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                           Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                           const std::vector<bool>& inner_nodes, const size_t node_shift)
            : _inner_nodes{inner_nodes}
            , _node_shift{node_shift}
            , _K_inner{K_inner}
            , _K_bound{K_bound} {
            for(size_t i = 0; i < DoF; ++i) {
                _inner[i].resize(_inner_nodes.size());
                _bound[i].resize(_inner_nodes.size());
            }
        }

        void fill(const size_t node) {
            for(size_t i = 0; i < DoF; ++i) {
                _inner_index[i] = _K_inner.outerIndexPtr()[DoF * (node - _node_shift) + i];
                _bound_index[i] = _K_bound.outerIndexPtr()[DoF * (node - _node_shift) + i];
                std::fill(_bound[i].begin(), _bound[i].end(), false);
                std::fill(std::next(_inner[i].begin(), DoF * node), _inner[i].end(), false);
            }
        }

        template<index_stage Stage>
        void stage(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                   std::array<std::vector<bool>, DoF>& inner,
                   std::array<size_t, DoF>& inner_index,
                   const size_t row, const size_t col) {
            const size_t i = row % DoF;
            if (!inner[i][col]) {
                if constexpr (Stage == index_stage::SHIFTS)
                    ++K.outerIndexPtr()[row - DoF * _node_shift + 1];
                if constexpr (Stage == index_stage::NONZERO)
                    K.innerIndexPtr()[inner_index[i]++] = col;
                inner[i][col] = true;
            }
        }

        template<index_stage Stage>
        void index(const size_t block_row, const size_t block_col) {
            for(size_t comp_row = 0; comp_row < DoF; ++comp_row)
                for(size_t comp_col = 0; comp_col < DoF; ++comp_col) {
                    const size_t row = DoF * block_row + comp_row,
                                 col = DoF * block_col + comp_col;
                    if (_inner_nodes[row] && _inner_nodes[col]) {
                        if (row <= col)
                            stage<Stage>(_K_inner, _inner, _inner_index, row, col);
                    } else if (row != col) {
                        if (!_inner_nodes[col])
                            stage<Stage>(_K_bound, _bound, _bound_index, row, col);
                    } else
                        stage<Stage>(_K_inner, _inner, _inner_index, row, col);
                }
        }
    };

    template<index_stage Stage, theory_t Theory>
    void mesh_index(const std::vector<bool>& inner_nodes) {
        const auto process_nodes = mesh().process_nodes();
        indexator ind{matrix_inner(), matrix_bound(), inner_nodes, process_nodes.front()};
#pragma omp parallel for default(none) shared(process_nodes) firstprivate(ind) schedule(dynamic)
        for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
            ind.fill(node);
            for(const I eL : mesh().elements(node)) {
                if constexpr (Theory == theory_t::LOCAL)
                    for(const size_t jL : std::ranges::iota_view{0u, mesh().container().nodes_count(eL)})
                        ind.template index<Stage>(node, mesh().container().node_number(eL, jL));
                if constexpr (Theory == theory_t::NONLOCAL)
                    for(const I eNL : mesh().neighbours(eL))
                        for(const size_t jNL : std::ranges::iota_view{0u, mesh().container().nodes_count(eNL)})
                            ind.template index<Stage>(node, mesh().container().node_number(eNL, jNL));
            }
        }
    }
    
    explicit finite_element_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);

    template<theory_t Theory>
    void create_matrix_portrait(const std::vector<bool>& is_inner);

    template<theory_t Theory, class Callback>
    void mesh_run(const Callback& callback) const;

    template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
    void calc_matrix(const std::vector<bool>& is_inner,
                     const theory_t theory,
                     const Influence_Function& influence_fun,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc);

public:
    virtual ~finite_element_matrix_2d() noexcept = default;

    const mesh::mesh_2d<T, I>& mesh() const;
    const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh_ptr() const;
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_inner();
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_bound();
    const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_inner() const;
    const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& matrix_bound() const;

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
const std::shared_ptr<mesh::mesh_2d<T, I>>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::mesh_ptr() const {
    return _mesh;
}

template<size_t DoF, class T, class I, class Matrix_Index>
Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_inner() {
    return _matrix_inner;
}

template<size_t DoF, class T, class I, class Matrix_Index>
Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_bound() {
    return _matrix_bound;
}

template<size_t DoF, class T, class I, class Matrix_Index>
const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_inner() const {
    return _matrix_inner;
}

template<size_t DoF, class T, class I, class Matrix_Index>
const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& finite_element_matrix_2d<DoF, T, I, Matrix_Index>::matrix_bound() const {
    return _matrix_bound;
}

template<size_t DoF, class T, class I, class Matrix_Index>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::clear() {
    _matrix_inner = Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>{};
    _matrix_bound = Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>{};
}

template<size_t DoF, class T, class I, class Matrix_Index>
template<theory_t Theory>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::create_matrix_portrait(const std::vector<bool>& is_inner) {
    mesh_index<index_stage::SHIFTS, Theory>(is_inner);
    utils::accumulate_shifts(matrix_inner());
    utils::accumulate_shifts(matrix_bound());
    std::cout << "Non-zero elements count: " << matrix_inner().nonZeros() + matrix_bound().nonZeros() << std::endl;
    utils::allocate_matrix(matrix_inner());
    utils::allocate_matrix(matrix_bound());
    mesh_index<index_stage::NONZERO, Theory>(is_inner);
}

template<size_t DoF, class T, class I, class Matrix_Index>
template<theory_t Theory, class Callback>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::mesh_run(const Callback& callback) const {
    const auto process_nodes = mesh().process_nodes();
#pragma omp parallel for default(none) shared(process_nodes) firstprivate(callback) schedule(dynamic)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        for(const I eL : mesh().elements(node)) {
            const size_t iL = mesh().global_to_local(eL, node);
            if constexpr (Theory == theory_t::LOCAL)
                for(const size_t jL : std::ranges::iota_view{0u, mesh().container().nodes_count(eL)})
                    callback(eL, iL, jL);
            if constexpr (Theory == theory_t::NONLOCAL)
                for(const I eNL : mesh().neighbours(eL))
                    for(const size_t jNL : std::ranges::iota_view{0u, mesh().container().nodes_count(eNL)})
                        callback(eL, eNL, iL, jNL);
        }
    }
}

template<size_t DoF, class T, class I, class Matrix_Index>
template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
void finite_element_matrix_2d<DoF, T, I, Matrix_Index>::calc_matrix(const std::vector<bool>& is_inner,
                                                                    const theory_t theory,
                                                                    const Influence_Function& influence_fun,
                                                                    const Integrate_Loc& integrate_rule_loc,
                                                                    const Integrate_Nonloc& integrate_rule_nonloc) {
    const auto assemble_predicate =
        [&is_inner](const size_t glob_row, const size_t glob_col, const theory_t theory) {
            for(const size_t row : std::views::iota(glob_row, glob_row + DoF))
                for(const size_t col : std::views::iota(glob_col, glob_col + DoF))
                    if (is_inner[row] && is_inner[col] ? row <= col     :
                        row != col                     ? !is_inner[col] : theory == theory_t::LOCAL)
                        return true;
            return false;
        };

    using block_t = std::array<T, DoF * DoF>;
    const auto assemble_block =
        [this, &is_inner, shift = DoF * mesh().process_nodes().front()]
        (const block_t& block, const size_t glob_row, const size_t glob_col, const theory_t theory) {
            for(auto [row, component] = std::pair{glob_row, block.cbegin()}; row < glob_row + DoF; ++row) {
                T* data_ptr = nullptr;
                for(size_t col = glob_col; col < glob_col + DoF; ++col, ++component)
                    if (is_inner[row] && is_inner[col]) {
                        if (row <= col) {
                            if (!data_ptr) data_ptr = &matrix_inner().coeffRef(row - shift, col);
                            else         ++data_ptr;
                            *data_ptr += *component;
                        }
                    } else if (row != col) {
                        if (!is_inner[col])
                            matrix_bound().coeffRef(row - shift, col) += *component;
                    } else if (theory == theory_t::LOCAL)
                        matrix_inner().coeffRef(row - shift, col) = T{1};
            }
        };

    mesh_run<theory_t::LOCAL>(
        [this, &integrate_rule_loc, &assemble_predicate, &assemble_block](const size_t e, const size_t i, const size_t j) {
            const size_t glob_row = DoF * mesh().container().node_number(e, i),
                         glob_col = DoF * mesh().container().node_number(e, j);
            if (assemble_predicate(glob_row, glob_col, theory_t::LOCAL))
                assemble_block(block_t{integrate_rule_loc(e, i, j)}, glob_row, glob_col, theory_t::LOCAL);
        });

    if (theory == theory_t::NONLOCAL)
        mesh_run<theory_t::NONLOCAL>(
            [this, &influence_fun, &integrate_rule_nonloc, &assemble_predicate, &assemble_block]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const size_t glob_row = DoF * mesh().container().node_number(eL,   iL),
                             glob_col = DoF * mesh().container().node_number(eNL, jNL);
                if (assemble_predicate(glob_row, glob_col, theory_t::NONLOCAL))
                    assemble_block(block_t{integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun)}, glob_row, glob_col, theory_t::NONLOCAL);
            });
}
*/

}

#endif
