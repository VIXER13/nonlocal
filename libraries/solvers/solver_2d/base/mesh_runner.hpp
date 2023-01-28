#ifndef NONLOCAL_MESH_RUNNER_HPP
#define NONLOCAL_MESH_RUNNER_HPP

#include <array>
#include <ranges>
#include <vector>

namespace nonlocal {

enum matrix_part : size_t { INNER, BOUND };

template<class T, class I>
using matrix_parts_t = std::array<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, 2>;

class mesh_runner_base {
    const std::vector<bool>& _is_inner;

protected:
    explicit mesh_runner_base(const std::vector<bool>& is_inner);

    template<class Callback>
    void filter(const size_t row, const size_t col, const Callback& callback);

public:
    virtual ~mesh_runner_base() noexcept = default;
};

mesh_runner_base::mesh_runner_base(const std::vector<bool>& is_inner)
    : _is_inner{is_inner} {}

template<class Callback>
void mesh_runner_base::filter(const size_t row, const size_t col, const Callback& callback) {
    if (_is_inner[col]) {
        if (row <= col && _is_inner[row])
            callback(INNER, row, col);
    } else if (row != col)
        callback(BOUND, row, col);
}

template<size_t DoF, class T, class I>
class indexator_base : public mesh_runner_base {
    using flags_t = std::array<std::array<std::vector<bool>, DoF>, 2>;
    flags_t _flags;
    matrix_parts_t<T, I>& _matrix;
    const size_t _node_shift;

protected:
    template<class Callback>
    static void check_flag(std::vector<bool>& flags, const size_t col, const Callback& callback);

    explicit indexator_base(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift);

    size_t node_shift() const noexcept;
    std::array<std::vector<bool>, DoF>& flags(const matrix_part part) noexcept;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix(const matrix_part part) noexcept;

public:
    ~indexator_base() noexcept override = default;

    void reset(const size_t node);
};

template<size_t DoF, class T, class I>
template<class Callback>
void indexator_base<DoF, T, I>::check_flag(std::vector<bool>& flags, const size_t col, const Callback& callback) {
    if (!flags[col]) {
        callback();
        flags[col] = true;
    }
}

template<size_t DoF, class T, class I>
indexator_base<DoF, T, I>::indexator_base(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift)
    : mesh_runner_base{is_inner}, _matrix{matrix}, _node_shift{node_shift} {
    for(auto& part : _flags)
        for(auto& dof : part)
            dof.resize(is_inner.size());
}

template<size_t DoF, class T, class I>
size_t indexator_base<DoF, T, I>::node_shift() const noexcept {
    return _node_shift;
}

template<size_t DoF, class T, class I>
std::array<std::vector<bool>, DoF>& indexator_base<DoF, T, I>::flags(const matrix_part part) noexcept {
    return _flags[part];
}

template<size_t DoF, class T, class I>
Eigen::SparseMatrix<T, Eigen::RowMajor, I>& indexator_base<DoF, T, I>::matrix(const matrix_part part) noexcept {
    return _matrix[part];
}

template<size_t DoF, class T, class I>
void indexator_base<DoF, T, I>::reset(const size_t node) {
    for(const size_t i : std::ranges::iota_view{0u, DoF}) {
        std::fill(_flags[BOUND][i].begin(), _flags[BOUND][i].end(), false);
        std::fill(std::next(_flags[INNER][i].begin(), DoF * node), _flags[INNER][i].end(), false);
    }
}

template<size_t DoF, class T, class I>
class shift_initializer final : public indexator_base<DoF, T, I> {
    using _base = indexator_base<DoF, T, I>;

public:
    explicit shift_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift);
    ~shift_initializer() noexcept override = default;

    void operator()(const size_t row_block, const size_t col_block);
};

template<size_t DoF, class T, class I>
shift_initializer<DoF, T, I>::shift_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift)
    : _base{matrix, is_inner, node_shift} {}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::operator()(const size_t row_block, const size_t col_block) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF})
            _base::filter(DoF * row_block + row_loc, DoF * col_block + col_loc, 
                [this, row_loc](const matrix_part part, const size_t row, const size_t col) {
                    _base::check_flag(_base::flags(part)[row_loc], col, [this, part, row]() {
                        ++_base::matrix(part).outerIndexPtr()[row - DoF * _base::node_shift() + 1];
                    });
            });
}

template<size_t DoF, class T, class I>
class index_initializer final : public indexator_base<DoF, T, I> {
    using _base = indexator_base<DoF, T, I>;
    using indices_t = std::array<std::array<size_t, DoF>, 2>;
    indices_t _current_indices = {};

public:
    explicit index_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift);
    ~index_initializer() noexcept override = default;

    void reset(const size_t node);

    void operator()(const size_t row_block, const size_t col_block);
};

template<size_t DoF, class T, class I>
index_initializer<DoF, T, I>::index_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift)
    : _base{matrix, is_inner, node_shift} {}

template<size_t DoF, class T, class I>
void index_initializer<DoF, T, I>::reset(const size_t node) {
    for(const matrix_part part : {matrix_part::INNER, matrix_part::BOUND})
        for(const size_t dof : std::ranges::iota_view{0u, DoF})
            _current_indices[part][dof] = _base::matrix(part).outerIndexPtr()[DoF * (node - _base::node_shift()) + dof];
    _base::reset(node);
}

template<size_t DoF, class T, class I>
void index_initializer<DoF, T, I>::operator()(const size_t row_block, const size_t col_block) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF})
            _base::filter(DoF * row_block + row_loc, DoF * col_block + col_loc, 
                [this, row_loc](const matrix_part part, const size_t row, const size_t col) {
                    _base::check_flag(_base::flags(part)[row_loc], col, [this, part, col, row_loc]() {
                        _base::matrix(part).innerIndexPtr()[_current_indices[part][row_loc]++] = col;
                    });
            });
}

template<class Callback>
void first_kind_indexator(const std::ranges::iota_view<size_t, size_t> rows, const std::vector<bool>& is_inner, const Callback& callback) {
    for(const size_t row : rows)
        if (!is_inner[row])
            callback(row);
}

}

#endif