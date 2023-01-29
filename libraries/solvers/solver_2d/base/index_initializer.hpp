#ifndef NONLOCAL_INDEX_INITIALIZER_HPP
#define NONLOCAL_INDEX_INITIALIZER_HPP

#include "indexator_base.hpp"
#include "matrix_separator_base.hpp"

namespace nonlocal {

template<size_t DoF, class T, class I>
class index_initializer final : public matrix_separator_base<T, I>
                              , public indexator_base<DoF> {
    using _matrix = matrix_separator_base<T, I>;
    using _indexator = indexator_base<DoF>;
    std::array<std::array<size_t, DoF>, 2> _current_indices = {};

    void run(const size_t row_block, const size_t col_block);

public:
    explicit index_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift);
    ~index_initializer() noexcept override = default;

    void reset(const size_t node);

    void operator()(const std::string&, const size_t, const size_t row_block, const size_t col_block);
    void operator()(const std::string&, const size_t, const size_t, const size_t row_block, const size_t col_block);
};

template<size_t DoF, class T, class I>
index_initializer<DoF, T, I>::index_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift)
    : _matrix{matrix, is_inner, node_shift}
    , _indexator{is_inner.size()} {}

template<size_t DoF, class T, class I>
void index_initializer<DoF, T, I>::reset(const size_t node) {
    for(const matrix_part part : {matrix_part::INNER, matrix_part::BOUND})
        for(const size_t dof : std::ranges::iota_view{0u, DoF})
            _current_indices[size_t(part)][dof] = _matrix::matrix(part).outerIndexPtr()[DoF * (node - _matrix::node_shift()) + dof];
    _indexator::reset(node);
}

template<size_t DoF, class T, class I>
void index_initializer<DoF, T, I>::run(const size_t row_block, const size_t col_block) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = DoF * row_block + row_loc;
            const size_t col = DoF * col_block + col_loc;
            if (const matrix_part part = _matrix::part(row, col); part != matrix_part::NO)
                _indexator::check_flag(_indexator::flags(part)[row_loc], col, [this, part, col, row_loc]() {
                    _matrix::matrix(part).innerIndexPtr()[_current_indices[size_t(part)][row_loc]++] = col;
                });
        }
}

template<size_t DoF, class T, class I>
void index_initializer<DoF, T, I>::operator()(const std::string&, const size_t, const size_t row_block, const size_t col_block) {
    run(row_block, col_block);
}

template<size_t DoF, class T, class I>
void index_initializer<DoF, T, I>::operator()(const std::string&, const size_t, const size_t, const size_t row_block, const size_t col_block) {
    run(row_block, col_block);
}

}

#endif