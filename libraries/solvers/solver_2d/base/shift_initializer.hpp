#ifndef NONLOCAL_SHIFT_INITIALIZER_HPP
#define NONLOCAL_SHIFT_INITIALIZER_HPP

#include "indexator_base.hpp"
#include "matrix_separator_base.hpp"

namespace nonlocal {

template<size_t DoF, class T, class I>
class shift_initializer final : public matrix_separator_base<T, I>
                              , public indexator_base<DoF> {
    using _matrix = matrix_separator_base<T, I>;
    using _indexator = indexator_base<DoF>;

    void run(const size_t row_block, const size_t col_block);

public:
    explicit shift_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift);
    ~shift_initializer() noexcept override = default;

    void operator()(const std::string&, const size_t, const size_t row_block, const size_t col_block);
    void operator()(const std::string&, const size_t, const size_t, const size_t row_block, const size_t col_block);
};

template<size_t DoF, class T, class I>
shift_initializer<DoF, T, I>::shift_initializer(matrix_parts_t<T, I>& matrix, const std::vector<bool>& is_inner, const size_t node_shift)
    : _matrix{matrix, is_inner, node_shift}
    , _indexator{is_inner.size()} {}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::run(const size_t row_block, const size_t col_block) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = DoF * row_block + row_loc;
            const size_t col = DoF * col_block + col_loc;
            if (const matrix_part part = _matrix::part(row, col); part != matrix_part::NO)
                _indexator::check_flag(_indexator::flags(part)[row_loc], col, [this, part, row]() {
                    ++_matrix::matrix(part).outerIndexPtr()[row - DoF * _matrix::node_shift() + 1];
                });
        }
}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::operator()(const std::string&, const size_t, const size_t row_block, const size_t col_block) {
    run(row_block, col_block);
}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::operator()(const std::string&, const size_t, const size_t, const size_t row_block, const size_t col_block) {
    run(row_block, col_block);
}

}

#endif