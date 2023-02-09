#ifndef NONLOCAL_SHIFT_INITIALIZER_HPP
#define NONLOCAL_SHIFT_INITIALIZER_HPP

#include "mesh_container_2d.hpp"

#include "indexator_base.hpp"
#include "matrix_separator_base.hpp"

namespace nonlocal {

template<size_t DoF, class T, class I>
class shift_initializer final : public matrix_separator_base<T, I>
                              , public indexator_base<DoF> {
    using _matrix = matrix_separator_base<T, I>;
    using _indexator = indexator_base<DoF>;

    const mesh::mesh_container_2d<T, I>& _mesh;

    void run(const size_t row_glob, const size_t col_glob);

public:
    explicit shift_initializer(matrix_parts_t<T, I>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                               const std::vector<bool>& is_inner, const size_t node_shift);
    ~shift_initializer() noexcept override = default;

    void operator()(const std::string&, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<size_t DoF, class T, class I>
shift_initializer<DoF, T, I>::shift_initializer(matrix_parts_t<T, I>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                                                const std::vector<bool>& is_inner, const size_t node_shift)
    : _matrix{matrix, is_inner, node_shift}
    , _indexator{is_inner.size()}
    , _mesh{mesh} {}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::run(const size_t row_glob, const size_t col_glob) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            if (const matrix_part part = _matrix::part(row, col); part != matrix_part::NO)
                _indexator::check_flag(_indexator::flags(part)[row_loc], col, [this, part, row]() {
                    ++_matrix::matrix(part).outerIndexPtr()[row - DoF * _matrix::node_shift() + 1];
                });
        }
}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
    run(DoF * _mesh.node_number(e, i), DoF * _mesh.node_number(e, j));
}

template<size_t DoF, class T, class I>
void shift_initializer<DoF, T, I>::operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    run(DoF * _mesh.node_number(eL, iL), DoF * _mesh.node_number(eNL, jNL));
}

}

#endif