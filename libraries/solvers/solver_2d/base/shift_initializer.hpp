#pragma once

#include "mesh_container_2d.hpp"

#include "matrix_indexator_base.hpp"
#include "matrix_separator_base.hpp"

namespace nonlocal {

template<class T, class I, class J, size_t DoF>
class shift_initializer final : public matrix_separator_base<T, J>
                              , public matrix_indexator_base<DoF> {
    using _matrix = matrix_separator_base<T, J>;
    using _indexator = matrix_indexator_base<DoF>;

    const mesh::mesh_container_2d<T, I>& _mesh;

    void run(const size_t row_glob, const size_t col_glob);

public:
    explicit shift_initializer(finite_element_matrix<T, J>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                               const std::vector<bool>& is_inner, const size_t node_shift, const bool is_symmetric);
    ~shift_initializer() noexcept override = default;

    void operator()(const std::string&, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<class T, class I, class J, size_t DoF>
shift_initializer<T, I, J, DoF>::shift_initializer(finite_element_matrix<T, J>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                                                   const std::vector<bool>& is_inner, const size_t node_shift, const bool is_symmetric)
    : _matrix{matrix, is_inner, node_shift, is_symmetric}
    , _indexator{is_inner.size(), is_symmetric}
    , _mesh{mesh} {}

template<class T, class I, class J, size_t DoF>
void shift_initializer<T, I, J, DoF>::run(const size_t row_glob, const size_t col_glob) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            const matrix_part part = _matrix::part(row, col);
            _indexator::check_flag(part, row_loc, col, [this, part, row]() {
                ++_matrix::matrix(part).outerIndexPtr()[row - DoF * _matrix::node_shift() + 1];
            });
        }
}

template<class T, class I, class J, size_t DoF>
void shift_initializer<T, I, J, DoF>::operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
    run(DoF * _mesh.node_number(e, i), DoF * _mesh.node_number(e, j));
}

template<class T, class I, class J, size_t DoF>
void shift_initializer<T, I, J, DoF>::operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    run(DoF * _mesh.node_number(eL, iL), DoF * _mesh.node_number(eNL, jNL));
}

}