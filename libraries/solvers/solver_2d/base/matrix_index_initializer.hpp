#pragma once

#include "matrix_indexator_base.hpp"
#include "matrix_separator_base.hpp"

namespace nonlocal::solver_2d {

template<class T, class I, class J, size_t DoF>
class matrix_index_initializer final : public matrix_separator_base<T, J>
                                     , public matrix_indexator_base<DoF> {
    using _matrix = matrix_separator_base<T, J>;
    using _indexator = matrix_indexator_base<DoF>;

    const mesh::mesh_container_2d<T, I>& _mesh;
    std::array<std::array<size_t, DoF>, 2> _current_indices = {};

    void run(const size_t row_glob, const size_t col_glob);

public:
    explicit matrix_index_initializer(finite_element_matrix<T, J>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                                      const std::vector<bool>& is_inner, const size_t node_shift, const bool is_symmetric);
    ~matrix_index_initializer() noexcept override = default;

    void reset(const size_t node);

    void operator()(const std::string&, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<class T, class I, class J, size_t DoF>
matrix_index_initializer<T, I, J, DoF>::matrix_index_initializer(finite_element_matrix<T, J>& matrix, const mesh::mesh_container_2d<T, I>& mesh, 
                                                                 const std::vector<bool>& is_inner, const size_t node_shift, const bool is_symmetric)
    : _matrix{matrix, is_inner, node_shift, is_symmetric}
    , _indexator{is_inner.size(), is_symmetric}
    , _mesh{mesh} {}

template<class T, class I, class J, size_t DoF>
void matrix_index_initializer<T, I, J, DoF>::reset(const size_t node) {
    for(const matrix_part part : {matrix_part::INNER, matrix_part::BOUND})
        for(const size_t dof : std::ranges::iota_view{0u, DoF})
            _current_indices[size_t(part)][dof] = _matrix::matrix(part).outerIndexPtr()[DoF * (node - _matrix::node_shift()) + dof];
    _indexator::reset(node);
}

template<class T, class I, class J, size_t DoF>
void matrix_index_initializer<T, I, J, DoF>::run(const size_t row_glob, const size_t col_glob) {
    for(const size_t row_loc : std::ranges::iota_view{0u, DoF})
        for(const size_t col_loc : std::ranges::iota_view{0u, DoF}) {
            const size_t row = row_glob + row_loc;
            const size_t col = col_glob + col_loc;
            const matrix_part part = _matrix::part(row, col);
            _indexator::check_flag(part, row_loc, col, [this, part, col, row_loc]() {
                _matrix::matrix(part).innerIndexPtr()[_current_indices[size_t(part)][row_loc]++] = col;
            });
        }
}

template<class T, class I, class J, size_t DoF>
void matrix_index_initializer<T, I, J, DoF>::operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
    run(DoF * _mesh.node_number(e, i), DoF * _mesh.node_number(e, j));
}

template<class T, class I, class J, size_t DoF>
void matrix_index_initializer<T, I, J, DoF>::operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    run(DoF * _mesh.node_number(eL, iL), DoF * _mesh.node_number(eNL, jNL));
}

}