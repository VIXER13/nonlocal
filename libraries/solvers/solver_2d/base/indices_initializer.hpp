#pragma once

#include <mesh/mesh_2d/indexator_base.hpp>
#include <mesh/mesh_2d/mesh_container_2d.hpp>

namespace nonlocal::solver_2d {

template<std::floating_point T>
class indices_initializer final : public mesh::indexator_base {
    using _base = mesh::indexator_base;

    std::vector<bool> _included;
    size_t _current_shift = 0zu;
    metamath::linear::sparse_matrix_portrait<>& _portrait;
    const mesh::mesh_container_2d<T>& _mesh;

    void check_node(const size_t row, const size_t col) {
        if (_base::check(row, col) && !_included[col]) {
            _portrait.indices[_current_shift++] = col;
            _included[col] = true;
        }
    }

public:
    explicit indices_initializer(metamath::linear::sparse_matrix_portrait<>& portrait, 
                                 const mesh::mesh_container_2d<T>& mesh, const bool is_symmetric)
        : mesh::indexator_base{is_symmetric}
        , _included(portrait.cols(), false)
        , _portrait{portrait}
        , _mesh{mesh} {}

    void reset(const size_t node) override {
        _current_shift = _portrait.shifts[node];
        std::fill(std::next(_included.begin(), is_symmetric() ? node : 0), _included.end(), false);
    }

    void operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
        check_node(_mesh.node_number(e, i), _mesh.node_number(e, j));
    }

    void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
        check_node(_mesh.node_number(eL, iL), _mesh.node_number(eNL, jNL));
    }
};

}