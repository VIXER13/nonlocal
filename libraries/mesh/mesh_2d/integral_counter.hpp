#pragma once

#include "indexator_base.hpp"
#include "mesh_container_2d.hpp"

namespace nonlocal::mesh {

template<class T, class I>
class integral_counter : public indexator_base {
    std::vector<size_t>& _shifts;
    const mesh_container_2d<T, I>& _mesh;

    void check_node(const size_t row, const size_t col);

public:
    explicit integral_counter(std::vector<size_t>& shifts, const mesh_container_2d<T, I>& mesh, const bool is_symmetric);
    ~integral_counter() noexcept override = default;

    void reset(const size_t) override;

    void operator()(const std::string&, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<class T, class I>
integral_counter<T, I>::integral_counter(std::vector<size_t>& shifts, const mesh_container_2d<T, I>& mesh, const bool is_symmetric) 
    : indexator_base{is_symmetric}
    , _shifts{shifts}
    , _mesh{mesh} {}

template<class T, class I>
void integral_counter<T, I>::reset(const size_t) {}

template<class T, class I>
void integral_counter<T, I>::check_node(const size_t row, const size_t col) {
    if (!is_symmetric() || col >= row)
        ++_shifts[row];
}

template<class T, class I>
void integral_counter<T, I>::operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
    check_node(_mesh.node_number(e, i), _mesh.node_number(e, j));
}

template<class T, class I>
void integral_counter<T, I>::operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    check_node(_mesh.node_number(eL, iL), _mesh.node_number(eNL, jNL));
}

}