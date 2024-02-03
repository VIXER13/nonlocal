#ifndef NONLOCAL_NONZERO_COUNTER_HPP
#define NONLOCAL_NONZERO_COUNTER_HPP

#include "indexator_base.hpp"
#include "mesh_container_2d.hpp"

namespace nonlocal::mesh {

template<class T, class I>
class nonzero_counter : public indexator_base {
    std::vector<bool> _included;
    std::vector<size_t>& _shifts;
    const mesh_container_2d<T, I>& _mesh;

    void check_node(const size_t row, const size_t col);

public:
    explicit nonzero_counter(std::vector<size_t>& shifts, const mesh_container_2d<T, I>& mesh, const bool is_symmetric);
    ~nonzero_counter() noexcept override = default;

    void reset(const size_t node) override;

    void operator()(const std::string&, const size_t e, const size_t i, const size_t j);
    void operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL);
};

template<class T, class I>
nonzero_counter<T, I>::nonzero_counter(std::vector<size_t>& shifts, const mesh_container_2d<T, I>& mesh, const bool is_symmetric)
    : indexator_base{is_symmetric}
    , _included(shifts.size(), false)
    , _shifts{shifts}
    , _mesh{mesh} {}

template<class T, class I>
void nonzero_counter<T, I>::reset(const size_t node) {
    std::fill(std::next(_included.begin(), is_symmetric() ? node : 0), _included.end(), false);
}

template<class T, class I>
void nonzero_counter<T, I>::check_node(const size_t row, const size_t col) {
    if ((!is_symmetric() || col >= row) && !_included[col]) {
        _included[col] = true;
        ++_shifts[row];
    }
}

template<class T, class I>
void nonzero_counter<T, I>::operator()(const std::string&, const size_t e, const size_t i, const size_t j) {
    check_node(_mesh.node_number(e, i), _mesh.node_number(e, j));
}

template<class T, class I>
void nonzero_counter<T, I>::operator()(const std::string&, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
    check_node(_mesh.node_number(eL, iL), _mesh.node_number(eNL, jNL));
}

}

#endif