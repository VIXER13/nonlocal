#pragma once

#include "indexator_base.hpp"
#include "mesh_container_2d.hpp"

namespace nonlocal::mesh {

template<std::floating_point T, std::integral I, std::integral S = size_t>
class nonzero_counter : public indexator_base {
    using _base = indexator_base;

    std::vector<bool> _included;
    std::vector<S>& _shifts;
    const mesh_container_2d<T, I>& _mesh;

    void check_node(const size_t row, const size_t col) {
        if (_base::check(row, col) && !_included[col]) {
            _included[col] = true;
            ++_shifts[row + 1];
        }
    }

public:
    explicit nonzero_counter(std::vector<S>& shifts, const mesh_container_2d<T, I>& mesh, const bool is_symmetric)
        : indexator_base{is_symmetric}
        , _included(shifts.size(), false)
        , _shifts{shifts}
        , _mesh{mesh} {}
    ~nonzero_counter() noexcept override = default;

    void reset(const size_t node) override {
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