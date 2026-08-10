#pragma once

#include <mesh/mesh_2d/indexator_base.hpp>

namespace nonlocal::solver_2d {

template<class T, class Local_Integrator, class Nonlocal_Integrator>
class integrator final : public mesh::indexator_base {
    using _base = mesh::indexator_base;
    using entity_t = metamath::types::container_type_t<T>;
    using floating_point_t = metamath::types::container_type_t<entity_t>;

    Local_Integrator _local_integrator;
    Nonlocal_Integrator _nonlocal_integrator;
    metamath::linear::sparse_matrix<T>& _matrix;
    const mesh::mesh_container_2d<floating_point_t>& _mesh;

public:
    explicit integrator(metamath::linear::sparse_matrix<T>& matrix,
                        const mesh::mesh_container_2d<floating_point_t>& mesh, const bool is_symmetric,
                        Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator)
        : _base{is_symmetric}
        , _matrix{matrix}
        , _mesh{mesh}
        , _local_integrator{std::move(local_integrator)}
        , _nonlocal_integrator{std::move(nonlocal_integrator)} {}

    void reset(const size_t node) override {}

    void operator()(const std::string& group, const size_t e, const size_t i, const size_t j) {
        const size_t row = _mesh.node_number(e, i);
        const size_t col = _mesh.node_number(e, j);
        using namespace metamath::operators;
        if (_base::check(row, col))
            _matrix(row, col) += _local_integrator(group, e, i, j);
    }

    void operator()(const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
        const size_t row = _mesh.node_number(eL, iL);
        const size_t col = _mesh.node_number(eNL, jNL);
        if (_base::check(row, col)) {
            using namespace metamath::operators;
            auto& value = _matrix(row, col) += _nonlocal_integrator(group, eL, eNL, iL, jNL);
            if (eL == eNL)
                value += _local_integrator(group, eL, iL, jNL);
        }
    }
};

}