#pragma once

#include "problem_settings.hpp"
#include "indices_initializer.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/nonzero_counter.hpp>
#include <metamath/linear/sparse_matrix.hpp>

namespace nonlocal::solver_2d {

namespace details {

// TODO: temporary here, but after refactoring will be moved in to integrator.hpp

template<std::floating_point T, class Local_Integrator, class Nonlocal_Integrator>
class integrator final : public mesh::indexator_base {
    using _base = mesh::indexator_base;

    Local_Integrator _local_integrator;
    Nonlocal_Integrator _nonlocal_integrator;
    metamath::linear::sparse_matrix<T>& _matrix;
    const mesh::mesh_container_2d<T>& _mesh;

public:
    explicit integrator(metamath::linear::sparse_matrix<T>& matrix, const mesh::mesh_container_2d<T>& mesh, const bool is_symmetric,
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

template<class T>
class matrix_assembler_base {
    using nodes_sequence = std::variant<
        std::ranges::iota_view<size_t, size_t>,
        std::vector<size_t>
    >;

    metamath::linear::sparse_matrix<T> _matrix;
    const mesh::mesh_2d<T>& _mesh;

    template<class Runner>
    void mesh_run(const problem_settings& settings, Runner&& runner);

protected:
    explicit matrix_assembler_base(const mesh::mesh_2d<T>& mesh);

    size_t rows() const noexcept;

    void init_shifts(const problem_settings& settings);
    void init_indices(const problem_settings& settings, const bool sort_indices = true);
    template<class Local_Integrator, class Nonlocal_Integrator>
    void calc_coeffs(const problem_settings& settings, Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator);

public:
    nodes_sequence processing_nodes;

    virtual ~matrix_assembler_base() noexcept = default;

    const mesh::mesh_2d<T>& mesh() const noexcept;
    metamath::linear::sparse_matrix<T>& matrix() noexcept;
    const metamath::linear::sparse_matrix<T>& matrix() const noexcept;
};

template<class T>
matrix_assembler_base<T>::matrix_assembler_base(const mesh::mesh_2d<T>& mesh)
    :  _mesh{mesh}
    , processing_nodes{_mesh.process_nodes()} {}

template<class T>
const mesh::mesh_2d<T>& matrix_assembler_base<T>::mesh() const noexcept {
    return _mesh;
}

template<class T>
metamath::linear::sparse_matrix<T>& matrix_assembler_base<T>::matrix() noexcept {
    return _matrix;
}

template<class T>
const metamath::linear::sparse_matrix<T>& matrix_assembler_base<T>::matrix() const noexcept {
    return _matrix;
}

template<class T>
size_t matrix_assembler_base<T>::rows() const noexcept {
    return std::visit([](const auto& nodes) { return nodes.size(); }, processing_nodes);
}

template<class T>
template<class Runner>
void matrix_assembler_base<T>::mesh_run(const problem_settings& settings, Runner&& runner) {
    std::visit([this, &settings, &runner](const auto& nodes) {
        mesh::utils::mesh_run(mesh(), nodes, settings.theories, std::move(runner));
    }, processing_nodes);
}

template<class T>
void matrix_assembler_base<T>::init_shifts(const problem_settings& settings) {
    mesh_run(settings, mesh::nonzero_counter{matrix().portrait.shifts, mesh().container(), settings.is_symmetric()});
    matrix().portrait.accumulate_shifts();
    logger::info() << "Non-zero elements count: " << matrix().non_zeros() << std::endl;
}

template<class T>
void matrix_assembler_base<T>::init_indices(const problem_settings& settings, const bool sort_indices) {
    matrix().portrait.allocate_indices();
    matrix().allocate_values();
    logger::info() << "Matrix allocated successfully" << std::endl;
    mesh_run(settings, indices_initializer{matrix().portrait, mesh().container(), settings.is_symmetric()});
    if (sort_indices)
        matrix().portrait.sort_indices();
}

template<class T>
template<class Local_Integrator, class Nonlocal_Integrator>
void matrix_assembler_base<T>::calc_coeffs(const problem_settings& settings, Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator) {
    mesh_run(settings, details::integrator{matrix(), mesh().container(), settings.is_symmetric(), std::move(local_integrator), std::move(nonlocal_integrator)});
}

}