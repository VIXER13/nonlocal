#pragma once

#include "problem_settings.hpp"
#include "indices_initializer.hpp"
#include "integrator.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/nonzero_counter.hpp>
#include <metamath/linear/sparse_matrix.hpp>

namespace nonlocal::solver_2d {

template<class T>
class matrix_assembler_base {
    using entity_t = metamath::types::container_type_t<T>;
    using floating_point_t = metamath::types::container_type_t<entity_t>;
    using nodes_sequence = std::variant<
        std::ranges::iota_view<size_t, size_t>,
        std::vector<size_t>
    >;

    metamath::linear::sparse_matrix<T> _matrix;
    const mesh::mesh_2d<floating_point_t>& _mesh;

    template<class Runner>
    void mesh_run(const problem_settings& settings, Runner&& runner);

protected:
    explicit matrix_assembler_base(const mesh::mesh_2d<floating_point_t>& mesh);

    size_t rows() const noexcept;

    void init_shifts(const problem_settings& settings);
    void init_indices(const problem_settings& settings, const bool sort_indices = true);
    template<class Local_Integrator, class Nonlocal_Integrator>
    void calc_coeffs(const problem_settings& settings, Local_Integrator&& local_integrator, Nonlocal_Integrator&& nonlocal_integrator);

public:
    nodes_sequence processing_nodes;

    virtual ~matrix_assembler_base() noexcept = default;

    const mesh::mesh_2d<floating_point_t>& mesh() const noexcept;
    metamath::linear::sparse_matrix<T>& matrix() noexcept;
    const metamath::linear::sparse_matrix<T>& matrix() const noexcept;
};

template<class T>
matrix_assembler_base<T>::matrix_assembler_base(const mesh::mesh_2d<floating_point_t>& mesh)
    :  _mesh{mesh}
    , processing_nodes{_mesh.process_nodes()} {}

template<class T>
const mesh::mesh_2d<typename matrix_assembler_base<T>::floating_point_t>& matrix_assembler_base<T>::mesh() const noexcept {
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
    mesh_run(settings, integrator{matrix(), mesh().container(), settings.is_symmetric(), std::move(local_integrator), std::move(nonlocal_integrator)});
}

}