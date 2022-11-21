#ifndef NONLOCAL_MESH_1D_UTILS_HPP
#define NONLOCAL_MESH_1D_UTILS_HPP

#include "mesh_1d.hpp"

#include <fstream>
#include <filesystem>

namespace nonlocal::mesh::utils {

template<class T, class Vector>
T integrate(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The integral cannot be found because the vector size does not match the number of nodes."};
    T integral = T{0};
    const auto& el = mesh.element();
    for(const size_t segment : mesh.segments()) {
        T segment_integral = T{0};
        for(const size_t e : mesh.elements(segment))
            for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()})
                for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
                    segment_integral += el.weight(q) * el.qN(i, q) * x[mesh.node_number(e, i)];
        integral += segment_integral * mesh.jacobian(segment);
    }
    return integral;
}

template<class T, class Vector>
std::vector<T> gradient_in_qnodes(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    const auto& el = mesh.element();
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The gradient cannot be found because the vector size does not match the number of nodes."};
    size_t qshift = 0;
    std::vector<T> gradient(el.qnodes_count() * mesh.elements_count(), T{0});
    for(const size_t segment : mesh.segments()) {
        const T jacobian = mesh.jacobian(segment);
        for(const size_t e : mesh.elements(segment))
            for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
                for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()})
                    gradient[qshift] += el.qNxi(i, q) * x[mesh.node_number(e, i)];
                gradient[qshift] /= jacobian;
                ++qshift;
            }
    }
    return gradient;
}

template<class T, class Vector>
std::vector<T> from_qnodes_to_nodes(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    const auto& el = mesh.element();
    if (el.qnodes_count() * mesh.elements_count() != x.size())
        throw std::logic_error{"Cannot approximate node values because vector size does not match number of quadrature nodes"};
    size_t qshift = 0;
    std::vector<T> gradient(mesh.nodes_count(), T{0});
    for(const size_t segment : mesh.segments())
        for(const size_t e : mesh.elements(segment)) {
            for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()})
                gradient[mesh.node_number(e, i)] += x[qshift + el.nearest_qnode(i)];
            qshift += el.qnodes_count();
        }
    for(const size_t node : mesh.nodes())
        gradient[node] /= mesh.node_elements(node).count();
    return gradient;
}

template<class T, std::ranges::random_access_range Vector>
void save_as_csv(const mesh::mesh_1d<T>& mesh, const Vector& x, const std::filesystem::path& path) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The result cannot be saved because the mesh nodes number and elements in the vector do not match."};
    std::ofstream csv{path};
    for(const size_t segment : std::ranges::iota_view{0u, mesh.segments_count()}) {
        const auto nodes = mesh.nodes(segment);
        for(const size_t node : std::ranges::iota_view{nodes.front(), nodes.back()}) // exclude last node
            csv << mesh.node_coord(node) << ',' << x[node] << '\n';
    }
    csv << mesh.length() << ',' << x[x.size() - 1] << '\n';
}

}

#endif