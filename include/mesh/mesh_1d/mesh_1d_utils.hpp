#ifndef MESH_1D_UTILS_HPP
#define MESH_1D_UTILS_HPP

#include "mesh_1d.hpp"

#include <ranges>
#include <fstream>
#include <filesystem>

namespace nonlocal::mesh::utils {

template<class T, std::ranges::random_access_range Vector>
T integrate(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The integral cannot be found because the vector size does not match the number of nodes."};
    T integral = 0;
    const auto& el = mesh.element();
    for(const size_t e : std::ranges::iota_view(size_t{0}, mesh.elements_count()))
        for(const size_t i : std::ranges::iota_view(size_t{0}, el.nodes_count()))
            for(const size_t q : std::ranges::iota_view(size_t{0}, el.qnodes_count()))
                integral += el.weight(q) * el.qN(i, q) * x[mesh.node_number(e, i)];
    return integral * mesh.jacobian();
}

template<class T, std::ranges::random_access_range Vector>
std::vector<T> approximate_gradient_in_qnodes(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The gradient cannot be found because the vector size does not match the number of nodes."};
    const auto& el = mesh.element();
    size_t qshift = 0;
    std::vector<T> gradient(el.qnodes_count() * mesh.elements_count(), T{0});
    for(const size_t e : std::ranges::iota_view(size_t{0}, mesh.elements_count()))
        for(const size_t q : std::ranges::iota_view(size_t{0}, el.qnodes_count())) {
            for(const size_t i : std::ranges::iota_view(size_t{0}, el.nodes_count()))
                gradient[qshift] += el.qNxi(i, q) * x[mesh.node_number(e, i)];
            gradient[qshift] /= mesh.jacobian();
            ++qshift;
        }
    return gradient;
}

template<class T, std::ranges::random_access_range Vector>
std::vector<T> from_qnodes_to_nodes(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    const auto& el = mesh.element();
    if (el.qnodes_count() * mesh.elements_count() != x.size())
        throw std::logic_error{"Cannot approximate node values because vector size does not match number of quadrature nodes"};
    size_t qshift = 0;
    std::vector<T> gradient(mesh.nodes_count(), T{0});
    for(const size_t e : std::ranges::iota_view(size_t{0}, mesh.elements_count())) {
        for(const size_t i : std::ranges::iota_view(size_t{0}, el.nodes_count())) {
            const size_t q = qshift + el.nearest_qnode(i);
            gradient[mesh.node_number(e, i)] += x[qshift];
        }
        qshift += el.qnodes_count();
    }
    for(const size_t node : std::ranges::iota_view(size_t{0}, mesh.nodes_count())) {
        const auto elements = mesh.node_elements(node);
        size_t count = (elements.named.curr_element != std::numeric_limits<size_t>::max()) +
                       (elements.named.next_element != std::numeric_limits<size_t>::max());
        gradient[node] /= count;
    }
    return gradient;
}

template<class T, std::ranges::random_access_range Vector, class Influence_Function>
std::vector<T> calc_flux(const mesh::mesh_1d<T>& mesh, const Vector& x,
                         const T local_weight, const Influence_Function& influnece) {
    const std::vector<T> gradient = approximate_gradient_in_qnodes(mesh, x);
    std::vector<T> flux = from_qnodes_to_nodes(mesh, gradient);
    if (local_weight < MAX_NONLOCAL_WEIGHT<T>) {
        const auto& el = mesh.element();
        std::vector<T> nonlocal_flux(flux.size(), 0);
        for(const size_t eL : std::ranges::iota_view(size_t{0}, mesh.elements_count()))
            for(const size_t i : std::ranges::iota_view(size_t{0}, el.nodes_count())) {
                const size_t node = mesh.node_number(eL, i);
                const T node_coord = mesh.node_coord(node);
                for(const size_t eNL : std::ranges::iota_view(mesh.left_neighbour(eL), mesh.right_neighbour(eL))) {
                    const size_t qshift = eNL * el.qnodes_count();
                    for(const size_t q : std::ranges::iota_view(size_t{0}, el.qnodes_count())) {
                        const T influence_weight = el.weight(q) * mesh.jacobian() * influnece(node_coord, mesh.quad_coord(eNL, q));
                        nonlocal_flux[node] += influence_weight * gradient[qshift + q];
                    }
                }
            }

        const T nonlocal_weight = T{1} - local_weight;
        for(const size_t node : std::ranges::iota_view(size_t{0}, mesh.nodes_count())) {
            const auto elements = mesh.node_elements(node);
            size_t count = (elements.named.curr_element != std::numeric_limits<size_t>::max()) +
                           (elements.named.next_element != std::numeric_limits<size_t>::max());
            flux[node] *= local_weight;
            flux[node] += nonlocal_weight * nonlocal_flux[node] / count;
        }
    }
    return flux;
}

template<class T, class Vector>
std::vector<T> gradient(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"nodes_count() != x.size()"};
    std::vector dx(x.size(), T{0});
    for(const size_t e : std::ranges::iota_view(size_t{0}, mesh.elements_count()))
        for(const size_t i : std::ranges::iota_view(size_t{0}, mesh.element().nodes_count() - (e != mesh.elements_count()-1))) {
            const size_t node = mesh.node_number(e, i);
            const T coord = mesh.node_coord(node);
            for(const size_t j : std::ranges::iota_view(size_t{0}, mesh.element().nodes_count()))
                dx[node] += mesh.element().Nxi(j, coord) * x[mesh.node_number(e, j)] / mesh.jacobian();
        }
    return dx;
}

template<class T, std::ranges::random_access_range Vector>
void save_as_csv(const mesh::mesh_1d<T>& mesh, const Vector& x, const std::filesystem::path& path) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The result cannot be saved because the mesh nodes number and elements in the vector do not match."};
    std::ofstream csv{path};
    for(const size_t segment : std::ranges::iota_view{size_t{0}, mesh.segments_count()}) {
        const auto nodes = mesh.segment_nodes(segment);
        for(const size_t node : std::ranges::iota_view{nodes.front(), nodes.back()}) // exclude last node
            csv << mesh.node_coord(node) << ',' << x[node] << '\n';
    }
    csv << mesh.length() << ',' << x[x.size() - 1] << '\n';
}

}

#endif