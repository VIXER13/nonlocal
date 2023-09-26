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

template<class T, class Callback>
std::vector<T> approximation_in_qnodes(const mesh::mesh_1d<T>& mesh, const Callback& callback) {
    size_t qshift = 0;
    std::vector<T> approximation(mesh.element().qnodes_count() * mesh.elements_count(), T{0});
    for(const size_t segment : mesh.segments()) {
        const T jacobian = mesh.jacobian(segment);
        for(const size_t e : mesh.elements(segment))
            for(const size_t q : mesh.element().qnodes()) {
                for(const size_t i : mesh.element().nodes())
                    approximation[qshift] += callback(e, i, q);
                ++qshift;
            }
    }
    return approximation;
}

template<class T, class Vector>
std::vector<T> gradient_in_qnodes(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"The gradient approximation cannot be found because the vector size does not match the number of nodes."};
    return approximation_in_qnodes(mesh, [&mesh, &x](const size_t e, const size_t i, const size_t q) {
        const size_t segment = mesh.segment_number(e);
        const T jacobian = mesh.jacobian(segment);
        return mesh.element().qNxi(i, q) * x[mesh.node_number(e, i)] / jacobian;
    });
}

template<class T, class Vector>
std::vector<T> from_nodes_to_qnodes(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() > x.size())
        throw std::logic_error{"The approximation cannot be found because the vector size does not match the number of nodes."};
    return approximation_in_qnodes(mesh, [&mesh, &x](const size_t e, const size_t i, const size_t q) {
        return mesh.element().qN(i, q) * x[mesh.node_number(e, i)];
    });
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

template<class T>
void save_as_csv(const std::filesystem::path& path, const mesh::mesh_1d<T>& mesh,
                 const std::vector<std::pair<std::string, const std::vector<T>&>>& data,
                 const std::optional<std::streamsize> precision = std::nullopt) {
    for(const auto& [name, vec] : data)
        if (mesh.nodes_count() != vec.size())
            throw std::logic_error{"The result cannot be saved because the mesh nodes number "
                                   "and elements in the vector \"" + name + "\" do not match."};
    std::ofstream csv{path};
    csv.precision(precision ? *precision : std::numeric_limits<T>::max_digits10);
    csv << 'x' << (data.empty() ? '\n' : ',');
    for(const size_t j : std::ranges::iota_view{0u, data.size()})
        csv << data[j].first << (j == data.size() - 1 ? '\n' : ',');
    for(const size_t i : std::ranges::iota_view{0u, mesh.nodes_count()}) {
        csv << mesh.node_coord(i) << (data.empty() ? '\n' : ',');
        for(const size_t j : std::ranges::iota_view{0u, data.size()})
            csv << data[j].second[i] << (j == data.size() - 1 ? '\n' : ',');
    }
}

}

#endif