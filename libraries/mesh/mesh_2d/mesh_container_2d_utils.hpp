#pragma once

#include "su2_parser.hpp"

#include "nonlocal_constants.hpp"

namespace nonlocal::mesh::utils {

template<class T, class I>
std::vector<std::vector<I>> node_elements_2d(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::vector<I>> node_elements(mesh.nodes_count());
    for(const size_t e : mesh.elements_2d()) {
        for(const size_t node : mesh.nodes(e))
            node_elements[node].push_back(e);
    }
    for(std::vector<I>& elements : node_elements)
        elements.shrink_to_fit();
    return node_elements;
}

template<class T, class I>
std::vector<std::unordered_map<I, uint8_t>> global_to_local(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering(mesh.elements_2d_count() + mesh.elements_1d_count());
#pragma omp parallel for default(none) shared(mesh, global_to_local_numbering)
    for(size_t e = 0; e < global_to_local_numbering.size(); ++e) {
        const std::vector<I>& nodes = mesh.nodes(e);
        for(const size_t i : std::ranges::iota_view{0u, nodes.size()})
            global_to_local_numbering[e][nodes[i]] = i;
    }
    return global_to_local_numbering;
}

template<class T, class I, class Shift>
std::vector<I> quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh, const Shift& shift) {
    std::vector<I> quad_shifts(mesh.elements_2d_count() + 1);
    quad_shifts[0] = 0;
    for(const size_t e : mesh.elements_2d())
        quad_shifts[e + 1] = quad_shifts[e] + shift(e);
    return quad_shifts;
}

template<class T, class I>
std::vector<I> elements_quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh) {
    return quadrature_shifts_2d(mesh, [&mesh](const size_t e) { return mesh.element_2d(e).qnodes_count(); });
}

template<class T, class I>
std::vector<I> element_node_shits_quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh) {
    return quadrature_shifts_2d(mesh, [&mesh](const size_t e) { 
        const auto& el = mesh.element_2d(e);
        return el.nodes_count() * el.qnodes_count();
    });
}

template<template<class, size_t> class Output, class T, class I, class Functor>
std::vector<Output<T, 2>> approx_in_all_quad_nodes(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts, const Functor& functor) {
    if(mesh.elements_2d_count() + 1 != qshifts.size())
        throw std::logic_error{"The number of quadrature shifts and elements does not match."};
    std::vector<Output<T, 2>> data(qshifts.back());
    for(const size_t e : mesh.elements_2d()) {
        const auto element_data = mesh.element_2d_data(e);
        for(const size_t q : std::ranges::iota_view{0u, element_data.element.qnodes_count()})
            data[qshifts[e] + q] = functor(element_data, q);
    }
    return data;
}

template<class T, class I>
std::vector<std::array<T, 2>> approx_all_quad_nodes(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts) {
    return approx_in_all_quad_nodes<std::array>(mesh, qshifts, 
        [](const auto& element_data, const size_t q) { return element_data.quad_coord(q); });
}

template<class T, class I>
std::vector<metamath::types::square_matrix<T, 2>> approx_all_jacobi_matrices(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts) {
    return approx_in_all_quad_nodes<metamath::types::square_matrix>(mesh, qshifts, 
        [](const auto& element_data, const size_t q) { return element_data.jacobi_matrix(q); });
}

template<class T, class I>
std::vector<std::array<T, 2>> derivatives_in_quad(const mesh_container_2d<T, I>& mesh,
                                                  const std::vector<I>& quad_element_shifts,
                                                  const std::vector<I>& quad_nodes_shifts,
                                                  const std::vector<metamath::types::square_matrix<T, 2>>& jacobi_matrices) {
    if (mesh.elements_2d_count() + 1 != quad_element_shifts.size() || mesh.elements_2d_count() + 1 != quad_nodes_shifts.size())
        throw std::logic_error{"The number of quadrature shifts and elements does not match."};
    if (quad_element_shifts.back() != jacobi_matrices.size())
        throw std::logic_error{"The size of Jacobi matrices vector does not match with the quadratures nodes count."};
    std::vector<std::array<T, 2>> derivatives(quad_nodes_shifts.back());
#pragma omp parallel for default(none) shared(mesh, quad_element_shifts, quad_nodes_shifts, jacobi_matrices, derivatives)
    for(size_t e = 0; e < mesh.elements_2d_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()})
            for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
                const metamath::types::square_matrix<T, 2>& J = jacobi_matrices[quad_element_shifts[e] + q];
                derivatives[quad_nodes_shifts[e] + i * el.qnodes_count() + q] = {
                     el.qNxi(i, q) * J[1][1] - el.qNeta(i, q) * J[1][0],
                    -el.qNxi(i, q) * J[0][1] + el.qNeta(i, q) * J[0][0]
                };
            }
    }
    return derivatives;
}

template<class T, class I>
std::vector<std::array<T, 2>> approx_centers_of_elements(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::array<T, 2>> centers(mesh.elements_2d_count(), std::array<T, 2>{});
#pragma omp parallel for default(none) shared(centers, mesh)
    for(size_t e = 0; e < centers.size(); ++e)
        centers[e] = mesh.element_2d_data(e).center();
    return centers;
}

template<class I, class T>
std::vector<I> find_neighbours(const T radius, const std::vector<std::array<T, 2>>& nodes, const size_t node) {
    std::vector<I> neighbours;
    for(const size_t i : std::ranges::iota_view{0u, nodes.size()})
        if (metamath::functions::distance(nodes[node], nodes[i]) <= radius)
            neighbours.push_back(i);
    neighbours.shrink_to_fit();
    return neighbours;
}

template<class I, class T>
std::vector<std::vector<I>> find_neighbours(const T radius, const std::vector<std::array<T, 2>>& nodes) {
    std::vector<std::vector<I>> neighbours(nodes.size());
#pragma omp parallel for default(none) shared(neighbours, nodes)
    for(size_t i = 0; i < nodes.size(); ++i)
        neighbours[i] = find_neighbours<I>(radius, nodes, i);
    return neighbours;
}

template<class T, class I>
std::vector<std::vector<I>> find_neighbours(const T radius, const std::vector<std::array<T, 2>>& nodes, const std::unordered_set<I>& nodes_to_search) {
    std::vector<std::vector<I>> neighbours(nodes.size());
    for(const size_t node : nodes_to_search)
        neighbours[node] = find_neighbours<I>(radius, nodes, node);
    return neighbours;
}

template<class Stream, class T, class I>
void save_as_vtk(Stream& stream, const mesh_container_2d<T, I>& mesh) {
    static constexpr auto write_element = []<size_t K0, size_t... K>(Stream& stream, const std::vector<I>& element, const std::index_sequence<K0, K...>) {
        stream << element[K0];
        ((stream << ' ' << element[K]), ...);
    };

    stream << "# vtk DataFile Version 4.2\n"
              "Data\n"
              "ASCII\n"
              "DATASET UNSTRUCTURED_GRID\n";

    stream << "POINTS " << mesh.nodes_count() << ' ' << vtk_data_type<T> << '\n';
    for(const size_t i : mesh.nodes()) {
        const std::array<T, 2>& point = mesh.node_coord(i);
        stream << point[X] << ' ' << point[Y] << " 0\n";
    }

    const auto elements_2d = mesh.elements_2d();
    const auto reduces = [&mesh](const size_t sum, const size_t e) { return sum + mesh.nodes_count(e) + 1; };
    const size_t list_size = std::reduce(elements_2d.begin(), elements_2d.end(), size_t{0}, reduces);

    stream << "CELLS " << mesh.elements_2d_count() << ' ' << list_size << '\n';
    for(const size_t e : elements_2d) {
        stream << mesh.nodes_count(e) << ' ';
        switch (mesh.element_type_2d(e)) {
            case element_2d_t::TRIANGLE:
                write_element(stream, mesh.nodes(e), std::index_sequence<0, 1, 2>{});
            break;

            case element_2d_t::QUADRATIC_TRIANGLE:
                write_element(stream, mesh.nodes(e), std::index_sequence<0, 1, 2, 3, 4, 5>{});
            break;

            case element_2d_t::BILINEAR:
                write_element(stream, mesh.nodes(e), std::index_sequence<0, 1, 2, 3>{});
            break;

            case element_2d_t::QUADRATIC_SERENDIPITY:
                write_element(stream, mesh.nodes(e), std::index_sequence<0, 2, 4, 6, 1, 3, 5, 7>{});
            break;

            case element_2d_t::QUADRATIC_LAGRANGE:
                write_element(stream, mesh.nodes(e), std::index_sequence<0, 2, 4, 6, 1, 3, 5, 7, 8>{});
            break;
            
            default:
                throw std::domain_error{"Unknown element."};
        }
        stream << '\n';
    }

    stream << "CELL_TYPES " << mesh.elements_2d_count() << '\n';
    for(const size_t e : elements_2d)
        stream << size_t(mesh.get_elements_set().local_to_model_2d(mesh.element_type_2d(e))) << '\n';
}

template<class T, class I>
void save_as_vtk(const std::filesystem::path& path_to_save, const mesh_container_2d<T, I>& mesh) {
    std::ofstream vtk{path_to_save};
    save_as_vtk(vtk, mesh);
}

template<class T>
void save_scalars_to_vtk(std::ofstream& output, const std::string_view name, const std::vector<T>& x) {
    output << "SCALARS " << name << ' ' << mesh::vtk_data_type<T> << " 1\n"
           << "LOOKUP_TABLE default\n";
    for(const T val : x)
        output << val << '\n';
}

template<class T>
void save_vectors_to_vtk(std::ofstream& output, const std::string_view name, const std::array<std::vector<T>, 2>& vector) {
    output << "VECTORS " << name << ' ' << mesh::vtk_data_type<T> << '\n';
    for(const size_t i : std::ranges::iota_view{0u, vector[X].size()})
        output << vector[X][i] << ' ' << vector[Y][i] << " 0\n";
}

template<class T>
void save_tensors_to_vtk(std::ofstream& output, const std::string_view name, const std::array<std::vector<T>, 3>& tensor) {
    output << "TENSORS " << name << ' ' << mesh::vtk_data_type<T> << '\n';
    for(const size_t i : std::ranges::iota_view{0u, tensor[0].size()})
        output << tensor[0][i] << ' ' << tensor[2][i] << " 0\n"
               << tensor[2][i] << ' ' << tensor[1][i] << " 0\n"
               << "0 0 0\n\n";
}

template<class T, class I>
void save_as_csv(const std::filesystem::path& path, const mesh_container_2d<T, I>& mesh,
                 const std::vector<std::pair<std::string, const std::vector<T>&>>& data,
                 const std::optional<std::streamsize> precision = std::nullopt) {
    for(const auto& [name, vec] : data)
        if (mesh.nodes_count() != vec.size())
            throw std::logic_error{"The result cannot be saved because the mesh nodes number "
                                   "and elements in the vector \"" + name + "\" do not match."};
    std::ofstream csv{path};
    csv.precision(precision ? *precision : std::numeric_limits<T>::max_digits10);
    csv << "x,y" << (data.empty() ? '\n' : ',');
    for(const size_t j : std::ranges::iota_view{0u, data.size()})
        csv << data[j].first << (j == data.size() - 1 ? '\n' : ',');
    for(const size_t i : std::ranges::iota_view{0u, mesh.nodes_count()}) {
        const std::array<T, 2>& node = mesh.node_coord(i);
        csv << node[X] << ',' << node[Y] << (data.empty() ? '\n' : ',');
        for(const size_t j : std::ranges::iota_view{0u, data.size()})
            csv << data[j].second[i] << (j == data.size() - 1 ? '\n' : ',');
    }
}

}