#ifndef SAVE_RESULTS_HPP
#define SAVE_RESULTS_HPP

#include <string_view>
#include "mesh_2d.hpp"

namespace mesh {

template<class T, class I>
template<size_t K0, size_t... K>
void mesh_2d<T, I>::write_element(std::ofstream& mesh_file, const std::vector<I>& element) const {
    mesh_file << element[K0];
    ((mesh_file << ' ' << element[K]), ...);
}

template<class T, class I>
void mesh_2d<T, I>::save_as_vtk(std::ofstream& mesh_file) const {
    static constexpr std::string_view data_type = std::is_same_v<T, float> ? "float" : "double";

    mesh_file << "# vtk DataFile Version 4.2\n"
                 "Data\n"
                 "ASCII\n"
                 "DATASET UNSTRUCTURED_GRID\n";

    mesh_file << "POINTS " << nodes_count() << ' ' << data_type << '\n';
    for(size_t i = 0; i < nodes_count(); ++i) {
        const std::array<T, 2>& point = node(i);
        mesh_file << point[0] << ' ' << point[1] << " 0\n";
    }

    static constexpr auto accumulator = [](const uintmax_t sum, const std::vector<I>& element) { return sum + element.size() + 1; };
    const uintmax_t list_size = std::accumulate(_elements.cbegin(), _elements.cend(), uintmax_t{0}, accumulator);
    mesh_file << "CELLS " << elements_count() << ' ' << list_size << '\n';
    for(size_t e = 0; e < elements_count(); ++e) {
        const auto& el = element_2d(element_2d_type(e));
        mesh_file << el->nodes_count() << ' ';
        switch (element_2d_type(e)) {
            case element_2d_t::TRIANGLE:
                write_element<0, 1, 2>(mesh_file, _elements[e]);
            break;

            case element_2d_t::QUADRATIC_TRIANGLE:
                write_element<0, 1, 2, 3, 4, 5>(mesh_file, _elements[e]);
            break;

            case element_2d_t::BILINEAR:
                write_element<0, 1, 2, 3>(mesh_file, _elements[e]);
            break;

            case element_2d_t::QUADRATIC_SERENDIPITY:
                write_element<0, 2, 4, 6, 1, 3, 5, 7>(mesh_file, _elements[e]);
            break;

            case element_2d_t::QUADRATIC_LAGRANGE:
                write_element<0, 2, 4, 6, 1, 3, 5, 7, 8>(mesh_file, _elements[e]);
            break;
            
            default:
                throw std::domain_error{"Unknown element."};
        }
        mesh_file << '\n';
    }

    mesh_file << "CELL_TYPES " << elements_count() << '\n';
    for(size_t e = 0; e < elements_count(); ++e) {
        switch (element_2d_type(e)) {
            case element_2d_t::TRIANGLE:
                mesh_file << uintmax_t(vtk_element_number::TRIANGLE) << '\n';
            break;

            case element_2d_t::QUADRATIC_TRIANGLE:
                mesh_file << uintmax_t(vtk_element_number::QUADRATIC_TRIANGLE) << '\n';
            break;

            case element_2d_t::BILINEAR:
                mesh_file << uintmax_t(vtk_element_number::BILINEAR) << '\n';
            break;

            case element_2d_t::QUADRATIC_SERENDIPITY:
                mesh_file << uintmax_t(vtk_element_number::QUADRATIC_SERENDIPITY) << '\n';
            break;

            case element_2d_t::QUADRATIC_LAGRANGE:
                mesh_file << uintmax_t(vtk_element_number::QUADRATIC_LAGRANGE) << '\n';
            break;
            
            default:
                throw std::domain_error{"Unknown element."};
        }
    }
}

}

#endif