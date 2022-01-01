#ifndef MESH_SU2_PARSER_HPP
#define MESH_SU2_PARSER_HPP

#include "mesh_container_2d.hpp"
#include <iostream>
#include <ranges>

namespace nonlocal::mesh {

template<class T, class I>
template<size_t... K, class Stream>
void mesh_2d<T, I>::read_element(Stream& mesh_file, std::vector<I>& element) {
    element.resize(sizeof...(K));
    (mesh_file >> ... >> element[K]);
}

template<class T, class I>
template<class Stream>
void mesh_2d<T, I>::read_su2(Stream& mesh_file) {
    std::string pass;
    size_t count = 0;

    // 2D элементы
    mesh_file >> pass >> pass >> pass >> count;
    _elements.resize(count);
    _elements_2d_type.resize(count);
    for(size_t e = 0; e < count; ++e) {
        uintmax_t type = 0;
        mesh_file >> type;
        switch(vtk_element_number(type)) {
            case vtk_element_number::TRIANGLE:
                _elements_2d_type[e] = element_2d_t::TRIANGLE;
                read_element<0, 1, 2>(mesh_file, _elements[e]);
            break;

            case vtk_element_number::QUADRATIC_TRIANGLE:
               _elements_2d_type[e] = element_2d_t::QUADRATIC_TRIANGLE;
               read_element<0, 1, 2, 3, 4, 5>(mesh_file, _elements[e]);
            break;

            case vtk_element_number::BILINEAR:
                _elements_2d_type[e] = element_2d_t::BILINEAR;
                read_element<0, 1, 2, 3>(mesh_file, _elements[e]);
            break;

            case vtk_element_number::QUADRATIC_SERENDIPITY:
                _elements_2d_type[e] = element_2d_t::QUADRATIC_SERENDIPITY;
                read_element<0, 2, 4, 6, 1, 3, 5, 7>(mesh_file, _elements[e]);
            break;

            case vtk_element_number::QUADRATIC_LAGRANGE:
               _elements_2d_type[e] = element_2d_t::QUADRATIC_LAGRANGE;
               read_element<0, 2, 4, 6, 1, 3, 5, 7, 8>(mesh_file, _elements[e]);
            break;

            default:
                throw std::domain_error{"Unknown 2D element."};
        }
        mesh_file >> pass;
    }

    // Узлы
    mesh_file >> pass >> count;
    _nodes.resize(count);
    for(size_t node = 0; node < count; ++node)
        mesh_file >> _nodes[node][0] >> _nodes[node][1] >> pass;

    std::string name;
    mesh_file >> pass >> count;
    for(const size_t b : std::views::iota(size_t{0}, count)) {
        mesh_file >> pass >> name >> pass >> count;
        _boundaries_names.template emplace_back(std::move(name));
        auto& boundary = _boundaries[_boundaries_names.back()];
        auto& boundary_elements_type = _elements_1d_type[_boundaries_names.back()];
        boundary.resize(count);
        boundary_elements_type.resize(count);
        for(const size_t e : std::views::iota(size_t{0}, count)) {
            uintmax_t type = 0;
            mesh_file >> type;
            switch(vtk_element_number(type)) {
                case vtk_element_number::LINEAR:
                    boundary_elements_type[e] = element_1d_t::LINEAR;
                    read_element<0, 1>(mesh_file, boundary[e]);
                break;

                case vtk_element_number::QUADRATIC:
                    boundary_elements_type[e] = element_1d_t::QUADRATIC;
                    read_element<0, 2, 1>(mesh_file, boundary[e]);
                break;

                default:
                    throw std::domain_error{"Unknown 1D element."};
            }
        }
    }
}

}

#endif