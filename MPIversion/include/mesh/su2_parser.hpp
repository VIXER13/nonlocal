#ifndef MESH_SU2_PARSER_HPP
#define MESH_SU2_PARSER_HPP

#include "mesh_2d.hpp"
#include <iostream>

namespace mesh {

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
    for(size_t el = 0; el < count; ++el) {
        uintmax_t type = 0;
        mesh_file >> type;
        switch(vtk_element_number(type)) {
            case vtk_element_number::TRIANGLE:
                _elements_2d_type[el] = element_2d_t::TRIANGLE;
                read_element<0, 1, 2>(mesh_file, _elements[el]);
            break;

            case vtk_element_number::QUADRATIC_TRIANGLE:
               _elements_2d_type[el] = element_2d_t::QUADRATIC_TRIANGLE;
               read_element<0, 1, 2, 3, 4, 5>(mesh_file, _elements[el]);
            break;

            case vtk_element_number::BILINEAR:
                _elements_2d_type[el] = element_2d_t::BILINEAR;
                read_element<0, 1, 2, 3>(mesh_file, _elements[el]);
            break;

            case vtk_element_number::QUADRATIC_SERENDIPITY:
                _elements_2d_type[el] = element_2d_t::QUADRATIC_SERENDIPITY;
                read_element<0, 2, 4, 6, 1, 3, 5, 7>(mesh_file, _elements[el]);
            break;

            case vtk_element_number::QUADRATIC_LAGRANGE:
               _elements_2d_type[el] = element_2d_t::QUADRATIC_LAGRANGE;
               read_element<0, 2, 4, 6, 1, 3, 5, 7, 8>(mesh_file, _elements[el]);
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

    // Элементы на границах
    mesh_file >> pass >> count;
    _boundaries.resize(count);
    _elements_1d_type.resize(count);
    for(size_t b = 0; b < _boundaries.size(); ++b) {
        mesh_file >> pass >> pass;
        if(pass != "Group_Of_All_Edges") {
            mesh_file >> pass >> count;
            _boundaries[b].resize(count);
            _elements_1d_type[b].resize(count);
            for(size_t el = 0; el < count; ++el) {
                uintmax_t type = 0;
                mesh_file >> type;
                switch(vtk_element_number(type)) {
                    case vtk_element_number::LINEAR:
                        _elements_1d_type[b][el] = element_1d_t::LINEAR;
                        read_element<0, 1>(mesh_file, _boundaries[b][el]);
                    break;

                    case vtk_element_number::QUADRATIC:
                        _elements_1d_type[b][el] = element_1d_t::QUADRATIC;
                        read_element<0, 2, 1>(mesh_file, _boundaries[b][el]);
                    break;

                    default:
                        throw std::domain_error{"Unknown 1D element."};
                }
            }
        } else {
            _boundaries.pop_back();
            _elements_1d_type.pop_back();
        }
    }
}

}

#endif