#ifndef MESH_NEAREST_NEIGHBOR_HPP
#define MESH_NEAREST_NEIGHBOR_HPP

#include "utils.hpp"
#include <iostream>

namespace mesh {

// Ищет соседние элементы относительно центров других элементов.
// За соседа принимаем те элементы, центры которых попали в радиус.
template<class Type, class Index>
void mesh_2d<Type, Index>::find_elements_neighbors(const Type r) {
    const std::vector<std::array<Type, 2>> centres = utils::get_centres_of_elements(*this);
    _elements_neighbors.resize(elements_count());
    for(size_t elL = 0; elL < elements_count(); ++elL) {
        _elements_neighbors[elL].resize(0);
        _elements_neighbors[elL].reserve(elements_count());
        for(size_t elNL = 0; elNL < elements_count(); ++elNL)
            if(utils::distance(centres[elL], centres[elNL]) < r)
                _elements_neighbors[elL].push_back(elNL);
        _elements_neighbors[elL].shrink_to_fit();
    }
}

// Ищет соседние элементы относительно узлов сетки.
// За соседа принимаем те элементы, центры которых попали в радиус.
template<class Type, class Index>
void mesh_2d<Type, Index>::find_nodes_neighbors(const Type r) {
    const std::vector<std::array<Type, 2>> centres = utils::get_centres_of_elements(*this);
    _nodes_neighbors.resize(nodes_count());
    for(size_t node = 0; node < nodes_count(); ++node) {
        _nodes_neighbors[node].resize(0);
        _nodes_neighbors[node].reserve(elements_count());
        for(size_t el = 0; el < elements_count(); ++el) {
            if(utils::distance(this->node(node), centres[el]) < r)
                _nodes_neighbors[node].push_back(el);
        }
        _nodes_neighbors[node].shrink_to_fit();
    }
}

}

#endif