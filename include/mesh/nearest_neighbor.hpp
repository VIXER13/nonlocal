#ifndef MESH_NEAREST_NEIGHBOR_HPP
#define MESH_NEAREST_NEIGHBOR_HPP

#include "utils.hpp"
#include <iostream>

namespace mesh {

// Вычисление ориентированной площади треугольника на точках a, b, c.
template<class T>
T oriented_triangle_area(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c) noexcept {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

// Проверка отрезков AB и CD на предмет пересечения.
template<class T>
bool intersect(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c, const std::array<T, 2>& d) noexcept {
    static constexpr auto check_segment = [](T a, T b, T c, T d) {
        if (a > b) std::swap(a, b);
        if (c > d) std::swap(c, d);
        return std::max(a, c) <= std::min(b, d);
    };
    
    return check_segment(a[0], b[0], c[0], d[0]) && 
           check_segment(a[1], b[1], c[1], d[1]) &&
           oriented_triangle_area(a, b, c) * oriented_triangle_area(a, b, d) <= 0 &&
           oriented_triangle_area(c, d, a) * oriented_triangle_area(c, d, b) <= 0;
}

template<class Type, class Index>
bool mesh_2d<Type, Index>::check_intersect(const std::array<Type, 2>& A, const std::array<Type, 2>& B) const {
    for(size_t b = 0; b < boundary_groups_count(); ++b) {
        for(size_t el = 0; el < elements_count(b); ++el) {
            const auto& be = element_1d(element_1d_type(b, el));
            for(size_t i = 0; i < be->nodes_count()-1; ++i) {
                if(intersect(node(node_number(b, el, i)), node(node_number(b, el, i+1)), A, B))
                    return true;
            }
        }
    }
    return false;
}

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
            if(utils::distance(centres[elL], centres[elNL]) < r)// &&
               //!check_intersect(centres[elL], centres[elNL]))
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
            if(utils::distance(this->node(node), centres[el]) < r)// &&
               //!check_intersect(this->node(node), centres[el]))
                _nodes_neighbors[node].push_back(el);
        }
        _nodes_neighbors[node].shrink_to_fit();
    }
}

}

#endif