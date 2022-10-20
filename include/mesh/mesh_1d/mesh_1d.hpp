#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include "metamath.hpp"

#include <cstdlib>
#include <memory>
#include <optional>
#include <numeric>
#include <ranges>

namespace nonlocal::mesh {

// Левый и правый элемент относительно узла, а также локальный номер относительно элемента
// struct left_right_element final {
//     struct node_element final {
//         size_t node, element;
//     };
//     std::optional<node_element> left, right;
// };

template<class T>
struct section_length_and_elements final {
    T length = T{1};
    size_t elements = 1;
};

template<class T>
class mesh_1d final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

    using finite_element_1d = metamath::finite_element::element_1d_integrate_base<T>;
    using finite_element_1d_ptr = std::unique_ptr<finite_element_1d>;

    finite_element_1d_ptr _element;
    std::vector<section_length_and_elements<T>> _sections;

    static std::vector<section_length_and_elements<T>> accumulate(const std::vector<section_length_and_elements<T>>& sections);

public:
    explicit mesh_1d(finite_element_1d_ptr&& element,
                     const std::vector<section_length_and_elements<T>>& sections);

    const finite_element_1d& element() const noexcept;

    size_t sections_count() const noexcept;
    size_t elements_count() const noexcept;
    size_t elements_count(const size_t section) const noexcept;
    size_t nodes_count() const noexcept;
    size_t nodes_count(const size_t section) const noexcept;

    T length() const noexcept;
    T length(const size_t section) const noexcept;
    T step(const size_t section) const noexcept;    
    T jacobian(const size_t section) const noexcept;
    std::array<T, 2> bounds(const size_t section) const noexcept;
};

template<class T>
mesh_1d<T>::mesh_1d(finite_element_1d_ptr&& element,
                    const std::vector<section_length_and_elements<T>>& sections)
    : _element{std::move(element)}
    , _sections{accumulate(sections)} {}

template<class T>
std::vector<section_length_and_elements<T>> mesh_1d<T>::accumulate(const std::vector<section_length_and_elements<T>>& sections) {
    std::vector<section_length_and_elements<T>> accumulated = sections;
    for(const size_t i : std::ranges::iota_view{size_t{1}, accumulated.size()}) {
        accumulated[i].length += accumulated[i - 1].length;
        accumulated[i].elements += accumulated[i - 1].elements;
    }
    return accumulated;
}

template<class T>
const mesh_1d<T>::finite_element_1d& mesh_1d<T>::element() const noexcept {
    return *_element;
}

template<class T>
size_t mesh_1d<T>::sections_count() const noexcept {
    return _sections.size();
}

template<class T>
size_t mesh_1d<T>::elements_count() const noexcept {
    return _sections.back().elements;
}

template<class T>
size_t mesh_1d<T>::elements_count(const size_t section) const noexcept {
    return section ? _sections[section].elements - _sections[section - 1].elements : _sections[section].elements;
}

template<class T>
size_t mesh_1d<T>::nodes_count() const noexcept {
    return elements_count() * (element().nodes_count() - 1) + 1;
}

template<class T>
size_t mesh_1d<T>::nodes_count(const size_t section) const noexcept {
    return elements_count(section) * (element().nodes_count() - 1) + 1;
}

template<class T>
T mesh_1d<T>::length() const noexcept {
    return _sections.back().length;
}

template<class T>
T mesh_1d<T>::length(const size_t section) const noexcept {
    return section ? _sections[section].length - _sections[section - 1].length : _sections.front().length;
}

template<class T>
T mesh_1d<T>::step(const size_t section) const noexcept {
    return length(section) / elements_count(section);
}

template<class T>
T mesh_1d<T>::jacobian(const size_t section) const noexcept {
    using enum metamath::finite_element::side_1d;
    return step(section) / (element().boundary(RIGHT) - element().boundary(LEFT));
}

template<class T>
std::array<T, 2> mesh_1d<T>::bounds(const size_t section) const noexcept {
    return {section ? _sections[section - 1].length : T{0}, _sections[section].length};
}




// // Класс сетки для одномерных задач.
// // Выбрана гипотеза, что сетка равномерная, состоит из однородных элементов и все элементы и узлы пронумерованы слева направо.
// // Все данные о нумерации и местоположении узлов вычисляются и не хранятся.
// template<class T>
// class mesh_1d final {
//     static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

//     using finite_element_1d = metamath::finite_element::element_1d_integrate_base<T>;
//     using finite_element_1d_ptr = std::unique_ptr<finite_element_1d>;

//     static finite_element_1d_ptr make_default_element();

//     finite_element_1d_ptr _element = make_default_element();
//     std::array<T, 2> _section = {T{-1}, T{1}};
//     size_t _elements_count = 1, _nodes_count = 2;
//     T _step = (_section.back() - _section.front()) / _elements_count;
//     T _jacobian = T{2} / _step;
//     std::vector<T> _quad_coord_loc;
//     size_t _neighbours_count;

// public:
//     explicit mesh_1d(finite_element_1d_ptr&& element, const size_t elements_count = 1, const std::array<T, 2> section = {T{-1}, T{1}});

//     const finite_element_1d& element() const noexcept;
//     const std::array<T, 2>& section() const noexcept;
//     size_t elements_count() const noexcept;
//     size_t nodes_count() const noexcept;
//     T jacobian() const noexcept;
//     T node_coord(const size_t node) const noexcept;
//     T quad_coord(const size_t e, const size_t q) const noexcept;
//     size_t node_number(const size_t e, const size_t i) const noexcept;
//     bool is_boundary_node(const size_t node) const noexcept;

//     curr_next_elements node_elements(const size_t node) const noexcept; // Возвращает номера элементов и локальные номера узлов на элементах

//     void calc_neighbours_count(const T r) noexcept;
//     size_t left_neighbour(const size_t e) const noexcept;
//     size_t right_neighbour(const size_t e) const noexcept;
// };

// template<class T>
// typename mesh_1d<T>::finite_element_1d_ptr mesh_1d<T>::make_default_element() {
//     using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, 1>;
//     using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, 1>;
//     return std::make_unique<element_1d>(quadrature{});
// }

// template<class T>
// mesh_1d<T>::mesh_1d(finite_element_1d_ptr&& element, const size_t elements_count, const std::array<T, 2> section)
//     : _element{std::move(element)}
//     , _elements_count{elements_count}
//     , _section{section}
//     , _nodes_count{_elements_count * (_element->nodes_count() - 1) + 1}
//     , _step{(_section.back() - _section.front()) / _elements_count}
//     , _jacobian{_step / (_element->boundary(metamath::finite_element::side_1d::RIGHT) - _element->boundary(metamath::finite_element::side_1d::LEFT))}
//     , _quad_coord_loc(_element->qnodes_count()) {
//     for(size_t q = 0; q < _quad_coord_loc.size(); ++q)
//         _quad_coord_loc[q] = (_element->quadrature().node(q) - _element->boundary(metamath::finite_element::side_1d::LEFT)) * _jacobian;
// }

// template<class T>
// const typename mesh_1d<T>::finite_element_1d& mesh_1d<T>::element() const noexcept { return *_element; }

// template<class T>
// const std::array<T, 2>& mesh_1d<T>::section() const noexcept { return _section; }

// template<class T>
// size_t mesh_1d<T>::elements_count() const noexcept { return _elements_count; }

// template<class T>
// size_t mesh_1d<T>::nodes_count() const noexcept { return _nodes_count; }

// template<class T>
// T mesh_1d<T>::jacobian() const noexcept { return _jacobian; }

// template<class T>
// T mesh_1d<T>::node_coord(const size_t node) const noexcept {
//     return section().front() + node * (section().back() - section().front()) / nodes_count();
// }

// template<class T>
// T mesh_1d<T>::quad_coord(const size_t e, const size_t q) const noexcept { return _step * e + _quad_coord_loc[q]; }

// template<class T>
// size_t mesh_1d<T>::node_number(const size_t e, const size_t i) const noexcept { return e * (element().nodes_count()-1) + i; }

// template<class T>
// bool mesh_1d<T>::is_boundary_node(const size_t node) const noexcept { return !node || node == nodes_count()-1; }

// template<class T>
// curr_next_elements mesh_1d<T>::node_elements(const size_t node) const noexcept {
//     if (!node)
//         return {.named = {.curr_element = 0, .curr_loc_number = 0}};
//     if (node == nodes_count()-1)
//         return {.named = {.curr_element = elements_count()-1, .curr_loc_number = _element->nodes_count()-1}};

//     const auto div = std::div(int64_t(node), int64_t(_element->nodes_count() - 1));
//     if (div.rem)
//         return {.named = {.curr_element = size_t(div.quot), .curr_loc_number = size_t(div.rem)}};
//     return {.named = {.curr_element = size_t(div.quot-1), .curr_loc_number = _element->nodes_count()-1,
//                       .next_element = size_t(div.quot),   .next_loc_number = 0}};
// }

// template<class T>
// void mesh_1d<T>::calc_neighbours_count(const T r) noexcept { _neighbours_count = r / _step + 1; }

// template<class T>
// size_t mesh_1d<T>::left_neighbour(const size_t e) const noexcept {
//     return e > _neighbours_count ? e - _neighbours_count : 0;
// }

// template<class T>
// size_t mesh_1d<T>::right_neighbour(const size_t e) const noexcept {
//     const size_t right = e + _neighbours_count;
//     return right > elements_count() ? elements_count() : right;
// }

}

#endif