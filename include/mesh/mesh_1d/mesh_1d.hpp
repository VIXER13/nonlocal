#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include "metamath.hpp"
#include <cstdlib>
#include <memory>

namespace mesh {

// Вспомогательная структура, необходимая для поузлового обхода сетки.
// Хранит в себе пару номеров элементов, которым принадлежит узел
// и локальные номера этого узла относительно каждого из элементов.
// Для удобства, имеет два варианта просмотра, именованный и в виде массива.
union curr_next_elements final {
    struct _curr_next_elements final {
        size_t curr_element    = std::numeric_limits<size_t>::max(),
               curr_loc_number = std::numeric_limits<size_t>::max(),
               next_element    = std::numeric_limits<size_t>::max(),
               next_loc_number = std::numeric_limits<size_t>::max();
    } named = {};
    std::array<std::array<size_t, 2>, 2> arr;
};

// Класс сетки для одномерных задач.
// Выбрана гипотеза, что сетка равномерная, состоит из однородных элементов и все элементы и узлы пронумерованы слева направо.
// Все данные о нумерации и местоположении узлов вычисляются и не хранятся.
template<class T, class I = int32_t>
class mesh_1d final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    static_assert(std::is_integral_v<I>, "The I must be integral.");

    using Finite_Element_1D_Ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

    static Finite_Element_1D_Ptr make_default_element();

    Finite_Element_1D_Ptr _element = make_default_element();
    std::array<T, 2> _section = {T{-1}, T{1}};
    size_t _elements_count = 1, _nodes_count = 2;
    T _step = (_section.back() - _section.front()) / _elements_count,
      _jacobian = T{2} / _step;
    std::vector<T> _quad_coord_loc;
    size_t _neighbours_count;

public:
    explicit mesh_1d(Finite_Element_1D_Ptr&& element, const size_t elements_count = 1, const std::array<T, 2> section = {T{-1}, T{1}});

    const Finite_Element_1D_Ptr& element() const noexcept;
    const std::array<T, 2>& section() const noexcept;
    size_t elements_count() const noexcept;
    size_t nodes_count() const noexcept;
    T jacobian() const noexcept;
    T quad_coord(const size_t e, const size_t q) const noexcept;
    size_t node_number(const size_t e, const size_t i) const noexcept;
    bool is_boundary_node(const size_t node) const noexcept;

    curr_next_elements node_elements(const size_t node) const noexcept; // Возвращает номера элементов и локальные номера узлов на элементах

    void calc_neighbours_count(const T r) noexcept;
    size_t left_neighbour(const size_t e) const noexcept;
    size_t right_neighbour(const size_t e) const noexcept;
};

template<class T, class I>
typename mesh_1d<T, I>::Finite_Element_1D_Ptr mesh_1d<T, I>::make_default_element() {
    using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss1>;
    using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::linear>;
    return std::make_unique<element_1d>(quadrature{});
}

template<class T, class I>
mesh_1d<T, I>::mesh_1d(Finite_Element_1D_Ptr&& element, const size_t elements_count, const std::array<T, 2> section)
    : _element{std::move(element)}
    , _elements_count{elements_count}
    , _section{section}
    , _nodes_count{_elements_count * (_element->nodes_count() - 1) + 1}
    , _step{(_section.back() - _section.front()) / _elements_count}
    , _jacobian{_step / (_element->boundary(metamath::finite_element::side_1d::RIGHT) - _element->boundary(metamath::finite_element::side_1d::LEFT))}
    , _quad_coord_loc(_element->qnodes_count()) {
    for(size_t q = 0; q < _quad_coord_loc.size(); ++q)
        _quad_coord_loc[q] = (_element->quadrature()->node(q)[0] - _element->boundary(metamath::finite_element::side_1d::LEFT)) * _jacobian;
}

template<class T, class I>
const typename mesh_1d<T, I>::Finite_Element_1D_Ptr& mesh_1d<T, I>::element() const noexcept { return _element; }

template<class T, class I>
const std::array<T, 2>& mesh_1d<T, I>::section() const noexcept { return _section; }

template<class T, class I>
size_t mesh_1d<T, I>::elements_count() const noexcept { return _elements_count; }

template<class T, class I>
size_t mesh_1d<T, I>::nodes_count() const noexcept { return _nodes_count; }

template<class T, class I>
T mesh_1d<T, I>::jacobian() const noexcept { return _jacobian; }

template<class T, class I>
T mesh_1d<T, I>::quad_coord(const size_t e, const size_t q) const noexcept { return _step * e + _quad_coord_loc[q]; }

template<class T, class I>
size_t mesh_1d<T, I>::node_number(const size_t e, const size_t i) const noexcept { return e * (element()->nodes_count()-1) + i; }

template<class T, class I>
bool mesh_1d<T, I>::is_boundary_node(const size_t node) const noexcept { return !node || node == nodes_count()-1; }

template<class T, class I>
curr_next_elements mesh_1d<T, I>::node_elements(const size_t node) const noexcept {
    if (!node)
        return {.named = {.curr_element = 0, .curr_loc_number = 0}};
    if (node == nodes_count()-1)
        return {.named = {.curr_element = elements_count()-1, .curr_loc_number = _element->nodes_count()-1}};

    const auto div = std::div(int64_t(node), int64_t(_element->nodes_count() - 1));
    if (div.rem)
        return {.named = {.curr_element = size_t(div.quot), .curr_loc_number = size_t(div.rem)}};
    return {.named = {.curr_element = size_t(div.quot-1), .curr_loc_number = _element->nodes_count()-1,
                      .next_element = size_t(div.quot),   .next_loc_number = 0}};
}

template<class T, class I>
void mesh_1d<T, I>::calc_neighbours_count(const T r) noexcept { _neighbours_count = r / _step + 1; }

template<class T, class I>
size_t mesh_1d<T, I>::left_neighbour(const size_t e) const noexcept {
    return e > _neighbours_count ? e - _neighbours_count : 0;
}

template<class T, class I>
size_t mesh_1d<T, I>::right_neighbour(const size_t e) const noexcept {
    const size_t right = e + _neighbours_count;
    return right > elements_count() ? elements_count() : right;
}

}

#endif