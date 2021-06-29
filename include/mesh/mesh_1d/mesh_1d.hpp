#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include "metamath.hpp"
#include <cstdlib>
#include <memory>

namespace mesh {

union curr_next_elements final {
    struct _curr_next_elements final {
        size_t curr_element    = std::numeric_limits<size_t>::max(),
               curr_loc_number = std::numeric_limits<size_t>::max(),
               next_element    = std::numeric_limits<size_t>::max(),
               next_loc_number = std::numeric_limits<size_t>::max();
    } named = {};
    std::array<std::array<size_t, 2>, 2> arr;
};

template<class T, class I = int32_t>
class mesh_1d final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    static_assert(std::is_integral_v<I>, "The I must be integral.");

    using Finite_Element_1D_Ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

    static Finite_Element_1D_Ptr make_default_element();

    Finite_Element_1D_Ptr _element = make_default_element();
    std::array<T, 2> _section = {T{-1}, T{1}};
    size_t _elements_count = 1, _nodes_count = 2;

public:
    explicit mesh_1d() = default;
    explicit mesh_1d(Finite_Element_1D_Ptr&& element, const size_t elements_count = 1, const std::array<T, 2> section = {T{-1}, T{1}});

    const Finite_Element_1D_Ptr& element() const;
    const std::array<T, 2>& section() const;
    size_t elements_count() const;
    size_t nodes_count() const;

    // Возвращает номера элементов и локальные номера узлов на элементах
    curr_next_elements node_elements(const size_t node) const;

    //std::array<size_t, 2> node_elements(const size_t node) const;
    size_t node_begin(const size_t e) const;
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
    , _nodes_count{_elements_count * (_element->nodes_count() - 1) + 1} {}

template<class T, class I>
const typename mesh_1d<T, I>::Finite_Element_1D_Ptr& mesh_1d<T, I>::element() const { return _element; }

template<class T, class I>
const std::array<T, 2>& mesh_1d<T, I>::section() const { return _section; }

template<class T, class I>
size_t mesh_1d<T, I>::elements_count() const { return _elements_count; }

template<class T, class I>
size_t mesh_1d<T, I>::nodes_count() const { return _nodes_count; }

template<class T, class I>
curr_next_elements mesh_1d<T, I>::node_elements(const size_t node) const {
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
size_t mesh_1d<T, I>::node_begin(const size_t e) const { return e * (_element->nodes_count() - 1); }

}

#endif