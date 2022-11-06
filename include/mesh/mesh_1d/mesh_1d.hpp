#ifndef MESH_1D_HPP
#define MESH_1D_HPP

#include "metamath.hpp"

#include <cstdlib>
#include <memory>
#include <numeric>

namespace nonlocal::mesh {

struct left_right_element final {
    struct node_element final {
        metamath::types::optional_integer<size_t> element = std::nullopt;
        size_t node = 0;

        constexpr node_element() noexcept = default;
        constexpr node_element(const std::nullopt_t) noexcept {}
        constexpr node_element(const size_t e, const size_t i) noexcept
            : element{e}
            , node{i} {}

        node_element& operator=(const std::nullopt_t) noexcept {
            element = std::nullopt;
            node = 0;
            return *this;
        }

        constexpr explicit operator bool() const noexcept {
            return bool(element);
        }
    } left, right;

    constexpr std::array<node_element, 2> to_array() const noexcept {
        return {left, right};
    }
};

template<class T>
struct segment_data final {
    T length = T{1};
    size_t elements = 1;
};

template<class T>
class mesh_1d final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

    using finite_element_1d = metamath::finite_element::element_1d_integrate_base<T>;
    using finite_element_1d_ptr = std::unique_ptr<finite_element_1d>;

    finite_element_1d_ptr _element;
    std::vector<segment_data<T>> _segments;
    std::vector<size_t> _neighbours_count;

    static std::vector<segment_data<T>> accumulate(const std::vector<segment_data<T>>& segments);

public:
    explicit mesh_1d(finite_element_1d_ptr&& element, const std::vector<segment_data<T>>& segments);

    const finite_element_1d& element() const noexcept;

    size_t segments_count() const noexcept;
    size_t elements_count() const noexcept;
    size_t elements_count(const size_t segment) const noexcept;
    size_t nodes_count() const noexcept;
    size_t nodes_count(const size_t segment) const noexcept;
    size_t segment_number(const size_t e) const noexcept;
    size_t node_number(const size_t e, const size_t i) const noexcept;
    std::ranges::iota_view<size_t, size_t> segment_elements(const size_t segment) const noexcept;
    std::ranges::iota_view<size_t, size_t> segment_nodes(const size_t segment) const noexcept;
    std::ranges::iota_view<size_t, size_t> neighbours(const size_t e) const noexcept;
    left_right_element node_elements(const size_t node) const noexcept;

    T length() const noexcept;
    T length(const size_t segment) const noexcept;
    T step(const size_t segment) const noexcept;    
    T jacobian(const size_t segment) const noexcept;
    std::array<T, 2> bounds(const size_t segment) const noexcept;

    void find_neighbours(const std::vector<T>& radii);
};

template<class T>
mesh_1d<T>::mesh_1d(finite_element_1d_ptr&& element, const std::vector<segment_data<T>>& segments)
    : _element{std::move(element)}
    , _segments{accumulate(segments)}
    , _neighbours_count(_segments.size(), 1) {}

template<class T>
std::vector<segment_data<T>> mesh_1d<T>::accumulate(const std::vector<segment_data<T>>& segments) {
    std::vector<segment_data<T>> accumulated = segments;
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
size_t mesh_1d<T>::segments_count() const noexcept {
    return _segments.size();
}

template<class T>
size_t mesh_1d<T>::elements_count() const noexcept {
    return _segments.back().elements;
}

template<class T>
size_t mesh_1d<T>::elements_count(const size_t segment) const noexcept {
    return segment ? _segments[segment].elements - _segments[segment - 1].elements : _segments[segment].elements;
}

template<class T>
size_t mesh_1d<T>::nodes_count() const noexcept {
    return elements_count() * (element().nodes_count() - 1) + 1;
}

template<class T>
size_t mesh_1d<T>::nodes_count(const size_t segment) const noexcept {
    return elements_count(segment) * (element().nodes_count() - 1) + 1;
}

template<class T>
size_t mesh_1d<T>::segment_number(const size_t e) const noexcept {
    size_t segment = 0;
    for(; _segments[segment + 1].elements < e; ++segment);
    return segment;
}

template<class T>
size_t mesh_1d<T>::node_number(const size_t e, const size_t i) const noexcept {
    return e * (element().nodes_count() - 1) + i;
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::segment_elements(const size_t segment) const noexcept {
    return {segment ? _segments[segment - 1].elements : 0, _segments[segment].elements};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::segment_nodes(const size_t segment) const noexcept {
    const size_t elements_nodes_count = element().nodes_count() - 1;
    const std::ranges::iota_view<size_t, size_t> elements = segment_elements(segment);
    return {elements_nodes_count * elements.front(), elements_nodes_count * (elements.back() + 1) + 1};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::neighbours(const size_t e) const noexcept {
    const size_t segment = segment_number(e);
    const size_t left_bound  = segment ? _segments[segment - 1].elements : 0;
    const size_t right_bound = _segments[segment].elements;
    const size_t left  = e > _neighbours_count[segment] ? e - _neighbours_count[segment] : 0;
    const size_t right = e + _neighbours_count[segment];
    return {
        left  < left_bound  ? left_bound  : left,
        right > right_bound ? right_bound : right
    };
}

template<class T>
left_right_element mesh_1d<T>::node_elements(const size_t node) const noexcept {
    using data = left_right_element::node_element;

    if (!node)
        return {.left  = data{0, 0},
                .right = std::nullopt};

    const size_t elements_nodes_count = element().nodes_count() - 1;
    if (node == nodes_count() - 1)
        return {.left  = std::nullopt, 
                .right = data{elements_count() - 1, elements_nodes_count}};

    const auto [e, i] = std::div(int64_t(node), int64_t(elements_nodes_count));
    if (i)
        return {.left  = data(e, i),
                .right = std::nullopt};
    return {.left  = data(e - 1, elements_nodes_count),
            .right = data(e,     0)};
}

template<class T>
T mesh_1d<T>::length() const noexcept {
    return _segments.back().length;
}

template<class T>
T mesh_1d<T>::length(const size_t segment) const noexcept {
    return segment ? _segments[segment].length - _segments[segment - 1].length : _segments.front().length;
}

template<class T>
T mesh_1d<T>::step(const size_t segment) const noexcept {
    return length(segment) / elements_count(segment);
}

template<class T>
T mesh_1d<T>::jacobian(const size_t segment) const noexcept {
    using enum metamath::finite_element::side_1d;
    return step(segment) / (element().boundary(RIGHT) - element().boundary(LEFT));
}

template<class T>
std::array<T, 2> mesh_1d<T>::bounds(const size_t segment) const noexcept {
    return {segment ? _segments[segment - 1].length : T{0}, _segments[segment].length};
}

template<class T>
void mesh_1d<T>::find_neighbours(const std::vector<T>& radii) {
    if (radii.size() != _neighbours_count.size())
        throw std::runtime_error{"The number of radii does not match the number of segments."};
    for(const size_t segment : std::ranges::iota_view{size_t{0}, segments_count()})
        _neighbours_count[segment] = radii[segment] / step(segment) + 1;
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
//     std::array<T, 2> _segment = {T{-1}, T{1}};
//     size_t _elements_count = 1, _nodes_count = 2;
//     T _step = (_segment.back() - _segment.front()) / _elements_count;
//     T _jacobian = T{2} / _step;
//     std::vector<T> _quad_coord_loc;
//     size_t _neighbours_count;

// public:
//     explicit mesh_1d(finite_element_1d_ptr&& element, const size_t elements_count = 1, const std::array<T, 2> segment = {T{-1}, T{1}});

//     const finite_element_1d& element() const noexcept;
//     const std::array<T, 2>& segment() const noexcept;
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
// mesh_1d<T>::mesh_1d(finite_element_1d_ptr&& element, const size_t elements_count, const std::array<T, 2> segment)
//     : _element{std::move(element)}
//     , _elements_count{elements_count}
//     , _segment{segment}
//     , _nodes_count{_elements_count * (_element->nodes_count() - 1) + 1}
//     , _step{(_segment.back() - _segment.front()) / _elements_count}
//     , _jacobian{_step / (_element->boundary(metamath::finite_element::side_1d::RIGHT) - _element->boundary(metamath::finite_element::side_1d::LEFT))}
//     , _quad_coord_loc(_element->qnodes_count()) {
//     for(size_t q = 0; q < _quad_coord_loc.size(); ++q)
//         _quad_coord_loc[q] = (_element->quadrature().node(q) - _element->boundary(metamath::finite_element::side_1d::LEFT)) * _jacobian;
// }

// template<class T>
// const typename mesh_1d<T>::finite_element_1d& mesh_1d<T>::element() const noexcept { return *_element; }

// template<class T>
// const std::array<T, 2>& mesh_1d<T>::segment() const noexcept { return _segment; }

// template<class T>
// size_t mesh_1d<T>::elements_count() const noexcept { return _elements_count; }

// template<class T>
// size_t mesh_1d<T>::nodes_count() const noexcept { return _nodes_count; }

// template<class T>
// T mesh_1d<T>::jacobian() const noexcept { return _jacobian; }

// template<class T>
// T mesh_1d<T>::node_coord(const size_t node) const noexcept {
//     return segment().front() + node * (segment().back() - segment().front()) / nodes_count();
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