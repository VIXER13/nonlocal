#ifndef NONLOCAL_MESH_1D_HPP
#define NONLOCAL_MESH_1D_HPP

#include "metamath.hpp"

#include <cstdlib>
#include <memory>
#include <numeric>

#include <iostream>

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

    constexpr size_t count() const noexcept {
        return bool(left) + bool(right);
    }

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
    std::ranges::iota_view<size_t, size_t> segments() const noexcept;
    std::ranges::iota_view<size_t, size_t> elements() const noexcept;
    std::ranges::iota_view<size_t, size_t> nodes() const noexcept;
    std::ranges::iota_view<size_t, size_t> elements(const size_t segment) const noexcept;
    std::ranges::iota_view<size_t, size_t> nodes(const size_t segment) const noexcept;
    std::ranges::iota_view<size_t, size_t> neighbours(const size_t e) const noexcept;
    left_right_element node_elements(const size_t node) const noexcept;

    T length() const noexcept;
    T length(const size_t segment) const noexcept;
    T element_length(const size_t segment) const noexcept;
    T jacobian(const size_t segment) const noexcept;
    T node_coord(const size_t node) const noexcept;
    T qnode_coord(const size_t e, const size_t q) const noexcept;
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
    for(const size_t i : std::ranges::iota_view{1u, accumulated.size()}) {
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
    while(_segments[segment].elements <= e)
        ++segment;
    return segment;
}

template<class T>
size_t mesh_1d<T>::node_number(const size_t e, const size_t i) const noexcept {
    return e * (element().nodes_count() - 1) + i;
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::segments() const noexcept {
    return {size_t{0}, segments_count()};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::elements() const noexcept {
    return {size_t{0}, elements_count()};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::nodes() const noexcept {
    return {size_t{0}, nodes_count()};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::elements(const size_t segment) const noexcept {
    return {segment ? _segments[segment - 1].elements : 0, _segments[segment].elements};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::nodes(const size_t segment) const noexcept {
    const size_t elements_nodes_count = element().nodes_count() - 1;
    const std::ranges::iota_view<size_t, size_t> elements = mesh_1d<T>::elements(segment);
    return {elements_nodes_count * elements.front(), elements_nodes_count * (elements.back() + 1) + 1};
}

template<class T>
std::ranges::iota_view<size_t, size_t> mesh_1d<T>::neighbours(const size_t e) const noexcept {
    const size_t segment = segment_number(e);
    const size_t left_bound  = segment ? _segments[segment - 1].elements : 0;
    const size_t right_bound = _segments[segment].elements;
    const size_t left  = e > _neighbours_count[segment] ? e - _neighbours_count[segment] : 0;
    const size_t right = e + _neighbours_count[segment] + 1;
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
T mesh_1d<T>::element_length(const size_t segment) const noexcept {
    return length(segment) / elements_count(segment);
}

template<class T>
T mesh_1d<T>::jacobian(const size_t segment) const noexcept {
    using enum metamath::finite_element::side_1d;
    return element_length(segment) / (element().boundary(RIGHT) - element().boundary(LEFT));
}

template<class T>
T mesh_1d<T>::node_coord(const size_t node) const noexcept {
    size_t segment = 0, accumulated_nodes = nodes_count(segment) - 1;
    while(node >= accumulated_nodes && segment < segments_count() - 1)
        accumulated_nodes += nodes_count(++segment) - 1;
    const T left_bound = segment ? _segments[segment - 1].length : T{0};
    const T distance_between_nodes = length(segment) / (nodes_count(segment) - 1);
    const size_t node_number_in_segment = segment ? node + nodes_count(segment) - accumulated_nodes - 1 : node;
    return left_bound + node_number_in_segment * distance_between_nodes;
}

template<class T>
T mesh_1d<T>::qnode_coord(const size_t e, const size_t q) const noexcept {
    using enum metamath::finite_element::side_1d;
    const size_t segment = segment_number(e);
    const auto [left_bound, _] = bounds(segment);
    const size_t first_segment_element = elements(segment).front();
    const T qnode_coord_loc = (element().quadrature().node(q) - element().boundary(LEFT)) * jacobian(segment);
    return left_bound + element_length(segment) * (e - first_segment_element) + qnode_coord_loc;
}

template<class T>
std::array<T, 2> mesh_1d<T>::bounds(const size_t segment) const noexcept {
    return {segment ? _segments[segment - 1].length : T{0}, _segments[segment].length};
}

template<class T>
void mesh_1d<T>::find_neighbours(const std::vector<T>& radii) {
    if (radii.size() != _neighbours_count.size())
        throw std::runtime_error{"The number of radii does not match the number of segments."};
    for(const size_t segment : std::ranges::iota_view{0u, segments_count()})
        _neighbours_count[segment] = std::round(radii[segment] / element_length(segment));
}

}

#endif