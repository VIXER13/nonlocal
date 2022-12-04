#ifndef NONLOCAL_MESH_CONTAINER_2D_HPP
#define NONLOCAL_MESH_CONTAINER_2D_HPP

#include "vtk_elements_set.hpp"

#include <filesystem>
#include <fstream>
#include <ranges>
#include <numeric>

namespace nonlocal::mesh {

template<class T, class I>
class mesh_container_2d final {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
    static_assert(std::is_integral_v<I>, "The I must be integral.");

    std::unique_ptr<elements_set<T>> _elements_set; // TODO: make elements_set copyable
    std::vector<std::array<T, 2>> _nodes;
    std::vector<std::vector<I>> _elements;
    std::vector<uint8_t> _elements_types;
    std::vector<std::string> _groups_names_1d;
    std::vector<std::string> _groups_names_2d;
    std::unordered_map<std::string, std::ranges::iota_view<size_t, size_t>> _elements_groups;
    size_t _elements_2d_count = 0u;

    // TODO: move parser to other class or function
    template<size_t... K, class Stream>
    std::vector<I> read_element(Stream& mesh_file);
    template<class Stream>
    auto read_elements_2d(Stream& mesh_file);
    template<class Stream>
    auto read_nodes(Stream& mesh_file);
    template<class Stream>
    auto read_elements_1d(Stream& mesh_file);
    template<class Stream>
    void read_su2(Stream& mesh_file);

public:
    struct element_data_2d final {
        const mesh_container_2d& mesh;
        const std::vector<I>& nodes;
        const element_integrate_2d<T>& element;
        
        std::array<T, 2> quad_coord(const size_t q) const;
        metamath::types::square_matrix<T, 2> jacobi_matrix(const size_t q) const;
    };

    explicit mesh_container_2d(const std::filesystem::path& path_to_mesh);

    const std::vector<std::string>& groups_names_1d() const noexcept;
    const std::vector<std::string>& groups_names_2d() const noexcept;
    size_t groups_1d_count() const noexcept;
    size_t groups_2d_count() const noexcept;
    
    size_t elements_count() const noexcept;
    size_t elements_count(const std::string& group_name) const;
    size_t elements_1d_count() const;
    size_t elements_2d_count() const;

    std::ranges::iota_view<size_t, size_t> elements() const noexcept;
    std::ranges::iota_view<size_t, size_t> elements(const std::string& group_name) const;
    std::ranges::iota_view<size_t, size_t> elements_1d() const noexcept;
    std::ranges::iota_view<size_t, size_t> elements_2d() const noexcept;

    size_t nodes_count() const noexcept;
    size_t nodes_count(const size_t element) const;
    size_t node_number(const size_t element, const size_t i) const;

    const std::vector<I>& nodes(const size_t element) const;
    const std::array<T, 2>& node_coord(const size_t node) const;

    const elements_set<T>& get_elements_set() const;

    const element_integrate_1d<T>& element_1d(const size_t element) const;
    const element_integrate_2d<T>& element_2d(const size_t element) const;

    element_data_2d element_2d_data(const size_t element) const;

    void clear();
    void read_from_file(const std::filesystem::path& path_to_mesh);
    void renumbering(const std::vector<size_t>& permutation);
};

template<class T, class I>
std::array<T, 2> mesh_container_2d<T, I>::element_data_2d::quad_coord(const size_t q) const {
    std::array<T, 2> coord = {};
    using namespace metamath::functions;
    for(const size_t i : std::ranges::iota_view{0u, element.nodes_count()})
        coord += mesh.node_coord(nodes[i]) * element.qN(i, q);
    return coord;
}

template<class T, class I>
metamath::types::square_matrix<T, 2> mesh_container_2d<T, I>::element_data_2d::jacobi_matrix(const size_t q) const {
    metamath::types::square_matrix<T, 2> J = {};
    for(const size_t i : std::ranges::iota_view{0u, element.nodes_count()}) {
        const std::array<T, 2> derivative = {element.qNxi (i, q), element.qNeta(i, q)};
        using namespace metamath::functions;
        J[0] += mesh.node_coord(nodes[i])[0] * derivative;
        J[1] += mesh.node_coord(nodes[i])[1] * derivative;
    }
    return J;
}

template<class T, class I>
mesh_container_2d<T, I>::mesh_container_2d(const std::filesystem::path& path_to_mesh) {
    read_from_file(path_to_mesh);
}

template<class T, class I>
const std::vector<std::string>& mesh_container_2d<T, I>::groups_names_1d() const noexcept {
    return _groups_names_1d;
}

template<class T, class I>
const std::vector<std::string>& mesh_container_2d<T, I>::groups_names_2d() const noexcept {
    return _groups_names_2d;
}

template<class T, class I>
size_t mesh_container_2d<T, I>::groups_1d_count() const noexcept {
    return groups_names_1d().size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::groups_2d_count() const noexcept {
    return groups_names_2d().size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::elements_count() const noexcept {
    return _elements.size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::elements_count(const std::string& group_name) const {
    return _elements_groups.at(group_name).size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::elements_1d_count() const {
    return elements_count() - elements_2d_count();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::elements_2d_count() const {
    return _elements_2d_count;
}

template<class T, class I>
std::ranges::iota_view<size_t, size_t> mesh_container_2d<T, I>::elements() const noexcept {
    return {0u, elements_count()};
}

template<class T, class I>
std::ranges::iota_view<size_t, size_t> mesh_container_2d<T, I>::elements(const std::string& group_name) const {
    return _elements_groups.at(group_name);
}

template<class T, class I>
std::ranges::iota_view<size_t, size_t> mesh_container_2d<T, I>::elements_1d() const noexcept {
    return {elements_2d_count(), elements_count()};
}

template<class T, class I>
std::ranges::iota_view<size_t, size_t> mesh_container_2d<T, I>::elements_2d() const noexcept {
    return {0u, elements_2d_count()};
}

template<class T, class I>
size_t mesh_container_2d<T, I>::nodes_count() const noexcept {
    return _nodes.size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::nodes_count(const size_t element) const {
    return nodes(element).size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::node_number(const size_t element, const size_t i) const {
    return nodes(element)[i];
}

template<class T, class I>
const std::vector<I>& mesh_container_2d<T, I>::nodes(const size_t element) const {
    return _elements[element];
}

template<class T, class I>
const std::array<T, 2>& mesh_container_2d<T, I>::node_coord(const size_t node) const {
    return _nodes[node];
}

template<class T, class I>
const elements_set<T>& mesh_container_2d<T, I>::get_elements_set() const {
    return *_elements_set;
}

template<class T, class I>
const element_integrate_1d<T>& mesh_container_2d<T, I>::element_1d(const size_t element) const {
    return get_elements_set().element_1d(_elements_types[element]);
}

template<class T, class I>
const element_integrate_2d<T>& mesh_container_2d<T, I>::element_2d(const size_t element) const {
    return get_elements_set().element_2d(_elements_types[element]);
}

template<class T, class I>
mesh_container_2d<T, I>::element_data_2d mesh_container_2d<T, I>::element_2d_data(const size_t element) const {
    return {.mesh = *this, .nodes = nodes(element), .element = element_2d(element)};
}

template<class T, class I>
void mesh_container_2d<T, I>::clear() {
    _elements_set = nullptr;
    _nodes.clear();
    _nodes.shrink_to_fit();
    _elements.clear();
    _elements.shrink_to_fit();
    _elements_types.clear();
    _elements_types.shrink_to_fit();
    _groups_names_1d.clear();
    _groups_names_1d.shrink_to_fit();
    _groups_names_2d.clear();
    _groups_names_2d.shrink_to_fit();
    _elements_groups.clear();
    _elements_2d_count = 0u;
}

template<class T, class I>
void mesh_container_2d<T, I>::read_from_file(const std::filesystem::path& path_to_mesh) {
    const std::string extension = path_to_mesh.extension().string();
    if (extension == ".su2") {
        clear();
        _elements_set = std::make_unique<vtk_elements_set<T>>();
        std::ifstream mesh_file{path_to_mesh};
        read_su2(mesh_file);
        return;
    }
    throw std::domain_error{"Unable to read mesh with extension " + extension};
}

template<class T, class I>
void mesh_container_2d<T, I>::renumbering(const std::vector<size_t>& permutation) {
    if (permutation.size() != nodes_count())
        throw std::runtime_error{"Permutation size does not match the mesh nodes number"};

    std::vector<std::array<T, 2>> nodes(_nodes.size());
    for(const size_t i : std::ranges::iota_view{0u, nodes_count()})
        nodes[permutation[i]] = _nodes[i];
    _nodes = std::move(nodes);

    for(std::vector<I>& nodes : _elements)
        for(I& node : nodes)
            node = permutation[node];
}

}

#include "su2_parser.hpp"

#endif