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

    std::unique_ptr<elements_set<T>> _elements_set = std::make_unique<vtk_elements_set<T>>(); // TODO: remove static

    std::vector<std::array<T, 2>> _nodes;
    std::vector<std::vector<I>> _elements;
    std::vector<uint8_t> _elements_types;
    std::vector<std::string> _groups_names_1d;
    std::vector<std::string> _groups_names_2d;
    std::unordered_map<std::string, std::ranges::iota_view<size_t, size_t>> _elements_groups;
    size_t _elements_2d_count = 0u;

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

    size_t nodes_count() const noexcept;
    size_t nodes_count(const size_t e) const;
    size_t node_number(const size_t e, const size_t i) const;

    const std::array<T, 2>& node_coord(const size_t node) const;

    const elements_set<T>& get_elements_set() const;

    const element_integrate_1d<T>& element_1d(const size_t e) const;
    const element_integrate_2d<T>& element_2d(const size_t e) const;

    void clear();
    void read_from_file(const std::filesystem::path& path_to_mesh);
};

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
size_t mesh_container_2d<T, I>::nodes_count() const noexcept {
    return _nodes.size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::nodes_count(const size_t e) const {
    return _elements[e].size();
}

template<class T, class I>
size_t mesh_container_2d<T, I>::node_number(const size_t e, const size_t i) const {
    return _elements[e][i];
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
const element_integrate_1d<T>& mesh_container_2d<T, I>::element_1d(const size_t e) const {
    return get_elements_set().element_1d(_elements_types[e]);
}

template<class T, class I>
const element_integrate_2d<T>& mesh_container_2d<T, I>::element_2d(const size_t e) const {
    return get_elements_set().element_2d(_elements_types[e]);
}

template<class T, class I>
void mesh_container_2d<T, I>::clear() {
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
        std::ifstream mesh_file{path_to_mesh};
        read_su2(mesh_file);
        return;
    }
    throw std::domain_error{"Unable to read mesh with extension " + extension};
}

}

#endif