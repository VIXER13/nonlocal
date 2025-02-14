#pragma once

#include "mesh_container_2d.hpp"
#include "vtk_elements_set.hpp"

namespace nonlocal::mesh {

template<class T, class I>
class mesh_parser<T, I, mesh_format::SU2> final {
    mesh_container_2d<T, I>& _mesh;

    template<size_t... K>
    static std::vector<I> read_element(std::stringstream& stream);
    static std::tuple<size_t, std::vector<I>> read_element(const std::string& element);

    template<class Stream>
    static std::unordered_set<std::string> read_elements(Stream& mesh_file, const size_t count, const bool is_group);
    template<class Stream>
    static std::unordered_set<std::string> read_elements_2d(Stream& mesh_file);
    template<class Stream>
    static std::vector<std::array<T, 2>> read_nodes(Stream& mesh_file);
    template<class Stream>
    static std::unordered_map<std::string, std::unordered_set<std::string>> read_elements_groups(Stream& mesh_file);

    static void filter_elements(std::unordered_set<std::string>& elements_2d,
                                const std::unordered_map<std::string, std::unordered_set<std::string>>& elements_groups);

    bool is_1d_group(const std::unordered_set<std::string>& element_group) const;
    size_t read_group(const std::string& group, const std::unordered_set<std::string>& elements, size_t element_shift);

public:
    explicit mesh_parser(mesh_container_2d<T, I>& mesh) noexcept;

    template<class Stream>
    void parse(Stream& mesh_file);
};

template<class T, class I>
mesh_parser<T, I, mesh_format::SU2>::mesh_parser(mesh_container_2d<T, I>& mesh) noexcept
    : _mesh{mesh} {}

template<class T, class I>
template<size_t... K>
std::vector<I> mesh_parser<T, I, mesh_format::SU2>::read_element(std::stringstream& stream) {
    std::vector<I> element(sizeof...(K));
    (stream >> ... >> element[K]);
    return element;
}

template<class T, class I>
std::tuple<size_t, std::vector<I>> mesh_parser<T, I, mesh_format::SU2>::read_element(const std::string& element) {
    size_t type = 0;
    std::stringstream stream{element};
    stream >> type;
    switch(vtk_element_number(type)) {
        case vtk_element_number::LINEAR:
            return {type, read_element<0, 1>(stream)};
        case vtk_element_number::QUADRATIC:
            return {type, read_element<0, 2, 1>(stream)};
        case vtk_element_number::TRIANGLE:
            return {type, read_element<0, 1, 2>(stream)};
        case vtk_element_number::QUADRATIC_TRIANGLE:
            return {type, read_element<0, 1, 2, 3, 4, 5>(stream)};
        case vtk_element_number::BILINEAR:
            return {type, read_element<0, 1, 2, 3>(stream)};
        case vtk_element_number::QUADRATIC_SERENDIPITY:
            return {type, read_element<0, 2, 4, 6, 1, 3, 5, 7>(stream)};
        case vtk_element_number::QUADRATIC_LAGRANGE:
           return {type, read_element<2, 0, 6, 8, 1, 3, 7, 5, 4>(stream)};
        default:
            throw std::domain_error{"Unknown element type."};
    }
}

template<class T, class I>
template<class Stream>
std::unordered_set<std::string> mesh_parser<T, I, mesh_format::SU2>::read_elements(Stream& mesh_file, const size_t count, const bool is_group) {
    std::string element;
    std::unordered_set<std::string> elements(count);
    for(const size_t e : std::ranges::iota_view{0u, count}) {
        std::getline(mesh_file, element);
        element.resize(element.rfind(' ')); // Remove element number from string
        if (element.back() == ' ')
            element.pop_back(); // remove final space if there is one
        elements.emplace(std::move(element));
    }
    return elements;
}


template<class T, class I>
template<class Stream>
std::unordered_set<std::string> mesh_parser<T, I, mesh_format::SU2>::read_elements_2d(Stream& mesh_file) {
    std::string _;
    size_t elements_count = 0;
    mesh_file >> _ >> _ >> _ >> elements_count;
    std::getline(mesh_file, _); // read the line to the end
    static constexpr bool Is_Group = false;
    return read_elements(mesh_file, elements_count, Is_Group);
}

template<class T, class I>
template<class Stream>
std::vector<std::array<T, 2>> mesh_parser<T, I, mesh_format::SU2>::read_nodes(Stream& mesh_file) {
    size_t nodes_count = 0;
    std::string _;
    mesh_file >> _ >> nodes_count;
    std::vector<std::array<T, 2>> nodes(nodes_count);
    for(auto& [x, y] : nodes)
        mesh_file >> x >> y >> _;
    return nodes;
}

template<class T, class I>
template<class Stream>
std::unordered_map<std::string, std::unordered_set<std::string>> 
mesh_parser<T, I, mesh_format::SU2>::read_elements_groups(Stream& mesh_file) {
    std::string _;
    size_t groups_count = 0;
    mesh_file >> _ >> groups_count;
    std::unordered_map<std::string, std::unordered_set<std::string>> groups;
    for(const size_t group : std::ranges::iota_view{0u, groups_count}) {
        std::string group_name;
        size_t elements_count = 0;
        mesh_file >> _ >> group_name >> _ >> elements_count;
        std::getline(mesh_file, _); // read the line to the end
        static constexpr bool Is_Group = true;
        groups[group_name] = read_elements(mesh_file, elements_count, Is_Group);
    }
    return groups;
}

template<class T, class I>
void mesh_parser<T, I, mesh_format::SU2>::filter_elements(
    std::unordered_set<std::string>& elements_2d,
    const std::unordered_map<std::string, std::unordered_set<std::string>>& elements_groups) {
    for(const auto& [group, elements] : elements_groups)
        for(const auto& element : elements) 
            if (const auto it = elements_2d.find(element); it != elements_2d.end()) 
                elements_2d.erase(it);
}

template<class T, class I>
bool mesh_parser<T, I, mesh_format::SU2>::is_1d_group(const std::unordered_set<std::string>& element_group) const {
    // We assume that the group contains only one dimension elements
    const auto [type, _] = read_element(*element_group.begin()); // We do a test reading of the first element
    return _mesh.get_elements_set().is_element_1d(type);
}

template<class T, class I>
size_t mesh_parser<T, I, mesh_format::SU2>::read_group(const std::string& group, const std::unordered_set<std::string>& elements, size_t element_shift) {
    const auto& elements_set = _mesh.get_elements_set();
    _mesh._elements_groups[group] = std::ranges::iota_view{element_shift, element_shift + elements.size()};
    for(const auto& element : elements) {
        std::tie(_mesh._elements_types[element_shift], _mesh._elements[element_shift]) = read_element(element);
        _mesh._elements_types[element_shift] = elements_set.is_element_1d(_mesh._elements_types[element_shift]) ?
                                               uint8_t(elements_set.model_to_local_1d(_mesh._elements_types[element_shift])) :
                                               uint8_t(elements_set.model_to_local_2d(_mesh._elements_types[element_shift]));
        ++element_shift;
    }
    return element_shift;
}

template<class T, class I>
template<class Stream>
void mesh_parser<T, I, mesh_format::SU2>::parse(Stream& mesh_file) {
    _mesh._elements_set = std::make_unique<vtk_elements_set<T>>();
    auto elements_2d = read_elements_2d(mesh_file);
    _mesh._elements_2d_count = elements_2d.size();
    _mesh._nodes = read_nodes(mesh_file);
    auto elements_groups = read_elements_groups(mesh_file);
    if (elements_groups.contains(Default_Group_Name))
        throw std::domain_error{"The group name cannot be named \"" + Default_Group_Name + "\", as it is a reserved group name."};
    filter_elements(elements_2d, elements_groups);
    if (!elements_2d.empty())
        elements_groups[Default_Group_Name] = std::move(elements_2d);

    constexpr auto elements_counter = [](const size_t sum, const auto& group) noexcept { return sum + group.second.size(); };
    _mesh._elements.resize(std::accumulate(elements_groups.begin(), elements_groups.end(), size_t{0}, elements_counter));
    _mesh._elements_types.resize(_mesh._elements.size());

    size_t elements_2d_shift = 0;
    size_t elements_1d_shift = _mesh._elements_2d_count;
    for(const auto& [group, elements] : elements_groups) {
        if (is_1d_group(elements)) {
            _mesh._groups_1d.insert(group);
            elements_1d_shift = read_group(group, elements, elements_1d_shift);
        } else {
            _mesh._groups_2d.insert(group);
            elements_2d_shift = read_group(group, elements, elements_2d_shift);
        }
    }
}

}