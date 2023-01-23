#ifndef NONLOCAL_MESH_SU2_PARSER_HPP
#define NONLOCAL_MESH_SU2_PARSER_HPP

namespace nonlocal::mesh {

template<class T, class I>
template<size_t... K, class Stream>
std::vector<I> mesh_container_2d<T, I>::read_element(Stream& mesh_file) {
    std::vector<I> element(sizeof...(K));
    (mesh_file >> ... >> element[K]);
    return element;
}

template<class T, class I>
template<class Stream>
std::vector<I> mesh_container_2d<T, I>::read_element(Stream& mesh_file, const size_t type) {
    switch(vtk_element_number(type)) {
        case vtk_element_number::LINEAR:
            return read_element<0, 1>(mesh_file);

        case vtk_element_number::QUADRATIC:
            return read_element<0, 2, 1>(mesh_file);

        case vtk_element_number::TRIANGLE:
            return read_element<0, 1, 2>(mesh_file);

        case vtk_element_number::QUADRATIC_TRIANGLE:
            return read_element<0, 1, 2, 3, 4, 5>(mesh_file);

        case vtk_element_number::BILINEAR:
            return read_element<0, 1, 2, 3>(mesh_file);

        case vtk_element_number::QUADRATIC_SERENDIPITY:
            return read_element<0, 2, 4, 6, 1, 3, 5, 7>(mesh_file);

        // TODO: fix parse quadratic lagrange elements
        //case vtk_element_number::QUADRATIC_LAGRANGE:
        //   return read_element<0, 2, 4, 6, 1, 3, 5, 7, 8>(mesh_file);

        default:
            throw std::domain_error{"Unknown element type."};
    }
}

template<class T, class I>
template<class Stream>
auto mesh_container_2d<T, I>::read_elements_2d(Stream& mesh_file) {
    std::string pass;
    size_t elements_count = 0;
    mesh_file >> pass >> pass >> pass >> elements_count;
    std::vector<std::vector<I>> elements_2d(elements_count);
    std::vector<uint8_t> elements_types_2d(elements_count);
    for(const size_t e : std::ranges::iota_view{0u, elements_count}) {
        size_t type = 0;
        mesh_file >> type;
        elements_types_2d[e] = type;
        elements_2d[e] = read_element(mesh_file, type);
        mesh_file >> pass;
    }
    return std::make_tuple(std::move(elements_2d), std::move(elements_types_2d));
}

template<class T, class I>
template<class Stream>
auto mesh_container_2d<T, I>::read_nodes(Stream& mesh_file) {
    size_t nodes_count = 0;
    std::string pass;
    mesh_file >> pass >> nodes_count;
    std::vector<std::array<T, 2>> nodes(nodes_count);
    for(auto& [x, y] : nodes)
        mesh_file >> x >> y >> pass;
    return nodes;
}

template<class T, class I>
template<class Stream>
auto mesh_container_2d<T, I>::read_elements_groups(Stream& mesh_file) {
    std::string pass;
    size_t groups_count = 0;
    mesh_file >> pass >> groups_count;
    std::unordered_set<std::string> groups_names;
    std::unordered_map<std::string, std::vector<std::vector<I>>> elements_in_groups(groups_count);
    std::unordered_map<std::string, std::vector<uint8_t>> types_in_groups(groups_count);
    for(const size_t group : std::ranges::iota_view{0u, groups_count}) {
        std::string group_name;
        size_t elements_count = 0;
        mesh_file >> pass >> group_name >> pass >> elements_count;

        auto& elements = elements_in_groups[group_name] = {};
        auto& types = types_in_groups[group_name] = {};
        groups_names.emplace(std::move(group_name));
        types.resize(elements_count);
        elements.resize(elements_count);

        for(const size_t e : std::ranges::iota_view{0u, elements_count}) {
            size_t type = 0;
            mesh_file >> type;
            types[e] = type;
            elements[e] = read_element(mesh_file, type);
        }
    }
    return std::make_tuple(std::move(groups_names), std::move(elements_in_groups), std::move(types_in_groups));
}

template<class T, class I>
template<class Stream>
void mesh_container_2d<T, I>::read_su2(Stream& mesh_file) {
    auto [elements_2d, elements_types_2d] = read_elements_2d(mesh_file);
    _elements_2d_count = elements_2d.size();
    _nodes = read_nodes(mesh_file);
    auto [groups_names, elements_in_groups, elements_types] = read_elements_groups(mesh_file);

    size_t elements_2d_shift = 0;
    size_t elements_shift = _elements_2d_count;
    for(const auto& [group, types] : elements_types)
        if (get_elements_set().is_element_1d(types.front())) {
            _elements_groups[group] = std::ranges::iota_view{elements_shift, elements_shift + types.size()};
            elements_shift += types.size();
            _groups_1d.insert(group);
        } else {
            _elements_groups[group] = std::ranges::iota_view{elements_2d_shift, elements_2d_shift + types.size()};
            elements_2d_shift += types.size();
            _groups_2d.insert(group);
        }
    if (elements_2d_shift > elements_2d.size())
        throw std::domain_error{"Problem with parsing groups: some groups of 2D elements overlap each other."};

    _elements.resize(elements_shift);
    _elements_types.resize(elements_shift);
    elements_2d_shift = 0;
    elements_shift = _elements_2d_count;
    for(const auto& [group, range] : _elements_groups) {
        auto& types = elements_types[group];
        auto& elements = elements_in_groups[group];
        if (get_elements_set().is_element_1d(types.front())) {
            for(const size_t e : std::ranges::iota_view{0u, range.size()}) {
                _elements[range.front() + e] = std::move(elements[e]);
                _elements_types[range.front() + e] = uint8_t(get_elements_set().model_to_local_1d(types[e]));
            }
            elements_shift += range.size();
        } else {
            for(const size_t e : std::ranges::iota_view{0u, range.size()}) {
                _elements[range.front() + e] = std::move(elements[e]);
                _elements_types[range.front() + e] = uint8_t(get_elements_set().model_to_local_2d(types[e]));
            }
            elements_2d_shift += range.size();
        }
    }

    if (elements_2d_shift < elements_2d.size()) {
        const auto default_range = std::ranges::iota_view{elements_2d_shift, elements_2d.size()};
        _elements_groups["Default"] = default_range;
        for(const size_t current_element : default_range)
            for(const size_t e : std::ranges::iota_view{0u, _elements_2d_count}) {
                if (elements_2d[e].empty())
                    continue;
                const auto can_inserted = [&element = elements_2d[e]](const std::vector<I>& el) { return el == element; };
                if (elements_2d_shift == 0 || std::any_of(_elements.begin(), std::next(_elements.begin(), elements_2d_shift), can_inserted)) {
                    _elements[current_element] = std::move(elements_2d[e]);
                    _elements_types[current_element] = uint8_t(get_elements_set().model_to_local_2d(elements_types_2d[e]));
                    break;
                }
            }
    }
}

}

#endif