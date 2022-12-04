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
auto mesh_container_2d<T, I>::read_elements_2d(Stream& mesh_file) {
    size_t elements_count = 0;
    std::string pass;
    mesh_file >> pass >> pass >> pass >> elements_count;
    std::vector<std::vector<I>> elements_2d(elements_count);
    std::vector<uint8_t> elements_types_2d(elements_count);
    for(const size_t e : std::ranges::iota_view{0u, elements_count}) {
        size_t type = 0;
        mesh_file >> type;
        elements_types_2d[e] = get_elements_set().model_to_local_2d(type);
        switch(vtk_element_number(type)) {
            case vtk_element_number::TRIANGLE:
                elements_2d[e] = read_element<0, 1, 2>(mesh_file);
            break;

            case vtk_element_number::QUADRATIC_TRIANGLE:
               elements_2d[e] = read_element<0, 1, 2, 3, 4, 5>(mesh_file);
            break;

            case vtk_element_number::BILINEAR:
                elements_2d[e] = read_element<0, 1, 2, 3>(mesh_file);
            break;

            case vtk_element_number::QUADRATIC_SERENDIPITY:
                elements_2d[e] = read_element<0, 2, 4, 6, 1, 3, 5, 7>(mesh_file);
            break;

            case vtk_element_number::QUADRATIC_LAGRANGE:
               elements_2d[e] = read_element<0, 2, 4, 6, 1, 3, 5, 7, 8>(mesh_file);
            break;

            default:
                throw std::domain_error{"Unknown 2D element."};
        }
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
auto mesh_container_2d<T, I>::read_elements_1d(Stream& mesh_file) {
    std::string pass;
    size_t groups_count = 0;
    mesh_file >> pass >> groups_count;
    std::vector<std::string> groups_names_1d(groups_count);
    std::vector<std::vector<std::vector<I>>> elements_1d(groups_count);
    std::vector<std::vector<uint8_t>> elements_types_1d(groups_count);
    for(const size_t b : std::ranges::iota_view{0u, groups_count}) {
        size_t elements_count = 0;
        mesh_file >> pass >> groups_names_1d[b] >> pass >> elements_count;
        elements_1d[b].resize(elements_count);
        elements_types_1d[b].resize(elements_count);
        for(const size_t e : std::ranges::iota_view{0u, elements_count}) {
            size_t type = 0;
            mesh_file >> type;
            elements_types_1d[b][e] = get_elements_set().model_to_local_1d(type);
            switch(vtk_element_number(type)) {
                case vtk_element_number::LINEAR:
                    elements_1d[b][e] = read_element<0, 1>(mesh_file);
                break;

                case vtk_element_number::QUADRATIC:
                    elements_1d[b][e] = read_element<0, 2, 1>(mesh_file);
                break;

                default:
                    throw std::domain_error{"Unknown 1D element."};
            }
        }
    }
    return std::make_tuple(std::move(groups_names_1d), std::move(elements_1d), std::move(elements_types_1d));
}

template<class T, class I>
template<class Stream>
void mesh_container_2d<T, I>::read_su2(Stream& mesh_file) {
    auto [elements_2d, elements_types_2d] = read_elements_2d(mesh_file);
    auto nodes = read_nodes(mesh_file);
    auto [groups_names_1d, elements_1d, elements_types_1d] = read_elements_1d(mesh_file);

    _nodes = std::move(nodes);
    _groups_names_1d = std::move(groups_names_1d);
    _groups_names_2d = {"Default"};
    _elements_2d_count = elements_2d.size();
    _elements_groups["Default"] = {0u, _elements_2d_count};
    size_t elements = _elements_2d_count;
    for(const size_t b : std::ranges::iota_view{0u, _groups_names_1d.size()}) {
        _elements_groups[_groups_names_1d[b]] = std::ranges::iota_view{elements, elements + elements_1d[b].size()};
        elements += elements_1d[b].size();
    }
    
    _elements.resize(elements);
    _elements_types.resize(elements);
    for(const size_t e : std::ranges::iota_view{0u, _elements_2d_count}) {
        _elements[e] = std::move(elements_2d[e]);
        _elements_types[e] = elements_types_2d[e];
    }

    size_t curr_element = _elements_2d_count;
    for(const size_t b : std::ranges::iota_view{0u, _groups_names_1d.size()})
        for(const size_t e : std::ranges::iota_view{0u, elements_1d[b].size()}) {
            _elements[curr_element] = std::move(elements_1d[b][e]);
            _elements_types[curr_element] = elements_types_1d[b][e];
            ++curr_element;
        }
}

}

#endif