#pragma once

#include "mesh_2d.hpp"

namespace nonlocal::mesh {

enum class diam_adding : uint8_t { NO, MAX, MIN, MEAN };

template<class T>
struct rectangle final {
    T left  = std::numeric_limits<T>::max();
    T down  = std::numeric_limits<T>::max();
    T right = -std::numeric_limits<T>::max();
    T up    = -std::numeric_limits<T>::max();

    T length() const noexcept {
        return right - left;
    }

    T width() const noexcept {
        return up - down;
    }

    void expand_corners(const std::array<T, 2>& node) noexcept {
        left  = std::min(left,  node[0]);
        down  = std::min(down,  node[1]);
        right = std::max(right, node[0]);
        up    = std::max(up,    node[1]);
    }
};

template<class T, class I>
rectangle<T> mesh_corners(const mesh_container_2d<T, I>& mesh) {
    rectangle<T> corners;
    for(const size_t node : mesh.nodes())
        corners.expand_corners(mesh.node_coord(node));
    return corners;
}

template<class T, class I>
rectangle<T> group_corners(const mesh_container_2d<T, I>& mesh, const std::string& group) {
    rectangle<T> corners;
    for(const size_t e : mesh.elements(group))
        for(const size_t node : mesh.nodes(e))
            corners.expand_corners(mesh.node_coord(node));
    return corners;
}

template<class T, class I>
T diam_between_elements(const mesh_2d<T, I>& mesh, 
                        const std::string& group, 
                        const std::vector<std::array<T, 2>>& centers,
                        const diam_adding diam) {
    if (diam == diam_adding::NO)
        return T{0};
    std::unordered_map<diam_adding, T> diams = {
        {diam_adding::MAX, T{0}},
        {diam_adding::MIN, std::numeric_limits<T>::max()},
        {diam_adding::MEAN, T{0}}
    };
    std::set<std::pair<I, I>> pairs;
    for (const size_t element : mesh.container().elements(group))
        for (const size_t node : mesh.container().nodes(element))
            for (const size_t adjacent_element : mesh.elements(node)) {
                const T distance = metamath::functions::powered_distance<2>(centers[element], centers[adjacent_element]);
                diams[diam_adding::MAX] = std::max(diams[diam_adding::MAX], distance);
                diams[diam_adding::MIN] = std::min(diams[diam_adding::MIN], distance);
                if (diam == diam_adding::MEAN && !pairs.contains({element, adjacent_element}) && 
                                                 !pairs.contains({adjacent_element, element})) {
                    pairs.insert({element, adjacent_element});
                    diams[diam_adding::MEAN] += std::sqrt(distance);
                }
            }
    const T result = diams.at(diam);
    return diam == diam_adding::MEAN ? result / pairs.size() : std::sqrt(result);
}

template<class T, class I>
std::unordered_map<std::string, T> search_radii(const mesh_2d<T, I>& mesh, 
                                                const std::unordered_map<std::string, T>& nonlocal_radii,
                                                const std::vector<std::array<T, 2>>& centers,
                                                const diam_adding diam) {
    std::unordered_map<std::string, T> radii;
    for(const auto& [group, radius] : nonlocal_radii)
        radii[group] = radius + diam_between_elements(mesh, group, centers, diam);
    return radii;
}

template<class T>
std::pair<size_t, size_t> subarea_id(const rectangle<T>& corners, const std::array<T, 2>& center, const std::array<T, 2>& radius) {
    return {
        (center[0] - corners.left) / radius[0],
        (center[1] - corners.down) / radius[1]
    };
}

template<class I, class T>
metamath::types::matrix<std::vector<I>> split_elements_by_subareas(const std::vector<std::array<T, 2>>& centers,
                                                                   const std::ranges::iota_view<size_t, size_t>& elements, 
                                                                   const rectangle<T>& corners,
                                                                   const T radius) {
    metamath::types::matrix<std::vector<I>> subareas(size_t(corners.length() / radius) + 1u, size_t(corners.width() / radius) + 1u);
    for(const size_t e : elements) {
        const auto [row, col] = subarea_id(corners, centers[e], {radius, radius});
        subareas(row, col).push_back(e);
    }
    return subareas;
}

template<class I>
std::vector<std::pair<size_t, size_t>> subareas_ids(const metamath::types::matrix<std::vector<I>>& subareas, const size_t row, const size_t col) {
    std::vector<std::pair<size_t, size_t>> result;
    result.reserve(9);
    for(const size_t i : std::ranges::iota_view{row ? row - 1 : row, row != subareas.rows() - 1 ? row + 2 : row + 1})
        for(const size_t j : std::ranges::iota_view{col ? col - 1 : col, col != subareas.cols() - 1 ? col + 2 : col + 1})
            result.push_back({i, j});
    return result;
}

template<class T, class I>
neighbours_t<T, I> find_neighbours(const mesh_2d<T, I>& mesh, const std::unordered_map<std::string, T>& nonlocal_radii, const diam_adding add_diam = diam_adding::MAX) {
    std::vector<std::vector<I>> neighbours(mesh.container().elements_2d_count());
    if (nonlocal_radii.empty())
        return {nonlocal_radii, neighbours};
    const std::unordered_set<I> process_elements = mesh.process_elements();
    const std::vector<std::array<T, 2>> centers = utils::approx_centers_of_elements(mesh.container());
    const std::unordered_map<std::string, T> radii = search_radii(mesh, nonlocal_radii, centers, add_diam);
    for(const auto& [group, radius] : radii) {
        if (radius <= T{0}) {
            logger::get().log(logger::log_level::WARNING)
                << "The search radius for the \"" << group << "\" group turned out to be less than 0" << std::endl;
            continue;
        }
        const auto elements_range = mesh.container().elements(group);
        const auto corners = group_corners(mesh.container(), group);
        const auto subareas = split_elements_by_subareas<I>(centers, elements_range, corners, radius);
// Parallelism is disabled because it results in a high degree of data fragmentation in memory. Enable at your own risk
//#pragma omp parallel for
        for(size_t k = 0; k < elements_range.size(); ++k) {
            const size_t eL = elements_range[k];
            if (!process_elements.contains(eL))
                continue;

            const auto [row, col] = subarea_id(corners, centers[eL], {radius, radius});
            const auto ids = subareas_ids(subareas, row, col);
            const auto reserver = [&subareas](const size_t sum, const std::pair<size_t, size_t>& id) {
                return sum + subareas(id.first, id.second).size();
            };
            neighbours[eL].reserve(std::accumulate(ids.begin(), ids.end(), size_t{0}, reserver));

            for(const auto [i, j] : ids)
                for(const size_t eNL : subareas(i, j))
                    if (metamath::functions::distance(centers[eL], centers[eNL]) <= radius)
                        neighbours[eL].push_back(eNL);
            neighbours[eL].shrink_to_fit();
        }
    }
    return {std::move(radii), std::move(neighbours)};
}

}