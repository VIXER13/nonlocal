#pragma once

#include "mesh_container_2d.hpp"

#include <metamath/types/traits.hpp>
#include <constants/nonlocal_constants.hpp>

namespace nonlocal::mesh::utils {

template<std::floating_point T>
using data_pair_t = std::variant<
    std::pair<std::string, std::reference_wrapper<const std::vector<T>>>,
    std::pair<std::array<std::string, 2>, std::reference_wrapper<const std::vector<std::array<T, 2>>>>,
    std::pair<std::array<std::string, 3>, std::reference_wrapper<const std::vector<std::array<T, 3>>>>
>;

template<class T>
std::string to_string(const T& names) {
    using U = std::remove_cvref_t<decltype(names)>;
    if constexpr (std::is_same_v<U, std::string>)
        return names;
    else {
        std::string result;
        static constexpr size_t N = std::tuple_size_v<T>;
        for(const size_t i : std::ranges::iota_view{0zu, N})
            result += names[i] + (i == N - 1 ? "" : ",");
        return result;
    }
}

template<std::floating_point T, std::integral I>
void validate(const mesh_container_2d<T, I>& mesh, const std::vector<data_pair_t<T>>& data) {
    for(const auto& pair : data)
        std::visit([&mesh](const auto& data) {
            if (const auto& [name, vec] = data; vec.get().size() != mesh.nodes_count())
                throw std::logic_error{"The result cannot be saved to csv because the mesh nodes number "
                                       "and elements in the vector \"" + to_string(name) + "\" do not match."};
        }, pair);
}

template<std::floating_point T, std::integral I>
void save_as_csv(const std::filesystem::path& path, const mesh_container_2d<T, I>& mesh,
                 const std::vector<data_pair_t<T>>& data, const std::optional<std::streamsize> precision = std::nullopt) {
    validate(mesh, data);
    std::ofstream csv{path};
    csv.precision(precision ? *precision : std::numeric_limits<T>::max_digits10);

    // Header
    csv << "x,y" << (data.empty() ? '\n' : ',');
    for(const size_t j : std::ranges::iota_view{0u, data.size()}) {
        std::visit([&csv](const auto& data) { csv << to_string(data.first); }, data[j]);
        csv << (j == data.size() - 1 ? '\n' : ',');
    }

    // Data
    for(const size_t node : std::ranges::iota_view{0u, mesh.nodes_count()}) {
        const auto& [x, y] = mesh.node_coord(node);
        csv << x << ',' << y << (data.empty() ? '\n' : ',');
        for(const size_t j : std::ranges::iota_view{0u, data.size()}) {
            std::visit([&csv, node](const auto& data) {
                using U = std::remove_cvref_t<decltype(data.second.get())>;
                const auto& vec = data.second.get();
                if constexpr (std::is_same_v<U, std::vector<T>>)
                    csv << vec[node];
                else if constexpr (std::is_same_v<U, std::vector<std::array<T, 2>>>)
                    csv << vec[node][X] << ',' << vec[node][Y];
                else if constexpr (std::is_same_v<U, std::vector<std::array<T, 3>>>)
                    csv << vec[node][X] << ',' << vec[node][Y] << ',' << vec[node][Z];
                else
                    static_assert(false, "Unsupported type for saving csv.");
            }, data[j]);
            csv << (j == data.size() - 1 ? '\n' : ',');
        }
    }
}

}