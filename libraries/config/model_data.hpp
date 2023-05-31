#ifndef NONLOCAL_MODEL_DATA_HPP
#define NONLOCAL_MODEL_DATA_HPP

#include "config_utils.hpp"

#include <type_traits>
#include <array>
#include <ranges>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
class model_data final {
    using radius_t = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;

    static radius_t read_radius(const nlohmann::json& data, const std::string& field) {
        std::array<T, Dimension> result;
        if (data.is_number())
            result.fill(data.get<T>());
        else if (data.is_array() && data.size() == Dimension)
            for(const size_t i : std::ranges::iota_view{0u, Dimension})
                result[i] = data[i].get<T>();
        else
            throw std::domain_error{"Field \"" + field + "\" must be an array with length " + std::to_string(Dimension)};

        if constexpr (std::is_same_v<radius_t, T>)
            return result.front();
        else
            return result;
    }

public:
    T local_weight = T{1};             // required
    radius_t nonlocal_radius = {T{0}}; // required
    radius_t search_radius = {T{0}};   // if skipped sets equal nonlocal_radius

    explicit constexpr model_data() noexcept = default;
    explicit model_data(const nlohmann::json& model) {
        check_required_fields(model, { "local_weight", "nonlocal_radius" });
        local_weight = model["local_weight"].get<T>();
        nonlocal_radius = read_radius(model["nonlocal_radius"], "nonlocal_radius");
        search_radius = !model.contains("search_radius") ? nonlocal_radius :
                        read_radius(model["search_radius"], "search_radius");
    }

    operator nlohmann::json() const {
        return {
            {"local_weight", local_weight},
            {"nonlocal_radius", nonlocal_radius},
            {"search_radius", search_radius}
        };
    }
};

}

#endif