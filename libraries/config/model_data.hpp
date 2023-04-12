#ifndef NONLOCAL_MODEL_DATA_HPP
#define NONLOCAL_MODEL_DATA_HPP

#include <json/value.h>

#include <type_traits>
#include <array>
#include <ranges>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
class model_data final {
    using radius_t = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;

    static radius_t read_radius(const Json::Value& arr, const std::string& field) {
        std::array<T, Dimension> result;
    
        if (arr.isDouble())
            result.fill(arr.template as<T>());
        else if (arr.isArray() && arr.size() == Dimension)
            for(const Json::ArrayIndex i : std::ranges::iota_view{0u, Dimension})
                result[i] = arr[i].template as<T>();
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
    explicit model_data(const Json::Value& model) {
        check_required_fields(model, { "local_weight", "nonlocal_radius" });
        local_weight = model["local_weight"].template as<T>();
        nonlocal_radius = read_radius(model["nonlocal_radius"], "nonlocal_radius");
        search_radius = !model.isMember("search_radius") ? nonlocal_radius :
                        read_radius(model["search_radius"], "search_radius");
    }

    Json::Value to_json() const {
        Json::Value result;
        result["local_weight"] = local_weight;
        if constexpr (std::is_same_v<radius_t, T>) {
            result["nonlocal_radius"] = nonlocal_radius;
            result["search_radius"] = search_radius;
        } else {
            Json::Value& nonloc = result["nonlocal_radius"] = Json::arrayValue;
            Json::Value& search = result["search_radius"] = Json::arrayValue;
            for(const size_t i : std::ranges::iota_view{0u, Dimension}) {
                nonloc.append(nonlocal_radius[i]);
                search.append(search_radius[i]);
            }
        }
        return result;
    }
};

}

#endif