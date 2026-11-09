#pragma once

#include "read_model.hpp"

namespace nonlocal::config {

template<std::floating_point T, template<class, auto...> class Physical>
std::vector<equation_parameters<1u, T, Physical>> read_materials_1d(
    const nlohmann::json& config, const std::string& path, const std::string& prefix,
    Physical<T> (*read_physical)(const nlohmann::json&, const std::string&)) {
    if (!config.is_array())
        throw std::domain_error{ "Field \"" + path + "\" must be an array of materials." };

    std::vector<equation_parameters<1u, T, Physical>> parameters;
    parameters.reserve(config.size());
    for (const size_t i : std::ranges::iota_view{ 0u, config.size() }) {
        const auto& material = config[i];
        const std::string material_path = append_access_sign(append_access_sign(path, i));
        check_required_fields(material, { "physical" }, material_path);
        const std::string model_field = get_model_field(material, material_path, prefix);
        parameters.push_back({ .model = model_field.empty() ? model_parameters<1u, T>{}
                                                            : read_model_1d<T>(material[model_field], material_path + model_field),
                               .physical = read_physical(material["physical"], material_path + "physical") });
    }
    return parameters;
}

} // namespace nonlocal::config
