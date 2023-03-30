#include "general_config_data.hpp"

#include <exception>

namespace {

template<class T>
T get(const Json::Value& type, const std::unordered_map<std::string, T>& types,
      const std::string_view error_not_string, const std::string_view error_invalid_type) {
    if (!type.isString())
        throw std::domain_error{error_not_string.data()};
    const std::string str = type.asString();
    if (const auto it = types.find(str); it != types.cend())
        return it->second;
    throw std::domain_error{error_invalid_type.data() + str};
}

}

namespace nonlocal::config {

save_data::save_data(const Json::Value& save) 
	: _folder{save.get("folder", ".").asString()} {
	if (const Json::Value& precision = save.get("precision", {}); precision.isIntegral())
		_precision = precision.asInt64();
	for (const std::string& field : save.getMemberNames())
		if (field != "folder" && field != "precision")
			_names[field] = save[field].asString();
}

std::optional<std::streamsize> save_data::precision() const noexcept {
	return _precision;
}

const std::filesystem::path& save_data::folder() const noexcept {
	return _folder;
}

bool save_data::contains(const std::string& key) const {
	return _names.contains(key);
}

std::string save_data::get_name(const std::string& key, const std::optional<std::string>& default_name) const {
	if (const auto it = _names.find(key); it != _names.cend())
		return it->second;
	else if (default_name)
		return *default_name;
	throw std::runtime_error{"Unknown key: " + key};
}

std::filesystem::path save_data::make_path(const std::string& name, const std::string& extension) const {
	return folder() / (name + extension);
}

std::filesystem::path save_data::path(const std::string& key, const std::string& extension,
									  const std::optional<std::string>& default_name) const {
	return make_path(get_name(key, default_name), extension);
}

Json::Value save_data::to_json() const {
	Json::Value result;
	result["folder"] = _folder.string();
	for(const auto& [key, name] : _names)
		result[key] = name;
	if (_precision)
		result["precision"] = *_precision;
	return result;
}

thermal::boundary_condition_t get_thermal_condition(const Json::Value& kind) {
    static const std::unordered_map<std::string, thermal::boundary_condition_t> kinds = {
        {"temperature", thermal::boundary_condition_t::TEMPERATURE},
        {"flux",        thermal::boundary_condition_t::FLUX},
        {"convection",  thermal::boundary_condition_t::CONVECTION},
        {"radiation",   thermal::boundary_condition_t::RADIATION},
        {"combined",    thermal::boundary_condition_t::COMBINED}
    };
    return kind.isIntegral() ?
           thermal::boundary_condition_t(kind.asUInt64()) :                     
           get(kind, kinds, "Boundary condition kind must be an integer or string.", "Unknown boundary condition type: ");
}

const std::string& get_thermal_condition(const thermal::boundary_condition_t kind) {
    static const std::unordered_map<thermal::boundary_condition_t, std::string> kinds {
        {thermal::boundary_condition_t::TEMPERATURE, "temperature"},
        {thermal::boundary_condition_t::FLUX,        "flux"},
        {thermal::boundary_condition_t::CONVECTION,  "convection"},
        {thermal::boundary_condition_t::RADIATION,   "radiation"},
        {thermal::boundary_condition_t::COMBINED,    "combined"}
    };
    if (const auto it = kinds.find(kind); it != kinds.cend())
        return it->second;
    throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(kind))};
}

size_t get_order(const Json::Value& order) {
    static const std::unordered_map<std::string, size_t> orders = {
        {"linear",    1},
        {"quadratic", 2},
        {"qubic",     3},
        {"quartic",   4},
        {"quintic",   5}
    };
    const size_t result = order.isIntegral() ?
                          order.asUInt64() :
                          get(order, orders, "Element and quadrature orders must be an integer or string.", "Invalid element or quadrature order: ");
    if (!result || result > 5)
        throw std::domain_error{"Invalid element or quadrature order: " + std::to_string(result)};
    return result;
}

const std::string& get_order(const size_t order) {
    static const std::array<std::string, 5> orders = {
        "linear", "quadratic", "qubic", "quartic", "quintic"
    };
    if (const size_t ind = order - 1; ind < orders.size())
        return orders[ind];
    throw std::domain_error{"Unknown order: " + std::to_string(order)};
}

material_t get_material(const Json::Value& material) {
    static const std::unordered_map<std::string, material_t> materials = {
        {"isotropic",   material_t::ISOTROPIC},
        {"orthotropic", material_t::ORTHOTROPIC},
        {"anisotropic", material_t::ANISOTROPIC}
    };
    return get(material, materials, "Material type must be a string.", "Invalid material type: ");
}

const std::string& get_material(const material_t material) {
    static const std::array<std::string, 3> materials = {
        "isotropic", "orthotropic", "anisotropic"
    };
    return materials[size_t(material)];
}

}