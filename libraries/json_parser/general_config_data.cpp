#include "general_config_data.hpp"

#include <exception>

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

}