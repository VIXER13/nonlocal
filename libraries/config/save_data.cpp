#include "save_data.hpp"

namespace nonlocal::config {

save_data::save_data(const nlohmann::json& save) 
	: _folder{save.value("folder", std::filesystem::current_path().string())} {
	if (const nlohmann::json precision = save.value("precision", nlohmann::json{}); precision.is_number_integer())
		_precision = precision.get<std::streamsize>();
	for(const auto& [field, value] : save.items())
		if (field != "folder" && field != "precision")
			_names[field] = value.get<std::string>();
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

save_data::operator nlohmann::json() const {
	nlohmann::json result = _names;
	result["folder"] = _folder.string();
	if (_precision)
		result["precision"] = std::make_unsigned_t<std::streamsize>(*_precision);
	return result;
}

}