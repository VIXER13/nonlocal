#include "general_config_data.hpp"

#include <exception>

namespace nonlocal::config {

save_data::save_data(const Json::Value& save) 
	: _folder{save.get("folder", ".").asString()} {
	for (const std::string& field : save.getMemberNames())
		if (field != "folder")
			_names[field] = save[field].asString();
}

bool save_data::contains(const std::string& key) const {
	return _names.contains(key);
}

const std::filesystem::path& save_data::folder() const noexcept {
	return _folder;
}

std::filesystem::path save_data::path(const std::string& key, const std::string& extension,
									  const std::optional<std::string>& default_name) const {
	if (const auto it = _names.find(key); it != _names.cend())
		return folder() / (it->second + extension);
	else if (default_name)
		return folder() / (*default_name + extension);
	throw std::domain_error{"Unknown key: " + key};
}

}