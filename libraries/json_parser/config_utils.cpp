#include "config_utils.hpp"

#include <fstream>
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

Json::Value read_json_from_file(const std::filesystem::path& config) {
	std::ifstream file{config};
    Json::Reader reader;
	Json::Value value;
	if (reader.parse(file, value))
        return value;
	throw std::runtime_error{"Invalid json: " + config.string() + ".\nError message: " + reader.getFormattedErrorMessages()};
}

void check_required_fields(const Json::Value& value, const std::vector<std::string>& required) {
	std::string error_message;
	for(const std::string& field : required)
		if (!value.isMember(field))
			error_message += "Field \"" + field + "\" missed.\n";
	if (!error_message.empty())
		throw std::domain_error{"Some required fields are missing:\n" + error_message};
}
    
}