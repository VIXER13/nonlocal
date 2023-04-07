#include "config_utils.hpp"

#include <fstream>
#include <exception>

namespace nonlocal::config {

Json::Value read_json(const std::filesystem::path& path) {
	std::ifstream file{path};
    Json::Reader reader;
	Json::Value value;
	if (reader.parse(file, value))
        return value;
	throw std::runtime_error{"Invalid json file: " + path.string() + ".\nError message: " + reader.getFormattedErrorMessages()};
}

Json::Value read_json(const char *const begin, const char *const end) {
	Json::Reader reader;
	Json::Value value;
	if (reader.parse(begin, end, value))
        return value;
	throw std::runtime_error{"Invalid json: " + reader.getFormattedErrorMessages()};
}

Json::Value read_json(const std::string& str) {
	return read_json(str.data(), std::next(str.data(), str.size()));
}

void save_json(const std::filesystem::path& path, const Json::Value& value) {
	Json::StyledStreamWriter writer;
	std::ofstream file{path};
	writer.write(file, value);
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