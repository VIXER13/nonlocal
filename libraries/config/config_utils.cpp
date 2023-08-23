#include "config_utils.hpp"

#include <iostream>
#include <fstream>
#include <exception>

namespace {

void check_fields(const nlohmann::json& value, const std::vector<std::string>& fields, const std::string& path, const bool is_required) {
    std::string message;
	for(const std::string& field : fields)
		if (!value.contains(field))
			message += "Field \"" + path + field + "\" missed.\n";
	if (!message.empty()) {
        message.pop_back(); // remove redudant '\n'
        if (is_required) throw std::domain_error{"Some required fields are missing:\n" + message};
        else std::cout << "INFO: Some optional fields are missing:\n" + message << std::endl;
    }
}

}

namespace nonlocal::config {

nlohmann::json parse_json(const std::filesystem::path& path) {
    std::ifstream file{path};
    return nlohmann::json::parse(file);
}

void dump_json(const nlohmann::json& value, const std::filesystem::path& path, const int indent, const char indent_char) {
    std::ofstream file{path};
    file << value.dump(indent, indent_char);
}

std::string append_access_sign(std::string path, const std::optional<size_t> index) {
    using namespace std::literals;
    return path += path.empty() ? ""s :
                   index        ? '[' + std::to_string(*index) + ']' : "."s;
}

void check_required_fields(const nlohmann::json& value, const std::vector<std::string>& required, const std::string& path) {
    check_fields(value, required, path, true);
}

void check_optional_fields(const nlohmann::json& value, const std::vector<std::string>& optional, const std::string& path) {
    check_fields(value, optional, path, false);
}

}