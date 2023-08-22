#include "config_utils.hpp"

#include <iostream>
#include <fstream>
#include <exception>

namespace {

bool has_value(const std::variant<std::string, size_t>& index) {
    if (std::holds_alternative<size_t>(index))
        return true;
    return !std::get<std::string>(index).empty();
}

std::string index_value(const std::variant<std::string, size_t>& index) {
    std::string result = "[";
    if (std::holds_alternative<size_t>(index))
        result += std::to_string(std::get<size_t>(index));
    else
        result += std::get<std::string>(index);
    return result + "]";
}

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

std::string append_access_sign(std::string path, const std::variant<std::string, size_t>& index) {
    return path += path.empty()     ? std::string{}      :
                   has_value(index) ? index_value(index) : std::string{'.'};
}

void check_required_fields(const nlohmann::json& value, const std::vector<std::string>& required, const std::string& path) {
    check_fields(value, required, path, true);
}

void check_optional_fields(const nlohmann::json& value, const std::vector<std::string>& optional, const std::string& path) {
    check_fields(value, optional, path, false);
}

}