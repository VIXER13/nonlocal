#include "config_utils.hpp"

#include <logger/logger.hpp>

#include <iostream>
#include <fstream>
#include <exception>

namespace {

template<bool Is_Optional>
void check_fields(const nlohmann::json& value, const std::vector<std::string>& fields, const std::string& path) {
    std::string message;
	for(const std::string& field : fields)
		if (!value.contains(field))
			message += "Field \"" + path + field + "\" is missed.\n";
	if (!message.empty()) {
        if constexpr (Is_Optional)
            logger::debug() << "Some optional fields are missing:\n" + message << std::flush;
        else {
            message.pop_back(); // remove redudant '\n'
            throw std::domain_error{"Some required fields are missing:\n" + message};
        }
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
    static constexpr bool Is_Optional = false;
    check_fields<Is_Optional>(value, required, path);
}

void check_optional_fields(const nlohmann::json& value, const std::vector<std::string>& optional, const std::string& path) {
    static constexpr bool Is_Optional = true;
    check_fields<Is_Optional>(value, optional, path);
}

}