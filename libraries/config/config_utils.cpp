#include "config_utils.hpp"

#include <fstream>
#include <exception>
#include <unordered_map>

namespace {

template<class T>
T get(const nlohmann::json& type, const std::unordered_map<std::string, T>& types,
      const std::string_view error_not_string, const std::string_view error_invalid_type) {
    if (!type.is_string())
        throw std::domain_error{error_not_string.data()};
    const std::string str = type.get<std::string>();;
    if (const auto it = types.find(str); it != types.cend())
        return it->second;
    throw std::domain_error{error_invalid_type.data() + str};
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

void check_required_fields(const nlohmann::json& value, const std::vector<std::string>& required) {
    std::string error_message;
	for(const std::string& field : required)
		if (!value.contains(field))
			error_message += "Field \"" + field + "\" missed.\n";
	if (!error_message.empty())
		throw std::domain_error{"Some required fields are missing:\n" + error_message};
}

size_t get_order(const nlohmann::json& order) {
    static const std::unordered_map<std::string, size_t> orders = {
        {"linear",    1},
        {"quadratic", 2},
        {"qubic",     3},
        {"quartic",   4},
        {"quintic",   5}
    };
    const size_t result = order.is_number_integer() ?
                          order.get<size_t>() :
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

}