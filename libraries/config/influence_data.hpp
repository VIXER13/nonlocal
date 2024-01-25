#ifndef NONLOCAL_CONFIG_INFLUENCE_DATA_HPP
#define NONLOCAL_CONFIG_INFLUENCE_DATA_HPP

#include "config_utils.hpp"

namespace nonlocal::config {

enum class influence_family_t : uint8_t {
    CUSTOM,
    CONSTANT,
    POLYNOMIAL,
    EXPONENTIAL
};

NLOHMANN_JSON_SERIALIZE_ENUM(influence_family_t, {
    {influence_family_t::CUSTOM, nullptr},
    {influence_family_t::CONSTANT, "constant"},
    {influence_family_t::POLYNOMIAL, "polynomial"},
    {influence_family_t::EXPONENTIAL, "exponential"}
})

template<std::floating_point T, size_t Dimension>
struct influence_data final {
    influence_family_t family = influence_family_t::POLYNOMIAL;
    size_t n = 2;
    size_t p = 2;
    size_t q = 1;
    std::string custom;

    explicit constexpr influence_data() noexcept = default;
    explicit influence_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_optional_fields(config, { "family" }, path_with_access);
        family = config["family"].get<influence_family_t>();
        if (family == influence_family_t::CUSTOM) {
            custom = config["custom"].get<std::string>();
        } else {
            if constexpr (Dimension != 1) {
                check_required_fields(config, { "n" }, path_with_access);
                n = config["n"].is_number_unsigned() ?
                    config["n"].get<size_t>() :
                    config["n"].is_string() && config["n"].get<std::string>() == "infinity" ?
                    std::numeric_limits<size_t>::max() : 
                    throw std::domain_error{"Unexpected n value."};
            }
            if (family != influence_family_t::CONSTANT) {
                check_required_fields(config, { "p", "q" }, path_with_access);
                p = config["p"].get<size_t>();
                q = config["q"].get<size_t>();
            }
        }
    }

    operator nlohmann::json() const {
        if (family == influence_family_t::CUSTOM)
            return {{"custom", custom}};
        nlohmann::json result = {{"family", family}};
        if (n == std::numeric_limits<size_t>::max())
            result["n"] = "infinity";
        else
            result["n"] = n;
        if (family != influence_family_t::CONSTANT) {
            result["p"] = p;
            result["q"] = q;
        }
        return result;
    }
};

}

#endif