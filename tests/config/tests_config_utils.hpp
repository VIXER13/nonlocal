#pragma once

#include <nlohmann/json.hpp>

namespace unit_tests {

template<class T, size_t N>
class mock_data final {
    nlohmann::json _value;

public:
    static constexpr std::string_view Prefix = "";

    explicit constexpr mock_data() noexcept = default;
    explicit constexpr mock_data(const nlohmann::json& value, const std::string&) 
        : _value(value) {}

    constexpr operator nlohmann::json() const {
        return _value;
    }
};

}