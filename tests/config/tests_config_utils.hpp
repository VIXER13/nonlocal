#ifndef UNIT_TESTS_CONFIG_UTILS_HPP
#define UNIT_TESTS_CONFIG_UTILS_HPP

#include <nlohmann/json.hpp>

namespace unit_tests {

template<class T, size_t N>
class mock_data final {
    nlohmann::json _value;

public:
    explicit constexpr mock_data() noexcept = default;
    explicit constexpr mock_data(const nlohmann::json& value, const std::string&) 
        : _value(value) {}

    constexpr operator nlohmann::json() const {
        return _value;
    }
};

}

#endif