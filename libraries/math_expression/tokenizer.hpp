#pragma once

#include <ostream>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <span>

namespace formula::utils {
struct token_t {
    enum class type_t : uint8_t {
        Unknown,
        Number,
        ParenthesisLeft,
        ParenthesisRight,
        Separator,
        Operator,
        Symbol,
    };

    type_t type = type_t::Unknown;
    std::string str;

    friend std::ostream& operator<<(std::ostream& os, token_t token);
};

std::vector<token_t> tokenize(std::string_view input);
std::unordered_map<std::string, std::size_t> get_variables(const std::span<token_t>& tokens);
}  // namespace formula::utils