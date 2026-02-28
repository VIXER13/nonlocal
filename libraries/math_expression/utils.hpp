#pragma once

#include <metamath/types/traits.hpp>

#include <cstdint>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

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
const std::unordered_map<std::string, std::size_t>& get_operator_priority();
std::unordered_map<std::string, std::size_t> get_variables(const std::span<token_t>& tokens);
void check_variables_admissibility(const std::unordered_map<std::string, std::size_t>& variables);

template<metamath::types::arithmetic T>
T get_number(const std::string& number, std::size_t* idx = nullptr, int base = 10) {
    if constexpr (std::is_same_v<T, float>)       
        return std::stof(number, idx);
    if constexpr (std::is_same_v<T, double>)      
        return std::stod(number, idx);
    if constexpr (std::is_same_v<T, long double>) 
        return std::stold(number, idx);
    if constexpr (std::is_same_v<T, int>)         
        return std::stoi(number, idx, base);
    if constexpr (std::is_same_v<T, long>)         
        return std::stol(number, idx, base);
    if constexpr (std::is_same_v<T, long long>)         
        return std::stoll(number, idx, base);
    if constexpr (std::is_same_v<T, unsigned long>) 
        return std::stoul(number, idx, base);
    if constexpr (std::is_same_v<T, unsigned long long>) 
        return std::stoull(number, idx, base);
    return std::stoi(number, idx, base); // bool, char 
}

}
