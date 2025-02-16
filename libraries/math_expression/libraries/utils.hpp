#pragma once

#include <string>
#include <concepts>

namespace formula::utils {

template <typename T> 
concept arithmetic = std::integral<T> || std::floating_point<T>;

std::string trim(std::string str);
std::string& delete_all(std::string& str, char symb);
std::size_t count_all(const std::string& s, char symb);

bool is_latin_str(const std::string& s);
bool is_number(const std::string& s);

template<arithmetic T>
T get_number(const std::string& number, std::size_t* idx = 0, int base = 10) {
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
