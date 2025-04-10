#include "utils.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace formula::utils {

std::string trim(std::string str)
{
    str.erase(str.find_last_not_of(' ')+1);    //suffixing spaces
    str.erase(0, str.find_first_not_of(' '));  //prefixing spaces
    return str;
}

std::string& delete_all(std::string& str, char symb) {
    str.erase(std::remove(str.begin(), str.end(), symb), str.end());
    return str;
}

std::size_t count_all(const std::string& s, char symb) {
    return std::accumulate(s.begin(), s.end(), std::size_t(0), [&](std::size_t res, char a){ return res + std::size_t(a == symb); });
}

bool is_latin_str(const std::string& s) {
    if (s.size() == 0 || s.size() > 1)
        throw std::domain_error{"Wrong format, latyn symbol contains of one symbol."};
    for (const char d : s)
        if (!std::isalpha(d)) 
            return false;
    return true;
};

bool is_number(const std::string& s) {
    if (const std::size_t dots = count_all(s, '.'); dots > 1)
        return false;
    for (const char d : s)
        if (d != '.' && !std::isdigit(d)) 
            return false;
    return true;
};

}