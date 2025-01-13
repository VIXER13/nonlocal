#pragma once

#include <string_view>

namespace metamath::utils {

template<class T>
constexpr std::string_view type_name() {
#if defined(__clang__)
    constexpr std::string_view prefix = "[T = ";
    constexpr std::string_view suffix = "]";
    constexpr std::string_view function = __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
    constexpr std::string_view prefix = "with T = ";
    constexpr std::string_view suffix = "; ";
    constexpr std::string_view function = __PRETTY_FUNCTION__;
#elif defined(__MSC_VER)
    constexpr std::string_view prefix = "GetTypeName<";
    constexpr std::string_view suffix = ">(void)";
    constexpr std::string_view function = __FUNCSIG__;
#else
#error Unsupported compiler
#endif
    const auto start = function.find(prefix) + prefix.size();
    const auto end = function.find(suffix);
    const auto size = end - start;
    return function.substr(start, size);
}

template<class T>
constexpr std::string_view type_name(const T&) {
    return type_name<T>();
}

}