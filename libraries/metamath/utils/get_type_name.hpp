#ifndef METAMATH_GET_TYPE_NAME_HPP
#define METAMATH_GET_TYPE_NAME_HPP

#include <string_view>

namespace metamath::utils {

template<class T>
constexpr std::string_view get_type_name() {
#if defined(__clang__)
    constexpr auto prefix = std::string_view{"[T = "};
    constexpr auto suffix = "]";
    constexpr auto function = std::string_view{__PRETTY_FUNCTION__};
#elif defined(__GNUC__)
    constexpr auto prefix = std::string_view{"with T = "};
    constexpr auto suffix = "; ";
    constexpr auto function = std::string_view{__PRETTY_FUNCTION__};
#elif defined(__MSC_VER)
    constexpr auto prefix = std::string_view{"GetTypeName<"};
    constexpr auto suffix = ">(void)";
    constexpr auto function = std::string_view{__FUNCSIG__};
#else
#   error Unsupported compiler
#endif
    const auto start = function.find(prefix) + prefix.size();
    const auto end = function.find(suffix);
    const auto size = end - start;
    return function.substr(start, size);
}

template<class T>
constexpr std::string_view get_type_name(const T&) {
    return get_type_name<T>();
}

}

#endif