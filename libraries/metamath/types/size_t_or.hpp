#pragma once

#include <variant>

namespace metamath::types {

template<class T>
using size_t_or = std::variant<size_t, T>;

template<class T>
T get_value(const size_t_or<T>& value) {
    if (std::holds_alternative<T>(value))
        return std::get<T>(value);
    return T(std::get<size_t>(value));
}

}
