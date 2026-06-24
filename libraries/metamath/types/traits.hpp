#pragma once

#include <concepts>
#include <memory>

namespace metamath::types {

template<class T> 
concept arithmetic = std::integral<T> || std::floating_point<T>;

template<class T>
concept copyable = requires(const T& v) { { v.copy() } -> std::convertible_to<std::unique_ptr<T>>; };

}