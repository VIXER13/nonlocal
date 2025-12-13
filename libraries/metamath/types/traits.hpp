#pragma once

#include <concepts>

namespace metamath::types {

template<class T> 
concept arithmetic = std::integral<T> || std::floating_point<T>;

}