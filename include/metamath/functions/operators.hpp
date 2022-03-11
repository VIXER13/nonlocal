#ifndef METAMATH_OPERATORS_HPP
#define METAMATH_OPERATORS_HPP

#include <type_traits>
#include <stddef.h>

namespace metamath::function {

template<class Container>
Container& operator+=(Container& lhs, const Container& rhs) {
    for(size_t i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

template<class Container, class T>
std::enable_if_t<std::is_arithmetic_v<T>, Container&> operator*=(Container& container, const T& val) noexcept {
    for(T& it : container)
        it *= val;
    return container;
}

template<class T, class Container>
std::enable_if_t<std::is_arithmetic_v<T>, Container> operator*(const T& val, const Container& container) {
    Container new_container = container;
    new_container *= val;
    return new_container;
}

template<class Container, class T>
std::enable_if_t<std::is_arithmetic_v<T>, Container> operator*(const Container& container, const T& val) {
    return val * container;
}

template<class Container, class T>
std::enable_if_t<std::is_arithmetic_v<T>, Container&> operator/=(Container& container, const T& val) noexcept {
    for(T& it : container)
        it /= val;
    return container;
}

template<class Container, class T>
std::enable_if_t<std::is_arithmetic_v<T>, Container> operator/(const Container& container, const T& val) {
    Container new_container = container;
    new_container /= val;
    return container;
}

}

#endif