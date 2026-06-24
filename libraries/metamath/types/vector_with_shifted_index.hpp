#pragma once

#include <vector>

namespace metamath::types {

template<class T>
struct vector_with_shifted_index final {
    std::vector<T> container;
    size_t shift = 0u;

    T& operator[](const size_t index) {
        return container[index - shift];
    }

    const T& operator[](const size_t index) const {
        return container[index - shift];
    }
};

}