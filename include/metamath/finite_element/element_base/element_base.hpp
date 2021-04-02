#ifndef FINITE_ELEMENT_BASE_HPP
#define FINITE_ELEMENT_BASE_HPP

#include <cstddef>

namespace metamath::finite_element {

class element_base {
public:
    virtual size_t nodes_count() const = 0; // Любой конечный элемент, вне зависимости от его размерности, имеет некоторое количество узлов.
    virtual ~element_base() noexcept = default;
};

}

#endif