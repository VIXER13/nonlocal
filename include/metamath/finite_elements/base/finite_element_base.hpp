#ifndef FINITE_ELEMENT_BASE_HPP
#define FINITE_ELEMENT_BASE_HPP

#include <cstddef>

namespace metamath::finite_element {

class element_base {
public:
    virtual size_t nodes_count() const = 0;
    virtual ~element_base() noexcept = default;
};

}

#endif