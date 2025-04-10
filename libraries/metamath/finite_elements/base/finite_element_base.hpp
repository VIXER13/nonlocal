#pragma once

#include <cstddef>
#include <ranges>

namespace metamath::finite_element {

class element_base {
public:
    virtual ~element_base() noexcept = default;
    virtual size_t nodes_count() const = 0;

    std::ranges::iota_view<size_t, size_t> nodes() const { return {0u, nodes_count()}; }
};

}