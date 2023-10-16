#ifndef FINITE_ELEMENT_QUADRATURE_BASE_HPP
#define FINITE_ELEMENT_QUADRATURE_BASE_HPP

#include <cstddef>
#include <type_traits>
#include <ranges>

namespace metamath::finite_element {

template<class T>
class quadrature_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

public:
    virtual ~quadrature_base() noexcept = default;

    virtual size_t nodes_count() const = 0;
    virtual T weight(const size_t i) const = 0;

    std::ranges::iota_view<size_t, size_t> nodes() const { return {0u, nodes_count()}; }
};

}

#endif