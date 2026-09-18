#pragma once

#include <metamath/types/traits.hpp>

namespace nonlocal::slae {

template<class T>
struct preconditioner_base {
    using entity_t = metamath::types::container_type_t<T>;
    using floating_point_t = metamath::types::container_type_t<entity_t>;

    virtual ~preconditioner_base() noexcept = default;
    virtual std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const = 0;
};

}