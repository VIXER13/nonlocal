#pragma once

#include <metamath/linear/sparse_matrix.hpp>
#include <metamath/types/traits.hpp>

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
struct preconditioner_base {
    using entity_t = metamath::types::container_type_t<T>;

    virtual ~preconditioner_base() noexcept = default;
    virtual void compute(metamath::linear::sparse_matrix<T, I, J>&& matrix) = 0;
    virtual std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const = 0;
};

}