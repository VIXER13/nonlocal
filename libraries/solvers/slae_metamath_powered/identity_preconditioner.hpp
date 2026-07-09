#pragma once

#include "preconditioner_base.hpp"

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
class identity_preconditioner final : public preconditioner_base<T, I, J> {
public:
    using typename preconditioner_base<T, I, J>::entity_t;

    void compute(metamath::linear::sparse_matrix<T, I, J>&&) override {}

    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        return rhs;
    }
};

}