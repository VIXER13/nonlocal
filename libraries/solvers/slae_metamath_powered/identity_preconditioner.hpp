#pragma once

#include "preconditioner_base.hpp"

namespace nonlocal::slae {

template<class T>
struct identity_preconditioner final : public preconditioner_base<T> {
    using typename preconditioner_base<T>::entity_t;

    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        return rhs;
    }
};

}