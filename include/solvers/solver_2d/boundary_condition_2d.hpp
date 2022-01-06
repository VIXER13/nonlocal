#ifndef NONLOCAL_BOUNDARY_CONDITION_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_2D_HPP

#include "../solvers_constants.hpp"
#include <functional>

namespace nonlocal {

template<class B, class T, size_t DoF, class Signature>
struct boundary_condition_2d_t final {
    static_assert(DoF > 0, "DoF must be greater than 0.");

    struct node final {
        B type = B(boundary_condition_t::SECOND_KIND);
        std::function<Signature> func = [](Signature) constexpr noexcept { return T{0}; };
    };

    std::conditional_t<DoF == 1, node, std::array<node, DoF>> cond;

    const node& operator[]([[maybe_unused]] const size_t i) const noexcept {
        if constexpr (DoF == 1)
            return cond;
        else
            return cond[i];
    }
};

template<class B, class T, size_t DoF>
using stationary_boundary_2d_t = boundary_condition_2d_t<B, T, DoF, T(const std::array<T, 2>&)>;

template<class B, class T, size_t DoF>
using nonstationary_boundary_2d_t = boundary_condition_2d_t<B, T, DoF, T(const T, const std::array<T, 2>&)>;

}

#endif