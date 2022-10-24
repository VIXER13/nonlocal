#ifndef NONLOCAL_EQUATION_PARAMETERS_HPP
#define NONLOCAL_EQUATION_PARAMETERS_HPP

#include "nonlocal_constants.hpp"

#include <array>
#include <functional>

namespace nonlocal {

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
class equation_parameters final {
    static_assert(Dimension > 0, "Dimension must be non-zero.");
    Physical<T, Args...> physical;
    struct final {
        using arg = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;
        std::function<T(const arg&, const arg&)> influence = [](const arg&, const arg&) constexpr noexcept { return T{0}; };
        T local_weight = T{1};
    } model;
};

}

#endif