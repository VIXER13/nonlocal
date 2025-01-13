#pragma once

#include "derivative.hpp"
#include "to_function.hpp"
#include "simplify.hpp"

namespace metamath::finite_element {

template<class T, size_t Parameters_Count, template<class, auto...> class Element_Type, auto... Args>
class derivative_element_basis_1d : public Element_Type<T, Args...> {
    using element_t = Element_Type<T, Args...>;

protected:
    using element_t::x;
    using element_t::basis;
    using element_t::nodes;

    static_assert(std::tuple_size_v<decltype(basis)> == nodes.size(), "The number of functions and nodes does not match.");

    static inline const std::array<std::function<T(const std::array<T, Parameters_Count>&)>, nodes.size()>
        basis_as_functions   = symbolic::to_function<T, Parameters_Count>(symbolic::simplify(basis)),
        differentiated_basis = symbolic::to_function<T, Parameters_Count>(symbolic::simplify(symbolic::derivative<x>(basis)));

    explicit derivative_element_basis_1d() = default;
    ~derivative_element_basis_1d() override = default;
};

}