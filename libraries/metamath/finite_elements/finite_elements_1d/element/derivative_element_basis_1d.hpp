#ifndef FINITE_ELEMENT_DERIVATIVE_BASIS_1D_HPP
#define FINITE_ELEMENT_DERIVATIVE_BASIS_1D_HPP

#include "derivative.hpp"
#include "to_function.hpp"
#include "simplify.hpp"

namespace metamath::finite_element {

template<class T, size_t Parameters_Count, size_t Derivative_Order, template<class, auto...> class Element_Type, auto... Args>
class derivative_element_basis_1d : public Element_Type<T, Args...> {
    using element_t = Element_Type<T, Args...>;
    using basis_function = std::function<T(const std::array<T, Parameters_Count>&)>;

protected:
    using element_t::x;
    using element_t::basis;
    using element_t::nodes;
    static_assert(std::tuple_size_v<decltype(basis)> == nodes.size(), "The number of functions and nodes does not match.");
    
private:
    using basis_functions = std::array<std::array<basis_function, nodes.size()>, Derivative_Order + 1>;

    template<size_t Order, class... E>
    static auto derivatives(basis_functions& result, const std::tuple<E...>& e) {
        if constexpr (Order <= Derivative_Order) {
            const auto expr = symbolic::simplify(symbolic::derivative<x>(e));
            result[Order] = symbolic::to_function<T, Parameters_Count>(expr);
            derivatives<Order + 1>(result, expr);
        }
    }

    template<class... E>
    static auto derivatives(const std::tuple<E...>& e) {
        basis_functions result;
        const auto expr = symbolic::simplify(e);
        result.front() = symbolic::to_function<T, Parameters_Count>(expr);
        derivatives<1>(result, expr);
        return result;
    }

protected:
    static inline const basis_functions basis_and_derivatives = derivatives(basis);

    explicit derivative_element_basis_1d() = default;
    ~derivative_element_basis_1d() override = default;
};

}

#endif