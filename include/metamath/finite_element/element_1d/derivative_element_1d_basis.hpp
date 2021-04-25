#ifndef DERIVATIVE_FINITE_ELEMENT_1D_BASIS_HPP
#define DERIVATIVE_FINITE_ELEMENT_1D_BASIS_HPP

#include "derivative.hpp"
#include "to_function.hpp"

namespace metamath::finite_element {

template<class T, template<class> class Element_Type, size_t Parameters_Count>
class derivative_element_1d_basis : public Element_Type<T> {
    using Element_Type<T>::xi;
    using Element_Type<T>::basis;
    using Element_Type<T>::nodes;

    static_assert(std::tuple_size<decltype(basis)>::value == nodes.size(), "The number of functions and nodes does not match.");

protected:
    static inline const std::array<std::function<T(const std::array<T, Parameters_Count>&)>, nodes.size()>
        N   = symdiff::to_function<T, Parameters_Count>(basis),
        Nxi = symdiff::to_function<T, Parameters_Count>(symdiff::derivative<xi>(basis));

    explicit derivative_element_1d_basis() = default;
    ~derivative_element_1d_basis() override = default;
};

}

#endif