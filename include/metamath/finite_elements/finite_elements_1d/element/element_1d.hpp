#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

#include "element_1d_base.hpp"
#include "symbolic/symbolic.hpp"

#include <functional>
#include <tuple>

namespace metamath::finite_element {

template<class T, template<class, auto...> class Element_Type, auto... Args>
class element_1d : public virtual element_1d_base<T>,
                   public Element_Type<T, Args...> {
    using element_t = Element_Type<T, Args...>;
    using element_t::nodes;
    using element_t::basis;
    using element_t::x;

    static_assert(std::tuple_size_v<decltype(basis)> == nodes.size(), "The number of functions and nodes does not match.");

protected:
    static inline const std::array<std::function<T(const std::array<T, 1>&)>, nodes.size()>
        _N   = symbolic::to_function<T, 1>(symbolic::simplify(basis)),
        _Nxi = symbolic::to_function<T, 1>(symbolic::simplify(symbolic::derivative<x>(basis)));

public:
    ~element_1d() override = default;

    size_t nodes_count() const override { return _N.size(); }

    T node(const size_t i) const override { return nodes[i]; }

    T N  (const size_t i, const T xi) const override { return _N  [i]({xi}); }
    T Nxi(const size_t i, const T xi) const override { return _Nxi[i]({xi}); }

    T boundary(const side_1d bound) const override { return element_t::boundary(bound); }
};

}

#endif