#pragma once

#include "element_2d_base.hpp"
#include "derivative_element_basis_2d.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Element_Type, auto...Args>
class element_2d : public virtual element_2d_base<T>,
                   public derivative_element_basis_2d<T, 2, Element_Type, Args...> {
    using derivative_base = derivative_element_basis_2d<T, 2, Element_Type, Args...>;

public:
    ~element_2d() override = default;

    size_t nodes_count() const override { return derivative_base::N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return derivative_base::nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& x) const override { return derivative_base::N   [i](x); }
    T Nxi (const size_t i, const std::array<T, 2>& x) const override { return derivative_base::Nxi [i](x); }
    T Neta(const size_t i, const std::array<T, 2>& x) const override { return derivative_base::Neta[i](x); }

    T boundary(const side_2d bound, const T x) const override { return derivative_base::boundary(bound, x); }
};

}