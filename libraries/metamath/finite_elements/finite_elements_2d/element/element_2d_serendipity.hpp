#pragma once

#include "element_2d.hpp"
#include "derivative_element_basis_2d.hpp"
#include "basis/serendipity_2d.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Element_Type, auto...Args>
class element_2d_serendipity : public virtual element_2d_base<T>,
                               public derivative_element_basis_2d<T, 3, Element_Type, Args...> {
    using derivative_base = derivative_element_basis_2d<T, 3, Element_Type, Args...>;

protected:
    using derivative_base::_p;

public:
    ~element_2d_serendipity() override = default;

    size_t nodes_count() const override { return derivative_base::N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return derivative_base::nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::N   [i]({xi[0], xi[1], _p}); }
    T Nxi (const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::Nxi [i]({xi[0], xi[1], _p}); }
    T Neta(const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::Neta[i]({xi[0], xi[1], _p}); }

    T boundary(const side_2d bound, const T x) const override { return derivative_base::boundary(bound, x); }
};

template<class T>
class element_2d<T, serendipity, 2> : public element_2d_serendipity<T, serendipity, 2> {
public:
    ~element_2d() override = default;
};

template<class T>
class element_2d<T, serendipity, 3> : public element_2d_serendipity<T, serendipity, 3> {
public:
    ~element_2d() override = default;
};

}