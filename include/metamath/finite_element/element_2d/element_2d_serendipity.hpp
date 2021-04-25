#ifndef FINITE_ELEMENT_2D_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_SERENDIPITY_HPP

#include "element_2d.hpp"
#include "basis/quadratic_serendipity.hpp"
#include "basis/qubic_serendipity.hpp"

namespace metamath::finite_element {

template<class T, template<class> class Element_Type>
class element_2d_serendipity : public virtual element_2d_base<T>,
                               public derivative_element_2d_basis<T, Element_Type, 3> {
    using derivative_base = derivative_element_2d_basis<T, Element_Type, 3>;
    using Element_Type<T>::_p;

public:
    ~element_2d_serendipity() override = default;

    size_t nodes_count() const override { return derivative_base::N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return Element_Type<T>::nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::N   [i]({xi[0], xi[1], _p}); }
    T Nxi (const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::Nxi [i]({xi[0], xi[1], _p}); }
    T Neta(const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::Neta[i]({xi[0], xi[1], _p}); }

    T boundary(const side_2d bound, const T x) const override { return Element_Type<T>::boundary(bound, x); }
};

// Специализация под квадратичные серендиповы элементы
template<class T>
class element_2d<T, quadratic_serendipity> : public element_2d_serendipity<T, quadratic_serendipity> {
public:
    ~element_2d() override = default;
};

// Специализация под кубические серендиповы элементы
template<class T>
class element_2d<T, qubic_serendipity> : public element_2d_serendipity<T, qubic_serendipity> {
public:
    ~element_2d() override = default;
};

}

#endif