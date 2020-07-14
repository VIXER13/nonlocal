#ifndef FINITE_ELEMENT_2D_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_SERENDIPITY_HPP

// В силу особенностей серендиповой аппроксимации, были выделены две специализации класса двумерных конечных элементов.

#include "element_base.hpp"
#include "element_2d_strategy/quadratic_serendipity.hpp"
#include "element_2d_strategy/qubic_serendipity.hpp"

namespace metamath::finite_element {

template<class T, template<class> class Element_Type>
class element_2d_serendipity : public virtual element_2d_base<T>,
                               public Element_Type<T> {
    static_assert(Element_Type<T>::N.size() == Element_Type<T>::nodes.size(),
                  "The number of functions and nodes does not match.");
    static_assert(Element_Type<T>::N.size() == Element_Type<T>::Nxi.size() &&
                  Element_Type<T>::N.size() == Element_Type<T>::Neta.size(),
                  "The number of functions and their derivatives does not match.");

public:
    size_t nodes_count() const override { return Element_Type<T>::N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return Element_Type<T>::nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& xi) const override { return Element_Type<T>::N   [i]({xi[0], xi[1], Element_Type<T>::_p}); }
    T Nxi (const size_t i, const std::array<T, 2>& xi) const override { return Element_Type<T>::Nxi [i]({xi[0], xi[1], Element_Type<T>::_p}); }
    T Neta(const size_t i, const std::array<T, 2>& xi) const override { return Element_Type<T>::Neta[i]({xi[0], xi[1], Element_Type<T>::_p}); }

    T boundary(const side_2d bound, const T x) const override { return Element_Type<T>::boundary(bound, x); }
};

// Специализация под квадратичные серендиповы элементы
template<class T>
class element_2d<T, quadratic_serendipity> : public element_2d_serendipity<T, quadratic_serendipity> {};

// Специализация под кубические серендиповы элементы
template<class T>
class element_2d<T, qubic_serendipity> : public element_2d_serendipity<T, qubic_serendipity> {};

}

#endif