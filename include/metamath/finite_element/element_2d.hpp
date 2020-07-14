#ifndef FINITE_ELEMENT_2D_HPP
#define FINITE_ELEMENT_2D_HPP

#include "element_base.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class T, template<class> class Element_Type>
class element_2d : public virtual element_2d_base<T>,
                   public Element_Type<T> {
    static_assert(Element_Type<T>::N.size() == Element_Type<T>::nodes.size(),
                  "The number of functions and nodes does not match.");
    static_assert(Element_Type<T>::N.size() == Element_Type<T>::Nxi.size() &&
                  Element_Type<T>::N.size() == Element_Type<T>::Neta.size(),
                  "The number of functions and their derivatives does not match.");

public:
    size_t nodes_count() const override { return Element_Type<T>::N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return Element_Type<T>::nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& xi) const override { return Element_Type<T>::N   [i](xi); }
    T Nxi (const size_t i, const std::array<T, 2>& xi) const override { return Element_Type<T>::Nxi [i](xi); }
    T Neta(const size_t i, const std::array<T, 2>& xi) const override { return Element_Type<T>::Neta[i](xi); }

    T boundary(const side_2d bound, const T x) const override { return Element_Type<T>::boundary(bound, x); }
};

}

#endif