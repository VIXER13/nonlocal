#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

#include "element_base.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class T, template<class> class Element_Type>
class element_1d : public virtual element_1d_base<T>,
                   public Element_Type<T> {
    static_assert(Element_Type<T>::N.size() == Element_Type<T>::nodes.size(),
                  "The number of functions and nodes does not match.");
    static_assert(Element_Type<T>::N.size() == Element_Type<T>::Nxi.size(),
                  "The number of functions and their derivatives does not match.");
                  
public:
    size_t nodes_count() const override { return Element_Type<T>::N.size(); }

    T node(const size_t i) const override { return Element_Type<T>::nodes[i]; }
        
    T N  (const size_t i, const T xi) const override { return Element_Type<T>::N  [i]({xi}); }
    T Nxi(const size_t i, const T xi) const override { return Element_Type<T>::Nxi[i]({xi}); }

    T boundary(const side_1d bound) const override { return Element_Type<T>::boundary(bound); }
};

}

#endif