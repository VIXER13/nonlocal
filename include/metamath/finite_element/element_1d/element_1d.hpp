#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

#include "element_1d_base.hpp"
#include "derivative_element_1d_basis.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type.
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class T, template<class> class Element_Type>
class element_1d : public virtual element_1d_base<T>,
                   public derivative_element_1d_basis<T, Element_Type, 1> {
    using derivative_base = derivative_element_1d_basis<T, Element_Type, 1>;

public:
    ~element_1d() override = default;

    size_t nodes_count() const override { return derivative_base::N.size(); }

    T node(const size_t i) const override { return Element_Type<T>::nodes[i]; }

    T N  (const size_t i, const T xi) const override { return derivative_base::N  [i]({xi}); }
    T Nxi(const size_t i, const T xi) const override { return derivative_base::Nxi[i]({xi}); }

    T boundary(const side_1d bound) const override { return Element_Type<T>::boundary(bound); }
};

}

#endif