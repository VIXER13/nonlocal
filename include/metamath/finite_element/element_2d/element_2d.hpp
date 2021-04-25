#ifndef FINITE_ELEMENT_2D_HPP
#define FINITE_ELEMENT_2D_HPP

#include "element_2d_base.hpp"
#include "derivative_element_2d_basis.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type.
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class T, template<class> class Element_Type>
class element_2d : public virtual element_2d_base<T>,
                   public derivative_element_2d_basis<T, Element_Type, 2> {
    using derivative_base = derivative_element_2d_basis<T, Element_Type, 2>;

public:
    ~element_2d() override = default;

    size_t nodes_count() const override { return derivative_base::N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return Element_Type<T>::nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::N   [i](xi); }
    T Nxi (const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::Nxi [i](xi); }
    T Neta(const size_t i, const std::array<T, 2>& xi) const override { return derivative_base::Neta[i](xi); }

    T boundary(const side_2d bound, const T x) const override { return Element_Type<T>::boundary(bound, x); }
};

}

#endif