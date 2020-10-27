#ifndef FINITE_ELEMENT_2D_HPP
#define FINITE_ELEMENT_2D_HPP

#include "element_base.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class T, template<class> class Element_Type>
class element_2d : public virtual element_2d_base<T>,
                   public Element_Type<T> {
    using Element_Type<T>::xi;
    using Element_Type<T>::eta;
    using Element_Type<T>::nodes;
    using Element_Type<T>::basis;

    static_assert(std::tuple_size<decltype(basis)>::value == nodes.size(), "The number of functions and nodes does not match.");

protected:
    static inline const std::array<std::function<T(const std::array<T, 2>&)>, nodes.size()>
        _N    = symdiff::to_function<T, 2>(basis),
        _Nxi  = symdiff::to_function<T, 2>(symdiff::derivative<xi>(basis)),
        _Neta = symdiff::to_function<T, 2>(symdiff::derivative<eta>(basis));

public:
    ~element_2d() override = default;

    size_t nodes_count() const override { return _N.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return nodes[i]; }

    T N   (const size_t i, const std::array<T, 2>& xi) const override { return _N   [i](xi); }
    T Nxi (const size_t i, const std::array<T, 2>& xi) const override { return _Nxi [i](xi); }
    T Neta(const size_t i, const std::array<T, 2>& xi) const override { return _Neta[i](xi); }

    T boundary(const side_2d bound, const T x) const override { return Element_Type<T>::boundary(bound, x); }
};

}

#endif