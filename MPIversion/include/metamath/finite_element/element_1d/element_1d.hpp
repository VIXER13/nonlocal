#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

#include "element_1d_base.hpp"
#include "derivative.hpp"
#include "to_function.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type.
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class T, template<class> class Element_Type>
class element_1d : public virtual element_1d_base<T>,
                   public Element_Type<T> {
    using Element_Type<T>::xi;
    using Element_Type<T>::basis;
    using Element_Type<T>::nodes;

    static_assert(std::tuple_size<decltype(basis)>::value == nodes.size(), "The number of functions and nodes does not match.");

protected:
    static inline const std::array<std::function<T(const std::array<T, 1>&)>, nodes.size()>
        _N   = symdiff::to_function<T, 1>(basis),
        _Nxi = symdiff::to_function<T, 1>(symdiff::derivative<xi>(basis));

public:
    ~element_1d() override = default;

    size_t nodes_count() const override { return _N.size(); }

    T node(const size_t i) const override { return nodes[i]; }

    T N  (const size_t i, const T xi) const override { return _N  [i]({xi}); }
    T Nxi(const size_t i, const T xi) const override { return _Nxi[i]({xi}); }

    T boundary(const side_1d bound) const override { return Element_Type<T>::boundary(bound); }
};

}

#endif