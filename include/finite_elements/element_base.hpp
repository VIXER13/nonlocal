#ifndef FINITE_ELEMENT_BASE_HPP
#define FINITE_ELEMENT_BASE_HPP

// В данном модуле описаны базовые интерфейсы классов конечных элементов.

#include "matrix.hpp"
#include "element_1d.hpp"
#include "element_2d.hpp"
#include "quadrature.hpp"

namespace finite_element {

class element_base {
public:
    virtual size_t nodes_count() const = 0; // Любой конечный элемент, вне зависимости от его размерности, имеет некоторое количество узлов.

    virtual ~element_base() noexcept = default;
};

template<class Type>
class element_1d_base : public element_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

public:
    virtual Type node(const size_t i) const = 0;

    virtual Type N  (const size_t i, const Type xi) const = 0; // Обращение к i-ой функции формы в точке xi.
    virtual Type Nxi(const size_t i, const Type xi) const = 0; // Аналогично для производной.

    virtual Type boundary(const side_1d bound) const = 0; // Геометрия элемента.
};

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса ElementType. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class Type, template<class> class Element_Type>
class element_1d : public virtual element_1d_base<Type>,
                   public Element_Type<Type> {
    static_assert(Element_Type<Type>::basic_N.size() == Element_Type<Type>::basic_Nxi.size() &&
                  Element_Type<Type>::basic_N.size() == Element_Type<Type>::nodes.size(),
                  "The number of functions and their derivatives does not match.");
public:
    size_t nodes_count() const override { return Element_Type<Type>::basic_N.size(); }

    Type node(const size_t i) const override { return Element_Type<Type>::nodes[i]; }
        
    Type N  (const size_t i, const Type xi) const override { return Element_Type<Type>::basic_N  [i](xi); }
    Type Nxi(const size_t i, const Type xi) const override { return Element_Type<Type>::basic_Nxi[i](xi); }

    Type boundary(const side_1d bound) const override { return Element_Type<Type>::boundary(bound); }
};

template<class Type>
class element_2d_base : public element_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

public:
    virtual const std::array<Type, 2>& node(const size_t i) const = 0;

    virtual Type N   (const size_t i, const Type xi, const Type eta) const = 0; // Обращение к i-ой функции формы в точке (xi, eta).
    virtual Type Nxi (const size_t i, const Type xi, const Type eta) const = 0; // Аналогично для производных.
    virtual Type Neta(const size_t i, const Type xi, const Type eta) const = 0;

    virtual Type boundary(const side_2d bound, const Type x) const = 0; // Геометрия элемента.
};

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class Type, template<class> class Element_Type>
class element_2d : public virtual element_2d_base<Type>,
                   public Element_Type<Type> {
    static_assert(Element_Type<Type>::basic_N.size() == Element_Type<Type>::basic_Nxi.size() &&
                  Element_Type<Type>::basic_N.size() == Element_Type<Type>::basic_Neta.size() &&
                  Element_Type<Type>::basic_N.size() == Element_Type<Type>::nodes.size(),
                  "The number of functions and their derivatives does not match.");

public:
    size_t nodes_count() const override { return Element_Type<Type>::basic_N.size(); }

    const std::array<Type, 2>& node(const size_t i) const override { return Element_Type<Type>::nodes[i]; }

    Type N   (const size_t i, const Type xi, const Type eta) const override { return Element_Type<Type>::basic_N   [i](xi, eta); }
    Type Nxi (const size_t i, const Type xi, const Type eta) const override { return Element_Type<Type>::basic_Nxi [i](xi, eta); }
    Type Neta(const size_t i, const Type xi, const Type eta) const override { return Element_Type<Type>::basic_Neta[i](xi, eta); }

    Type boundary(const side_2d bound, const Type x) const override { return Element_Type<Type>::boundary(bound, x); }
};

}

#endif